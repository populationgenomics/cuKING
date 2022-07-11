#include <absl/container/flat_hash_map.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <absl/status/status.h>
#include <absl/status/statusor.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/strip.h>
#include <absl/synchronization/blocking_counter.h>
#include <absl/time/time.h>
#include <arrow/filesystem/filesystem.h>
#include <arrow/result.h>
#include <arrow/status.h>
#include <google/cloud/storage/client.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <queue>
#include <string>
#include <thread>
#include <vector>

ABSL_FLAG(std::string, input_uri, "",
          "GCS URI containing the packed Parquet tables, e.g. "
          "gs://some/bucket/my_table.parquet");
ABSL_FLAG(std::string, output_uri, "",
          "The sparse relatedness matrix JSON output GCS URI, e.g. "
          "gs://some/bucket/relatedness.json");
ABSL_FLAG(size_t, num_reader_threads, 100,
          "How many threads to use for processing of Parquet partitions. This "
          "influences the amount of memory required.");
ABSL_FLAG(
    uint32_t, max_results, 100 << 20,
    "How many coefficients for related sample pairs to reserve memory for.");
ABSL_FLAG(
    float, king_coeff_threshold, 0.0442f,
    "Only store coefficients larger than this threshold. Defaults to 3rd "
    "degree or closer (see https://www.kingrelatedness.com/manual.shtml).");

namespace {

namespace gcs = google::cloud::storage;

// Returns ceil(a / b) for integers a, b.
template <typename T>
inline T CeilIntDiv(const T a, const T b) {
  return (a + b - 1) / b;
}

// Keeps track of time intervals.
class StopWatch {
 public:
  absl::Duration GetElapsedAndReset() {
    const auto now = absl::Now();
    const auto result = now - last_time_;
    last_time_ = now;
    return result;
  }

 private:
  absl::Time last_time_ = absl::Now();
};

inline absl::Status ToAbslStatus(absl::Status status) {
  return std::move(status);
}

inline absl::Status ToAbslStatus(const arrow::Status &status) {
  return absl::UnknownError(status.ToString());
}

// Unfortunately these status macros are still not part of Abseil. This is
// adapted from
// https://source.chromium.org/chromiumos/chromiumos/codesearch/+/main:src/platform2/missive/util/status_macros.h
// but allows converting from different Status implementations
// (e.g. absl::Status and arrow::Status) by providing an override for the
// ToAbslStatus converter.
#define RETURN_IF_ERROR(expr)                                                \
  do {                                                                       \
    /* Using _status below to avoid capture problems if expr is "status". */ \
    const auto _status = (expr);                                             \
    if (ABSL_PREDICT_FALSE(!_status.ok())) {                                 \
      return ToAbslStatus(_status);                                          \
    }                                                                        \
  } while (0)

// Internal helper for concatenating macro values.
#define STATUS_MACROS_CONCAT_NAME_INNER(x, y) x##y
#define STATUS_MACROS_CONCAT_NAME(x, y) STATUS_MACROS_CONCAT_NAME_INNER(x, y)

#define ASSIGN_OR_RETURN_IMPL(result, lhs, rexpr) \
  auto result = rexpr;                            \
  if (ABSL_PREDICT_FALSE(!result.ok())) {         \
    return ToAbslStatus(result.status());         \
  }                                               \
  lhs = *std::move(result);

#define ASSIGN_OR_RETURN(lhs, rexpr) \
  ASSIGN_OR_RETURN_IMPL(             \
      STATUS_MACROS_CONCAT_NAME(_status_or_value, __COUNTER__), lhs, rexpr)

// Custom deleter for RAII-style CUDA-managed array.
template <typename T>
struct CudaArrayDeleter {
  void operator()(T *const val) const { cudaFree(val); }
};

template <typename T>
using CudaArray = std::unique_ptr<T[], CudaArrayDeleter<T>>;

template <typename T>
CudaArray<T> NewCudaArray(const size_t size) {
  static_assert(std::is_pod<T>::value, "A must be a POD type.");
  T *buffer = nullptr;
  const auto err = cudaMallocManaged(&buffer, size * sizeof(T));
  if (err) {
    std::cerr << "Error: can't allocate CUDA memory: "
              << cudaGetErrorString(err) << std::endl;
    exit(1);
  }
  return CudaArray<T>(buffer, CudaArrayDeleter<T>());
}

__device__ float ComputeKing(const uint64_t *const het_i_entries,
                             const uint64_t *const hom_alt_i_entries,
                             const uint64_t *const het_j_entries,
                             const uint64_t *const hom_alt_j_entries,
                             const uint32_t num_entries) {
  // See https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.king.
  uint32_t num_het_i = 0, num_het_j = 0, num_both_het = 0, num_opposing_hom = 0;
  for (uint32_t k = 0; k < num_entries; ++k) {
    const uint64_t het_i = het_i_entries[k];
    const uint64_t hom_alt_i = hom_alt_i_entries[k];
    const uint64_t hom_ref_i = (~het_i) & (~hom_alt_i);

    const uint64_t het_j = het_j_entries[k];
    const uint64_t hom_alt_j = hom_alt_j_entries[k];
    const uint64_t hom_ref_j = (~het_j) & (~hom_alt_j);

    // Only count sites where both genotypes are defined.
    const uint64_t missing_mask = ~(het_i & hom_alt_i) & ~(het_j & hom_alt_j);

    num_het_i += __popcll(het_i & missing_mask);
    num_het_j += __popcll(het_j & missing_mask);
    num_both_het += __popcll(het_i & het_j & missing_mask);
    num_opposing_hom += __popcll(
        ((hom_ref_i & hom_alt_j) | (hom_ref_j & hom_alt_i)) & missing_mask);
  }

  // Return the "between-family" estimator.
  const uint32_t min_hets = num_het_i < num_het_j ? num_het_i : num_het_j;
  return 0.5f +
         (2.f * num_both_het - 4.f * num_opposing_hom - num_het_i - num_het_j) /
             (4.f * min_hets);
}

// Stores the KING coefficient for one pair of samples.
struct KingResult {
  uint32_t sample_i, sample_j;
  float coeff;
};

__global__ void ComputeKingKernel(const uint32_t num_samples,
                                  const uint32_t words_per_sample,
                                  const uint64_t *const bit_sets,
                                  const float coeff_threshold,
                                  const uint32_t max_results,
                                  KingResult *const results,
                                  uint32_t *const result_index) {
  const uint64_t index = uint64_t(blockIdx.x) * blockDim.x + threadIdx.x;
  const uint64_t i = index / num_samples;
  const uint64_t j = index % num_samples;
  if (i >= num_samples || i >= j) {
    return;
  }

  const uint32_t num_entries = words_per_sample / 2;
  const uint64_t offset_i = i * words_per_sample;
  const uint64_t offset_j = j * words_per_sample;
  const float coeff = ComputeKing(
      bit_sets + offset_i, bit_sets + offset_i + num_entries,
      bit_sets + offset_j, bit_sets + offset_j + num_entries, num_entries);

  if (coeff > coeff_threshold) {
    // Reserve a result slot atomically to avoid collisions.
    const uint32_t reserved = atomicAdd(result_index, 1u);
    if (reserved < max_results) {
      KingResult &result = results[reserved];
      result.sample_i = i;
      result.sample_j = j;
      result.coeff = coeff;
    }
  }
}

// Returns the bucket + blob name from a full gs:// URI.
absl::StatusOr<std::pair<std::string_view, std::string_view>> SplitGcsUri(
    std::string_view uri) {
  if (!absl::ConsumePrefix(&uri, "gs://")) {
    return absl::InvalidArgumentError(absl::StrCat("Unsupported URI: ", uri));
  }

  const size_t slash_pos = uri.find_first_of('/');
  if (slash_pos == std::string_view::npos) {
    return absl::InvalidArgumentError(
        absl::StrCat("Incomplete blob URI ", uri));
  }

  return std::make_pair(uri.substr(0, slash_pos), uri.substr(slash_pos + 1));
}

// Atomically clears a bit in a bit set.
inline void AtomicClearBit(uint64_t *const bit_set, uint64_t index) {
  uint64_t *const ptr = bit_set + (index >> 6);
  // C++20 adds std::atomic_ref, but we're compiling with C++17. We can use
  // relaxed memory ordering as bit set values don't depend on one another.
  __atomic_and_fetch(ptr, ~(uint64_t(1) << (index & size_t(63))),
                     __ATOMIC_RELAXED);
}

// Adapted from the Abseil thread pool.
class ThreadPool {
 public:
  explicit ThreadPool(const size_t num_threads) {
    for (size_t i = 0; i < num_threads; ++i) {
      threads_.push_back(std::thread(&ThreadPool::WorkLoop, this));
    }
  }

  ThreadPool(const ThreadPool &) = delete;
  ThreadPool &operator=(const ThreadPool &) = delete;

  ~ThreadPool() {
    {
      absl::MutexLock l(&mu_);
      for (size_t i = 0; i < threads_.size(); i++) {
        queue_.push(nullptr);  // Shutdown signal.
      }
    }
    for (auto &t : threads_) {
      t.join();
    }
  }

  // Schedule a function to be run on a ThreadPool thread immediately.
  void Schedule(std::function<void()> func) {
    assert(func != nullptr);
    absl::MutexLock l(&mu_);
    queue_.push(std::move(func));
  }

 private:
  bool WorkAvailable() const ABSL_EXCLUSIVE_LOCKS_REQUIRED(mu_) {
    return !queue_.empty();
  }

  void WorkLoop() {
    while (true) {
      std::function<void()> func;
      {
        absl::MutexLock l(&mu_);
        mu_.Await(absl::Condition(this, &ThreadPool::WorkAvailable));
        func = std::move(queue_.front());
        queue_.pop();
      }
      if (func == nullptr) {  // Shutdown signal.
        break;
      }
      func();
    }
  }

  absl::Mutex mu_;
  std::queue<std::function<void()>> queue_ ABSL_GUARDED_BY(mu_);
  std::vector<std::thread> threads_;
};

// Invokes the given function in parallel over a range. If any invocations
// returns an error, only one such error will be returned
// (non-deterministically).
absl::Status ParallelFor(ThreadPool *const thread_pool, const size_t begin,
                         const size_t end,
                         std::function<absl::Status(size_t index)> func) {
  absl::BlockingCounter blocking_counter(end - begin);
  absl::Mutex mu;
  absl::Status result ABSL_GUARDED_BY(mu);
  for (size_t i = begin; i < end; ++i) {
    thread_pool->Schedule([&, i] {
      auto status = func(i);
      if (!status.ok()) {
        const absl::MutexLock lock(&mu);
        result.Update(std::move(status));
      }
      blocking_counter.DecrementCount();
    });
  }
  blocking_counter.Wait();
  return result;
}

absl::Status Run() {
  // Validate flags.
  const auto input_uri = absl::GetFlag(FLAGS_input_uri);
  if (input_uri.empty()) {
    return absl::InvalidArgumentError("No input URI specified");
  }
  const auto output_uri = absl::GetFlag(FLAGS_output_uri);
  if (output_uri.empty()) {
    return absl::InvalidArgumentError("No output URI specified");
  }
  const auto num_reader_threads = absl::GetFlag(FLAGS_num_reader_threads);
  if (num_reader_threads == 0) {
    return absl::InvalidArgumentError("Invalid number of reader threads");
  }

  auto gcs_client =
      gcs::Client(google::cloud::Options{}.set<gcs::ConnectionPoolSizeOption>(
          num_reader_threads));

  // Find all Parquet table files.
  StopWatch stop_watch;
  std::cout << "Listing input files...";
  std::cout.flush();
  ASSIGN_OR_RETURN(auto bucket_and_path, SplitGcsUri(input_uri));

  std::cout << "Reading metadata...";
  std::cout.flush();
  std::string input_path;
  ASSIGN_OR_RETURN(const auto input_fs,
                   arrow::fs::FileSystemFromUri(input_uri, &input_path));
  ASSIGN_OR_RETURN(auto metadata_file, input_fs->OpenInputFile(absl::StrCat(
                                           input_path, "/metadata.json")));
  ASSIGN_OR_RETURN(const uint64_t metadata_file_size, metadata_file->GetSize());
  std::vector<uint8_t> metadata_buffer(metadata_file_size);
  ASSIGN_OR_RETURN(
      const uint64_t metadata_bytes_read,
      metadata_file->ReadAt(0, metadata_file_size, metadata_buffer.data()));
  if (metadata_bytes_read != metadata_file_size) {
    return absl::FailedPreconditionError(absl::StrCat(
        "Expected to read ", metadata_file_size, " metadata bytes, but got ",
        metadata_bytes_read, " bytes instead"));
  }
  const auto metadata = nlohmann::json::parse(
      metadata_buffer.begin(), metadata_buffer.end(),
      /* parser_callback_t */ nullptr, /* allow_exceptions */ false);
  if (metadata.is_discarded()) {
    return absl::FailedPreconditionError("Failed to parse metadata JSON");
  }
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  const std::vector<std::string_view> sample_ids = metadata["samples"];
  const uint32_t num_samples = sample_ids.size();
  std::cout << "Metadata contains " << num_samples << " samples with "
            << metadata["num_sites"] << " sites." << std::endl;

  std::cout << "Reading bit sets...";
  std::cout.flush();
  ASSIGN_OR_RETURN(auto bit_set_file, input_fs->OpenInputFile(absl::StrCat(
                                          input_path, "/bit_set.bin")));
  ASSIGN_OR_RETURN(const uint64_t bit_set_file_size, bit_set_file->GetSize());
  const uint32_t words_per_sample = metadata["words_per_sample"];
  if (bit_set_file_size !=
      uint64_t(num_samples) * words_per_sample * sizeof(uint64_t)) {
    return absl::FailedPreconditionError(
        absl::StrCat("Unexpected bit set file size: ", bit_set_file_size));
  }
  auto bit_set =
      NewCudaArray<uint64_t>(uint64_t(num_samples) * words_per_sample);
  ASSIGN_OR_RETURN(const uint64_t bit_set_bytes_read,
                   bit_set_file->ReadAt(0, bit_set_file_size, bit_set.get()));
  if (bit_set_bytes_read != bit_set_file_size) {
    return absl::FailedPreconditionError(absl::StrCat(
        "Expected to read ", bit_set_file_size, " bit set bytes, but got ",
        bit_set_bytes_read, " bytes instead"));
  }
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  std::cout << "Allocating memory for results...";
  std::cout.flush();
  const uint32_t kMaxResults = absl::GetFlag(FLAGS_max_results);
  auto results = NewCudaArray<KingResult>(kMaxResults);
  memset(results.get(), 0, sizeof(KingResult) * kMaxResults);
  // We just need a single value.
  auto result_index = NewCudaArray<uint32_t>(1);
  result_index[0] = 0;
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  std::cout << "Running KING CUDA kernel...";
  std::cout.flush();
  constexpr uint64_t kCudaBlockSize = 1024;
  const uint64_t kNumCudaBlocks =
      CeilIntDiv(uint64_t(num_samples) * num_samples, kCudaBlockSize);
  ComputeKingKernel<<<kNumCudaBlocks, kCudaBlockSize>>>(
      num_samples, words_per_sample, bit_set.get(),
      absl::GetFlag(FLAGS_king_coeff_threshold), kMaxResults, results.get(),
      result_index.get());

  // Wait for GPU to finish before accessing on host.
  cudaDeviceSynchronize();
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  // Free some memory for postprocessing.
  bit_set.reset();

  const uint32_t num_results = result_index[0];
  if (num_results > kMaxResults) {
    return absl::ResourceExhaustedError(
        "Could not store all results: try increasing the --max_results "
        "parameter.");
  }

  std::cout << "Found " << num_results
            << " coefficients above the cut-off threshold." << std::endl;

  std::cout << "Processing results...";
  std::cout.flush();
  std::vector<bool> related(num_samples);
  for (uint32_t i = 0; i < num_results; ++i) {
    const auto &result = results[i];
    related[result.sample_i] = related[result.sample_j] = true;
  }

  uint32_t num_related = 0;
  for (size_t i = 0; i < num_samples; ++i) {
    if (related[i]) {
      ++num_related;
    }
  }

  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;
  std::cout << num_related << " related samples found." << std::endl;

  std::cout << "Writing output...";
  std::cout.flush();

  // Create a map for JSON serialization.
  absl::flat_hash_map<std::string_view,
                      absl::flat_hash_map<std::string_view, float>>
      result_map;
  for (size_t i = 0; i < num_results; ++i) {
    const auto &result = results[i];
    result_map[sample_ids[result.sample_i]][sample_ids[result.sample_j]] =
        result.coeff;
  }

  std::string output_path;
  ASSIGN_OR_RETURN(const auto output_fs,
                   arrow::fs::FileSystemFromUri(output_uri, &output_path));
  ASSIGN_OR_RETURN(auto output_stream,
                   output_fs->OpenOutputStream(output_path));
  RETURN_IF_ERROR(output_stream->Write(nlohmann::json(result_map).dump(4)));
  RETURN_IF_ERROR(output_stream->Close());
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  return absl::OkStatus();
}

}  // namespace

int main(int argc, char **argv) {
  absl::ParseCommandLine(argc, argv);

  if (const auto status = Run(); !status.ok()) {
    std::cerr << std::endl << "Error: " << status << std::endl;
    return 1;
  }

  return 0;
}