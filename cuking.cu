#include <absl/container/flat_hash_map.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <absl/status/status.h>
#include <absl/status/statusor.h>
#include <absl/strings/match.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/strip.h>
#include <absl/synchronization/blocking_counter.h>
#include <absl/time/time.h>
#include <arrow/io/memory.h>
#include <google/cloud/storage/client.h>
#include <parquet/column_reader.h>
#include <parquet/file_reader.h>
#include <parquet/metadata.h>

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
ABSL_FLAG(
    std::string, requester_pays_project, "",
    "The user project to use for accessing Requester Pays buckets on GCS.");
ABSL_FLAG(size_t, num_reader_threads, 100,
          "How many threads to use for processing of Parquet partitions. This "
          "influences the amount of memory required.");
ABSL_FLAG(
    uint32_t, max_results, uint32_t(100) << 20,
    "How many coefficients for related sample pairs to reserve memory for.");
ABSL_FLAG(
    float, king_coeff_threshold, 0.0884f,
    "Only store coefficients larger than this threshold. Defaults to 2nd "
    "degree or closer (see https://www.kingrelatedness.com/manual.shtml).");

namespace {

namespace gcs = google::cloud::storage;

inline absl::Status ToAbslStatus(absl::Status status) {
  return std::move(status);
}

inline absl::Status ToAbslStatus(const google::cloud::Status &status) {
  return absl::UnknownError(status.message());
}

// Unfortunately these status macros are still not part of Abseil. This is
// adapted from
// https://source.chromium.org/chromiumos/chromiumos/codesearch/+/main:src/platform2/missive/util/status_macros.h
// but allows converting from different Status implementations (e.g.
// absl::Status and google::cloud::Status) by providing an override for the
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

__global__ void ComputeKingKernel(
    const uint32_t num_samples, const uint32_t words_per_sample,
    const uint64_t *const bit_sets, const float coeff_threshold,
    const uint32_t max_results, KingResult *const results,
    uint32_t *const result_index, uint32_t *const result_overflow) {
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
    } else {
      // result_index might overflow 32 bits, therefore keep a dedicated flag
      // that we ran out of space.
      atomicMax(result_overflow, 1);
    }
  }
}

// Atomically clears a bit in a bit set.
inline void AtomicClearBit(uint64_t *const bit_set, uint64_t index) {
  uint64_t *const ptr = bit_set + (index >> 6);
  // C++20 adds std::atomic_ref, but we're compiling with C++17. We can use
  // relaxed memory ordering as bit set values don't depend on one another.
  __atomic_and_fetch(ptr, ~(uint64_t(1) << (index & size_t(63))),
                     __ATOMIC_RELAXED);
}

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
  ASSIGN_OR_RETURN(const auto input_bucket_and_path, SplitGcsUri(input_uri));

  const auto output_uri = absl::GetFlag(FLAGS_output_uri);
  if (output_uri.empty()) {
    return absl::InvalidArgumentError("No output URI specified");
  }
  ASSIGN_OR_RETURN(const auto output_bucket_and_path, SplitGcsUri(output_uri));

  const auto num_reader_threads = absl::GetFlag(FLAGS_num_reader_threads);
  if (num_reader_threads == 0) {
    return absl::InvalidArgumentError("Invalid number of reader threads");
  }

  const gcs::UserProject requester_pays_project(
      absl::GetFlag(FLAGS_requester_pays_project));

  StopWatch stop_watch;
  std::cout << "Reading metadata...";
  std::cout.flush();
  const auto input_bucket = std::string(input_bucket_and_path.first);
  const auto input_path = std::string(input_bucket_and_path.second);
  auto gcs_client =
      gcs::Client(google::cloud::Options{}.set<gcs::ConnectionPoolSizeOption>(
          num_reader_threads));
  auto metadata_stream = gcs_client.ReadObject(
      input_bucket, absl::StrCat(input_path, "/metadata.json"),
      requester_pays_project);
  if (metadata_stream.bad()) {
    return absl::FailedPreconditionError(absl::StrCat(
        "Failed to read metadata: ", metadata_stream.status().message()));
  }
  nlohmann::json metadata;
  metadata_stream >> metadata;
  metadata_stream.Close();
  if (metadata.is_discarded()) {
    return absl::FailedPreconditionError("Failed to parse metadata JSON");
  }
  const auto &metadata_samples = metadata["samples"];
  const uint32_t num_samples = metadata_samples.size();
  std::vector<std::string> sample_ids;
  sample_ids.reserve(num_samples);
  for (const auto &sample_id : metadata_samples) {
    sample_ids.push_back(sample_id);
  }
  const size_t num_sites = metadata["num_sites"];
  metadata.clear();
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  // Allocate the memory required for the bit sets. To tranpose the Hail
  // MatrixTable, we concatenate each sample's genotypes across all partitions,
  // storing the results sequentially in a bit set. Each sample gets its own set
  // of words in the bit set, to facilitate comparing any two samples. Each
  // genotype is split into two bits (is_het followed by is_hom_var for each
  // sample), in distinct bit sets.
  const size_t words_per_sample = 2 * CeilIntDiv(num_sites, size_t(64));
  const size_t bit_set_size = words_per_sample * num_samples;
  std::cout << "Allocating "
            << CeilIntDiv(bit_set_size * sizeof(uint64_t), size_t(1) << 20)
            << " MiB of memory for bit set...";
  std::cout.flush();
  auto bit_set = NewCudaArray<uint64_t>(bit_set_size);
  // We initialize all bits to 1 as that indicates a missing value (i.e. is_het
  // and is_hom_var are both set). That's why below we only ever have to clear
  // bits (AtomicClearBit).
  memset(bit_set.get(), 0xFF, bit_set_size * sizeof(uint64_t));
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  // Find all Parquet table files.
  std::cout << "Listing input files...";
  std::cout.flush();
  std::vector<std::pair<std::string, size_t>> input_files;
  for (const auto &blob_metadata_or_status : gcs_client.ListObjects(
           input_bucket, gcs::Prefix(input_path), requester_pays_project)) {
    ASSIGN_OR_RETURN(const auto &blob_metadata, blob_metadata_or_status);
    if (!absl::EndsWith(blob_metadata.name(), ".parquet")) {
      continue;  // Skip files that are not Parquet tables.
    }
    input_files.emplace_back(blob_metadata.name(), blob_metadata.size());
  }
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;
  std::cout << "Found " << input_files.size() << " input files." << std::endl;

  // Read Parquet tables in parallel, populating the bit set.
  std::cout << "Processing Parquet tables...";
  std::cout.flush();
  ThreadPool thread_pool(num_reader_threads);
  std::atomic<size_t> num_processed(0);
  RETURN_IF_ERROR(
      ParallelFor(&thread_pool, 0, input_files.size(), [&](const size_t i) {
        // Make a copy of the GCS client for thread-safety.
        gcs::Client gcs_client_copy = gcs_client;
        const auto &input_file = input_files[i];
        auto input_stream =
            gcs_client_copy.ReadObject(input_bucket, input_file.first);
        if (input_stream.bad()) {
          return absl::FailedPreconditionError(
              absl::StrCat("Failed to read ", input_files[i].first, ": ",
                           input_stream.status().message()));
        }
        // Read the entire file into memory, to save roundtrips to GCS.
        std::vector<char> file_buffer(input_file.second);
        input_stream.read(file_buffer.data(), file_buffer.size());
        if (input_stream.bad()) {
          return absl::FailedPreconditionError(
              absl::StrCat("Failed to read ", input_file.first, ": ",
                           input_stream.status().message()));
        }
        input_stream.Close();

        auto file_reader = parquet::ParquetFileReader::Open(
            std::make_shared<arrow::io::BufferReader>(
                reinterpret_cast<const uint8_t *>(file_buffer.data()),
                file_buffer.size()));
        const auto file_metadata = file_reader->metadata();
        constexpr size_t kNumColumns = 3;
        if (file_metadata->num_columns() != kNumColumns) {
          return absl::FailedPreconditionError(absl::StrCat(
              "Expected ", kNumColumns, " columns, found ",
              file_metadata->num_columns(), " in ", input_file.first));
        }

        // Decompress all columns into memory, this way it's easier to process a
        // full row afterwards. We expect the following schema: row_idx(INT64),
        // col_idx(INT64), n_alt_alleles(INT32).
        const size_t num_rows = file_metadata->num_rows();
        std::vector<int64_t> row_idx_buffer, col_idx_buffer;
        std::vector<int32_t> n_alt_alleles_buffer;
        row_idx_buffer.resize(num_rows);
        col_idx_buffer.resize(num_rows);
        n_alt_alleles_buffer.resize(num_rows);
        size_t row_idx_offset = 0, col_idx_offset = 0, n_alt_alleles_offset = 0;
        const size_t num_row_groups = file_metadata->num_row_groups();
        for (size_t row_group = 0; row_group < num_row_groups; ++row_group) {
          auto row_group_reader = file_reader->RowGroup(row_group);
          // row_idx (INT64)
          {
            auto column_reader = row_group_reader->Column(0);
            if (column_reader->type() != parquet::Type::INT64) {
              return absl::FailedPreconditionError(
                  absl::StrCat("Expected INT64 type, found ",
                               parquet::TypeToString(column_reader->type()),
                               " in ", input_file.first));
            }
            parquet::Int64Reader *const int64_reader =
                static_cast<parquet::Int64Reader *>(column_reader.get());
            while (int64_reader->HasNext()) {
              // We can set def_levels and rep_levels to nullptr, as there are
              // no undefined values.
              int64_t values_read = 0;
              int64_reader->ReadBatch(
                  num_rows - row_idx_offset, /* def_levels */ nullptr,
                  /* rep_levels */ nullptr,
                  row_idx_buffer.data() + row_idx_offset, &values_read);
              row_idx_offset += values_read;
            }
          }
          // col_idx (INT64)
          {
            auto column_reader = row_group_reader->Column(1);
            if (column_reader->type() != parquet::Type::INT64) {
              return absl::FailedPreconditionError(
                  absl::StrCat("Expected INT64 type, found ",
                               parquet::TypeToString(column_reader->type()),
                               " in ", input_file.first));
            }
            parquet::Int64Reader *const int64_reader =
                static_cast<parquet::Int64Reader *>(column_reader.get());
            while (int64_reader->HasNext()) {
              // We can set def_levels and rep_levels to nullptr, as there are
              // no undefined values.
              int64_t values_read = 0;
              int64_reader->ReadBatch(
                  num_rows - col_idx_offset, /* def_levels */ nullptr,
                  /* rep_levels */ nullptr,
                  col_idx_buffer.data() + col_idx_offset, &values_read);
              col_idx_offset += values_read;
            }
          }
          // n_alt_alleles (INT32)
          {
            auto column_reader = row_group_reader->Column(2);
            if (column_reader->type() != parquet::Type::INT32) {
              return absl::FailedPreconditionError(
                  absl::StrCat("Expected INT32 type, found ",
                               parquet::TypeToString(column_reader->type()),
                               " in ", input_file.first));
            }
            parquet::Int32Reader *const int32_reader =
                static_cast<parquet::Int32Reader *>(column_reader.get());
            while (int32_reader->HasNext()) {
              // We can set def_levels and rep_levels to nullptr, as there are
              // no undefined values.
              int64_t values_read = 0;
              int32_reader->ReadBatch(
                  num_rows - n_alt_alleles_offset, /* def_levels */ nullptr,
                  /* rep_levels */ nullptr,
                  n_alt_alleles_buffer.data() + n_alt_alleles_offset,
                  &values_read);
              n_alt_alleles_offset += values_read;
            }
          }
        }

        // Update the bit set now that the whole table is in memory.
        for (size_t row = 0; row < num_rows; ++row) {
          const int32_t row_idx = row_idx_buffer[row];
          const int32_t col_idx = col_idx_buffer[row];
          const int32_t n_alt_alleles = n_alt_alleles_buffer[row];
          // Pointers to the beginning of the bit set for this sample.
          uint64_t *const is_het_ptr =
              bit_set.get() + col_idx * words_per_sample;
          uint64_t *const is_hom_var_ptr = is_het_ptr + words_per_sample / 2;
          switch (n_alt_alleles) {
            case 0:  // hom-ref
              AtomicClearBit(is_het_ptr, row_idx);
              AtomicClearBit(is_hom_var_ptr, row_idx);
              break;
            case 1:  // het
              AtomicClearBit(is_hom_var_ptr, row_idx);
              break;
            case 2:  // hom-var
              AtomicClearBit(is_het_ptr, row_idx);
              break;
            default:
              return absl::FailedPreconditionError(absl::StrCat(
                  "Invalid value for n_alt_alleles (", n_alt_alleles,
                  ") encountered in ", input_file.first));
          }
        }

        if ((++num_processed & ((size_t(1) << 10) - 1)) == 0) {
          std::cout << ".";  // Progress indicator.
          std::cout.flush();
        }
        return absl::OkStatus();
      }));
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  const uint32_t max_results = absl::GetFlag(FLAGS_max_results);
  std::cout << "Allocating "
            << CeilIntDiv(max_results * sizeof(KingResult), size_t(1) << 20)
            << " MiB of memory for results...";
  std::cout.flush();
  auto results = NewCudaArray<KingResult>(max_results);
  memset(results.get(), 0, sizeof(KingResult) * max_results);
  // Result counter and overflow flag.
  auto result_index_and_flag = NewCudaArray<uint32_t>(2);
  result_index_and_flag[0] = result_index_and_flag[1] = 0;
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  std::cout << "Running KING CUDA kernel for " << num_samples << " samples...";
  std::cout.flush();
  constexpr uint64_t kCudaBlockSize = 1024;
  const uint64_t kNumCudaBlocks =
      CeilIntDiv(uint64_t(num_samples) * num_samples, kCudaBlockSize);
  ComputeKingKernel<<<kNumCudaBlocks, kCudaBlockSize>>>(
      num_samples, words_per_sample, bit_set.get(),
      absl::GetFlag(FLAGS_king_coeff_threshold), max_results, results.get(),
      &result_index_and_flag[0], &result_index_and_flag[1]);

  // Wait for GPU to finish before accessing on host.
  cudaDeviceSynchronize();
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  if (result_index_and_flag[1] != 0) {
    return absl::ResourceExhaustedError(
        absl::StrCat("Could not store all results: try increasing the "
                     "--max_results parameter."));
  }

  const uint32_t num_results = result_index_and_flag[0];
  std::cout << "Processing " << num_results << " results...";
  std::cout.flush();
  // Free some memory for postprocessing.
  bit_set.reset();
  // Create a map for JSON serialization.
  absl::flat_hash_map<std::string_view,
                      absl::flat_hash_map<std::string_view, float>>
      result_map;
  for (uint32_t i = 0; i < num_results; ++i) {
    const auto &result = results[i];
    result_map[sample_ids[result.sample_i]][sample_ids[result.sample_j]] =
        result.coeff;
  }

  ASSIGN_OR_RETURN(
      const auto output_metadata,
      gcs_client.InsertObject(std::string(output_bucket_and_path.first),
                              std::string(output_bucket_and_path.second),
                              nlohmann::json(result_map).dump(),
                              requester_pays_project));

  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;
  std::cout << "Wrote " << CeilIntDiv(output_metadata.size(), size_t(1) << 20)
            << " MiB." << std::endl;

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