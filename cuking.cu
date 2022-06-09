#include <absl/base/thread_annotations.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <absl/synchronization/blocking_counter.h>
#include <absl/synchronization/mutex.h>
#include <absl/time/time.h>
#include <absl/types/span.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <queue>
#include <string>
#include <thread>
#include <vector>

ABSL_FLAG(std::string, sample_list, "",
          "A text file listing one .cuking input file path per line.");
ABSL_FLAG(size_t, sample_range_begin, 0,
          "The inclusive index of the first sample to consider in the sample "
          "list. Defaults to the beginning of the list.");
ABSL_FLAG(size_t, sample_range_end, static_cast<size_t>(-1),
          "The exclusive index of the last sample to consider in the sample "
          "list. Defaults to the full list.");
ABSL_FLAG(int, num_reader_threads, 100,
          "How many threads to use for parallel file reading.");
ABSL_FLAG(
    uint32_t, max_results, 100 << 20,
    "How many coefficients for related sample pairs to reserve memory for.");

namespace {
// Adapted from the Abseil thread pool.
class ThreadPool {
 public:
  explicit ThreadPool(const int num_threads) {
    assert(num_threads > 0);
    for (int i = 0; i < num_threads; ++i) {
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

struct ReadSamplesResult {
  uint32_t num_entries = 0;
  CudaArray<uint64_t> bit_sets;
};

// Reads and decompresses sample data from `paths`.
std::optional<ReadSamplesResult> ReadSamples(
    const absl::Span<std::string> &paths) {
  ThreadPool thread_pool(absl::GetFlag(FLAGS_num_reader_threads));
  absl::BlockingCounter blocking_counter(paths.size());
  absl::Mutex mutex;
  ReadSamplesResult result ABSL_GUARDED_BY(mutex);
  std::atomic<bool> success(true);
  for (size_t i = 0; i < paths.size(); ++i) {
    thread_pool.Schedule([&, i] {
      const auto &path = paths[i];
      // Determine the file size.
      std::error_code error_code;
      const size_t file_size = std::filesystem::file_size(path, error_code);
      if (error_code) {
        std::cerr << "Error: failed to determine size of \"" << path
                  << "\": " << error_code << std::endl;
        success = false;
        blocking_counter.DecrementCount();
        return;
      }

      std::ifstream in(path, std::ifstream::binary);
      if (!in) {
        std::cerr << "Error: failed to open \"" << path << "\"." << std::endl;
        success = false;
        blocking_counter.DecrementCount();
        return;
      }

      // Make sure the buffer is set and expected sizes match.
      const size_t num_entries = file_size / sizeof(uint64_t) / 2;
      {
        const absl::MutexLock lock(&mutex);
        if (result.num_entries == 0) {
          result.num_entries = num_entries;
          result.bit_sets =
              NewCudaArray<uint64_t>(num_entries * 2 * paths.size());
        } else if (result.num_entries != num_entries) {
          std::cerr << "Mismatch for number of entries encountered for \""
                    << path << "\": " << num_entries << " vs "
                    << result.num_entries << "." << std::endl;
          success = false;
          blocking_counter.DecrementCount();
          return;
        }
      }

      // Read the entire file.
      in.read(
          reinterpret_cast<char *>(result.bit_sets.get() + i * 2 * num_entries),
          file_size);
      if (!in) {
        std::cerr << "Error: failed to read \"" << path << "\"." << std::endl;
        success = false;
        blocking_counter.DecrementCount();
        return;
      }
      in.close();

      blocking_counter.DecrementCount();
    });
  }

  blocking_counter.Wait();
  return success ? std::move(result) : std::optional<ReadSamplesResult>();
}

__device__ float ComputeKing(const uint32_t num_entries,
                             const uint64_t *const het_i_entries,
                             const uint64_t *const hom_alt_i_entries,
                             const uint64_t *const het_j_entries,
                             const uint64_t *const hom_alt_j_entries) {
  // See https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.king.
  uint32_t num_het_i = 0, num_het_j = 0, num_both_het = 0, num_opposing_hom = 0;
  for (uint32_t k = 0; k < num_entries; ++k) {
    const uint64_t het_i = het_i_entries[k];
    const uint64_t hom_alt_i = hom_alt_i_entries[k];
    const uint64_t het_j = het_j_entries[k];
    const uint64_t hom_alt_j = hom_alt_j_entries[k];
    const uint64_t hom_ref_i = (~hom_alt_i) & (~het_i);
    const uint64_t hom_ref_j = (~hom_alt_j) & (~het_j);
    num_het_i += __popcll(het_i);
    num_het_j += __popcll(het_j);
    num_both_het += __popcll(het_i & het_j);
    num_opposing_hom +=
        __popcll((hom_ref_i & hom_alt_j) | (hom_ref_j & hom_alt_i));
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
                                  const uint32_t num_entries,
                                  const uint64_t *const bit_sets,
                                  const uint32_t max_results,
                                  KingResult *const results,
                                  uint32_t *const result_index) {
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  const int i = index / num_samples;
  const int j = index % num_samples;
  if (i >= num_samples || i >= j) {
    return;
  }

  const float coeff = ComputeKing(num_entries, bit_sets + i * 2 * num_entries,
                                  bit_sets + (i * 2 + 1) * num_entries,
                                  bit_sets + j * 2 * num_entries,
                                  bit_sets + (j * 2 + 1) * num_entries);

  // Only store results at 3rd degree or closer.
  // (https://www.kingrelatedness.com/manual.shtml).
  constexpr float kKingCutoff = 0.0442f;
  if (coeff > kKingCutoff) {
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

}  // namespace

int main(int argc, char **argv) {
  absl::ParseCommandLine(argc, argv);

  const auto &sample_list_file = absl::GetFlag(FLAGS_sample_list);
  if (sample_list_file.empty()) {
    std::cerr << "Error: no sample list file specified." << std::endl;
    return 1;
  }

  std::ifstream sample_list(sample_list_file);
  if (!sample_list) {
    std::cerr << "Error: failed to open sample list file." << std::endl;
    return 1;
  }
  std::string line;
  std::vector<std::string> sample_paths;
  while (std::getline(sample_list, line)) {
    if (line.empty()) {
      continue;
    }
    sample_paths.push_back(line);
  }

  const size_t sample_range_begin = absl::GetFlag(FLAGS_sample_range_begin);
  const size_t sample_range_end =
      std::min(absl::GetFlag(FLAGS_sample_range_end), sample_paths.size());
  if (sample_range_begin >= sample_range_end) {
    std::cerr << "Error: invalid sample range specified." << std::endl;
    return 1;
  }

  const size_t num_samples = sample_range_end - sample_range_begin;
  const auto sample_paths_span =
      absl::MakeSpan(sample_paths).subspan(sample_range_begin, num_samples);
  const auto samples = ReadSamples(sample_paths_span);
  if (!samples) {
    std::cerr << "Error: failed to read samples." << std::endl;
    return 1;
  }

  std::cout << "Read " << num_samples << " samples." << std::endl;

  const uint32_t kMaxResults = absl::GetFlag(FLAGS_max_results);
  auto results = NewCudaArray<KingResult>(kMaxResults);
  memset(results.get(), 0, sizeof(KingResult) * kMaxResults);
  auto result_index = NewCudaArray<uint32_t>(1);
  result_index[0] = 0;

  const absl::Time time_before = absl::Now();

  constexpr int kCudaBlockSize = 1024;
  const int kNumCudaBlocks =
      (num_samples * num_samples + kCudaBlockSize - 1) / kCudaBlockSize;
  ComputeKingKernel<<<kNumCudaBlocks, kCudaBlockSize>>>(
      num_samples, samples->num_entries, samples->bit_sets.get(), kMaxResults,
      results.get(), result_index.get());

  // Wait for GPU to finish before accessing on host.
  cudaDeviceSynchronize();

  const absl::Time time_after = absl::Now();

  const uint32_t num_results = result_index[0];
  if (num_results > kMaxResults) {
    std::cerr
        << "Error: could not store all results: try increasing --max_results"
        << std::endl;
    return 1;
  }

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

  std::cout << num_related << " related samples found." << std::endl;
  std::cout << "CUDA kernel time: " << (time_after - time_before) << std::endl;

  return 0;
}
