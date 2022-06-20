#include <absl/base/thread_annotations.h>
#include <absl/container/flat_hash_map.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <absl/synchronization/blocking_counter.h>
#include <absl/synchronization/mutex.h>
#include <absl/time/time.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <queue>
#include <string>
#include <thread>
#include <vector>

#include "gcs_client.h"

ABSL_FLAG(std::string, sample_map, "",
          "A JSON file mapping sample IDs to cuking input paths, e.g. "
          "gs://some/bucket/sample_map.json");
ABSL_FLAG(std::string, output, "",
          "The sparse matrix result JSON output path, e.g. "
          "gs://some/bucket/relatedness.json");
ABSL_FLAG(
    uint32_t, max_results, 100 << 20,
    "How many coefficients for related sample pairs to reserve memory for.");
ABSL_FLAG(int, num_reader_threads, 100,
          "How many threads to use for parallel file reading.");
ABSL_FLAG(float, coeff_threshold, 0.1f,
          "Only store coefficients larger than this threshold.");

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
  CudaArray<uint32_t> deltas;
};

// Reads and decompresses sample data from `paths`.
std::optional<ReadSamplesResult> ReadSamples(
    const std::vector<std::string> &paths,
    cuking::GcsClient *const gcs_client) {
  ThreadPool thread_pool(absl::GetFlag(FLAGS_num_reader_threads));
  absl::BlockingCounter blocking_counter(paths.size());
  absl::Mutex mutex;
  ReadSamplesResult result ABSL_GUARDED_BY(mutex);
  std::atomic<bool> success(true);
  for (size_t i = 0; i < paths.size(); ++i) {
    thread_pool.Schedule([&, i] {
      const auto &path = paths[i];

      auto content = gcs_client->Read(path);
      if (!content.ok()) {
        std::cerr << "Error: failed to read \"" << path
                  << "\": " << content.status() << std::endl;
        success = false;
        blocking_counter.DecrementCount();
        return;
      }

      // Make sure the buffer is set and expected sizes match.
      const size_t num_entries = content->size() / sizeof(uint32_t);
      {
        const absl::MutexLock lock(&mutex);
        if (result.num_entries == 0) {
          result.num_entries = num_entries;
          result.deltas = NewCudaArray<uint32_t>(num_entries * paths.size());
        } else if (result.num_entries != num_entries) {
          std::cerr << "Mismatch for number of entries encountered for \""
                    << path << "\": " << num_entries << " vs "
                    << result.num_entries << "." << std::endl;
          success = false;
          blocking_counter.DecrementCount();
          return;
        }
      }

      // Copy to the destination buffer.
      memcpy(reinterpret_cast<char *>(result.deltas.get() + i * num_entries),
             content->data(), content->size());

      blocking_counter.DecrementCount();
    });
  }

  blocking_counter.Wait();
  return success ? std::move(result) : std::optional<ReadSamplesResult>();
}

__device__ float ComputeRareShare(const uint32_t num_entries,
                                  const uint32_t *const deltas_i,
                                  const uint32_t *const deltas_j) {
  uint32_t index_i = 0, index_j = 0;
  uint64_t position_i = deltas_i[0], position_j = deltas_j[0];
  uint32_t num_shared = 0;
  while (true) {
    if (position_i == position_j) {
      ++num_shared;
      if (++index_i == num_entries || ++index_j == num_entries) {
        break;
      }
      position_i += deltas_i[index_i];
      position_j += deltas_j[index_j];
    } else if (position_i < position_j) {
      if (++index_i == num_entries) {
        break;
      }
      position_i += deltas_i[index_i];
    } else {
      if (++index_j == num_entries) {
        break;
      }
      position_j += deltas_j[index_j];
    }
  }

  return static_cast<float>(num_shared) / num_entries;
}

// Stores the coefficient for one pair of samples.
struct RareShareResult {
  uint32_t sample_i, sample_j;
  float coeff;
};

__global__ void ComputeRareShareKernel(const uint32_t num_samples,
                                       const uint32_t num_entries,
                                       const uint32_t *const deltas,
                                       const float coeff_threshold,
                                       const uint32_t max_results,
                                       RareShareResult *const results,
                                       uint32_t *const result_index) {
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  const int i = index / num_samples;
  const int j = index % num_samples;
  if (i >= num_samples || i >= j) {
    return;
  }

  const float coeff = ComputeRareShare(num_entries, deltas + i * num_entries,
                                       deltas + j * num_entries);

  if (coeff > coeff_threshold) {
    // Reserve a result slot atomically to avoid collisions.
    const uint32_t reserved = atomicAdd(result_index, 1u);
    if (reserved < max_results) {
      RareShareResult &result = results[reserved];
      result.sample_i = i;
      result.sample_j = j;
      result.coeff = coeff;
    }
  }
}

}  // namespace

int main(int argc, char **argv) {
  absl::ParseCommandLine(argc, argv);

  const auto &sample_map_file = absl::GetFlag(FLAGS_sample_map);
  if (sample_map_file.empty()) {
    std::cerr << "Error: no sample map file specified." << std::endl;
    return 1;
  }

  auto gcs_client =
      cuking::NewGcsClient(absl::GetFlag(FLAGS_num_reader_threads));
  auto sample_map_str = gcs_client->Read(sample_map_file);
  if (!sample_map_str.ok()) {
    std::cerr << "Error: failed to read sample map file: "
              << sample_map_str.status() << std::endl;
    return 1;
  }

  auto sample_map = nlohmann::json::parse(*sample_map_str);
  sample_map_str->clear();
  std::vector<std::string> sample_ids;
  std::vector<std::string> sample_paths;
  for (const auto &[id, path] : sample_map.items()) {
    sample_ids.push_back(id);
    sample_paths.push_back(path);
  }
  sample_map.clear();

  const size_t num_samples = sample_paths.size();
  auto samples = ReadSamples(sample_paths, gcs_client.get());
  if (!samples) {
    std::cerr << "Error: failed to read samples." << std::endl;
    return 1;
  }

  std::cout << "Read " << num_samples << " samples." << std::endl;

  const uint32_t kMaxResults = absl::GetFlag(FLAGS_max_results);
  auto results = NewCudaArray<RareShareResult>(kMaxResults);
  memset(results.get(), 0, sizeof(RareShareResult) * kMaxResults);
  auto result_index = NewCudaArray<uint32_t>(1);
  result_index[0] = 0;

  const absl::Time time_before = absl::Now();

  constexpr int kCudaBlockSize = 1024;
  const int kNumCudaBlocks =
      (num_samples * num_samples + kCudaBlockSize - 1) / kCudaBlockSize;
  ComputeRareShareKernel<<<kNumCudaBlocks, kCudaBlockSize>>>(
      num_samples, samples->num_entries, samples->deltas.get(),
      absl::GetFlag(FLAGS_coeff_threshold), kMaxResults, results.get(),
      result_index.get());

  // Wait for GPU to finish before accessing on host.
  cudaDeviceSynchronize();

  const absl::Time time_after = absl::Now();
  std::cout << "CUDA kernel time: " << (time_after - time_before) << std::endl;

  // Free some memory for postprocessing.
  samples->deltas.reset();

  const uint32_t num_results = result_index[0];
  if (num_results > kMaxResults) {
    std::cerr << "Error: could not store all results: try increasing the "
                 "--max_results parameter."
              << std::endl;
    return 1;
  }

  std::cout << "Found " << num_results
            << " coefficients above the cut-off threshold." << std::endl;

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

  // Create a map for JSON serialization.
  absl::flat_hash_map<std::string_view,
                      absl::flat_hash_map<std::string_view, float>>
      result_map;
  for (size_t i = 0; i < num_results; ++i) {
    const auto &result = results[i];
    result_map[sample_ids[result.sample_i]][sample_ids[result.sample_j]] =
        result.coeff;
  }

  if (const auto status = gcs_client->Write(absl::GetFlag(FLAGS_output),
                                            nlohmann::json(result_map).dump(4));
      !status.ok()) {
    std::cerr << "Failed to write output: " << status << std::endl;
    return 1;
  }

  return 0;
}
