#include <absl/base/thread_annotations.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <absl/synchronization/blocking_counter.h>
#include <absl/synchronization/mutex.h>
#include <absl/time/time.h>
#include <absl/types/span.h>
#include <zstd.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <queue>
#include <string>
#include <thread>
#include <vector>

#include "cuking.h"

ABSL_FLAG(std::string, sample_list, "",
          "A text file listing one .cuking.zst input file path per line.");
ABSL_FLAG(
    size_t, sample_range_begin, 0,
    "The inclusive index of the first sample to consider in the sample list.");
ABSL_FLAG(
    size_t, sample_range_end, 0,
    "The exclusive index of the last sample to consider in the sample list.");
ABSL_FLAG(int, num_reader_threads, 100,
          "How many threads to use for parallel file reading.");
ABSL_FLAG(bool, run_on_cpu, false,
          "Whether to run computation on CPU, useful for validation.");

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

// RAII-style CUDA-managed array.
template <typename T>
struct CudaArray {
  explicit CudaArray(const size_t num_elements) : num_elements(num_elements) {
    const auto err = cudaMallocManaged(&data, num_elements * sizeof(T));
    if (err) {
      std::cerr << "Error: can't allocate CUDA memory: "
                << cudaGetErrorString(err) << std::endl;
      exit(1);
    }
  }

  CudaArray &operator=(const CudaArray &) = delete;
  CudaArray(const CudaArray &) = delete;

  CudaArray(CudaArray &&other) {
    data = other.data;
    num_elements = other.num_elements;
    other.data = nullptr;
    other.num_elements = 0;
  }

  CudaArray &operator=(CudaArray &&other) {
    if (this != &other) {
      data = other.data;
      num_elements = other.num_elements;
      other.data = nullptr;
      other.num_elements = 0;
    }
    return *this;
  }

  ~CudaArray() {
    cudaFree(data);
    data = nullptr;
    num_elements = 0;
  }

  T *data = nullptr;
  size_t num_elements = 0;
};

struct Sample {
  const uint16_t *entries = nullptr;  // Not owned.
  uint32_t num_entries = 0;
  uint32_t num_hets = 0;
};

struct ReadSamplesResult {
  explicit ReadSamplesResult(size_t num_samples) : samples(num_samples) {}
  CudaArray<Sample> samples;
  std::vector<std::unique_ptr<CudaArray<uint16_t>>> buffers;
};

// Reads and decompresses sample data from `paths`, adding entries with
// corresponding index into `result`.  Returns false if any failures occurred.
bool ReadSamples(const absl::Span<std::string> &paths,
                 ThreadPool *const thread_pool,
                 ReadSamplesResult *const result) {
  result->samples = CudaArray<Sample>(paths.size());
  result->buffers.clear();
  result->buffers.resize(paths.size());

  absl::BlockingCounter blocking_counter(paths.size());
  std::atomic<bool> success(true);
  for (size_t i = 0; i < paths.size(); ++i) {
    thread_pool->Schedule([&, i] {
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

      // Read the entire file.
      std::vector<uint8_t> contents(file_size);
      in.read(reinterpret_cast<char *>(contents.data()), file_size);
      if (!in) {
        std::cerr << "Error: failed to read \"" << path << "\"." << std::endl;
        success = false;
        blocking_counter.DecrementCount();
        return;
      }
      in.close();

      // Validate the file header.
      cuking::FileHeader file_header;
      bool valid_header = true;
      if (contents.size() < sizeof(cuking::FileHeader)) {
        valid_header = false;
      } else {
        memcpy(&file_header, contents.data(), sizeof(cuking::FileHeader));
        for (size_t i = 0; i < sizeof(file_header.magic); ++i) {
          if (file_header.magic[i] != cuking::kExpectedMagic[i]) {
            valid_header = false;
            break;
          }
        }
      }

      if (!valid_header) {
        std::cerr << "Error: failed to validate header for \"" << path << "\"."
                  << std::endl;
        success = false;
        blocking_counter.DecrementCount();
        return;
      }

      // Decompress the contents.
      result->buffers[i] = std::make_unique<CudaArray<uint16_t>>(
          file_header.decompressed_size / sizeof(uint16_t));
      auto &sample = result->samples.data[i];
      sample.entries = result->buffers[i]->data;
      sample.num_entries = result->buffers[i]->num_elements;
      sample.num_hets = file_header.num_hets;
      const size_t zstd_result = ZSTD_decompress(
          result->buffers[i]->data, file_header.decompressed_size,
          contents.data() + sizeof(cuking::FileHeader),
          contents.size() - sizeof(cuking::FileHeader));
      if (ZSTD_isError(zstd_result) ||
          zstd_result != file_header.decompressed_size) {
        std::cerr << "Error: failed to decompress \"" << path << "\"."
                  << std::endl;
        success = false;
        blocking_counter.DecrementCount();
        return;
      }

      blocking_counter.DecrementCount();
    });
  }

  blocking_counter.Wait();
  return success;
}

__device__ __host__ inline uint32_t DecodeLocusDelta(const uint16_t value) {
  return value >> 2;
}

__device__ __host__ inline cuking::VariantCategory DecodeVariantCategory(
    const uint16_t value) {
  return static_cast<cuking::VariantCategory>(value & 3);
}

__device__ __host__ float ComputeKing(const Sample &sample_i,
                                      const Sample &sample_j) {
  // See https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.king.
  uint32_t num_both_het = 0, num_opposing_hom = 0;
  for (uint32_t index_i = 0, index_j = 0,
                pos_i = DecodeLocusDelta(sample_i.entries[0]),
                pos_j = DecodeLocusDelta(sample_j.entries[0]);
       pos_i != static_cast<uint32_t>(-1) ||
       pos_j != static_cast<uint32_t>(-1);) {
    if (pos_i < pos_j) {
      if (DecodeVariantCategory(sample_i.entries[index_i]) ==
          cuking::VariantCategory::kHomAlt) {
        ++num_opposing_hom;
      }

      if (++index_i < sample_i.num_entries) {
        pos_i += DecodeLocusDelta(sample_i.entries[index_i]);
      } else {
        pos_i = static_cast<uint32_t>(-1);
      }
    }

    if (pos_j < pos_i) {
      if (DecodeVariantCategory(sample_j.entries[index_j]) ==
          cuking::VariantCategory::kHomAlt) {
        ++num_opposing_hom;
      }

      if (++index_j < sample_j.num_entries) {
        pos_j += DecodeLocusDelta(sample_j.entries[index_j]);
      } else {
        pos_j = static_cast<uint32_t>(-1);
      }
    }

    if (pos_i == pos_j && pos_i != static_cast<uint32_t>(-1)) {
      if (DecodeVariantCategory(sample_i.entries[index_i]) ==
              cuking::VariantCategory::kHet &&
          DecodeVariantCategory(sample_j.entries[index_j]) ==
              cuking::VariantCategory::kHet) {
        ++num_both_het;
      }

      if (++index_i < sample_i.num_entries) {
        pos_i += DecodeLocusDelta(sample_i.entries[index_i]);
      } else {
        pos_i = static_cast<uint32_t>(-1);
      }

      if (++index_j < sample_j.num_entries) {
        pos_j += DecodeLocusDelta(sample_j.entries[index_j]);
      } else {
        pos_j = static_cast<uint32_t>(-1);
      }
    }
  }

  // Return the "between-family" estimator.
  const uint32_t min_hets = sample_i.num_hets < sample_j.num_hets
                                ? sample_i.num_hets
                                : sample_j.num_hets;
  return 0.5f + (2.f * num_both_het - 4.f * num_opposing_hom -
                 sample_i.num_hets - sample_j.num_hets) /
                    (4.f * min_hets);
}

__global__ void ComputeKingKernel(const Sample *const samples,
                                  float *const result) {
  const int i = blockIdx.x;
  const int stride = blockDim.x;
  const int num_samples = gridDim.x;
  for (int j = i + 1 + threadIdx.x; j < num_samples; j += stride) {
    result[i * num_samples + j] = ComputeKing(samples[i], samples[j]);
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
  std::string line;
  std::vector<std::string> sample_paths;
  while (std::getline(sample_list, line)) {
    if (line.empty()) {
      continue;
    }
    sample_paths.push_back(line);
  }

  const size_t sample_range_begin = absl::GetFlag(FLAGS_sample_range_begin);
  const size_t sample_range_end = absl::GetFlag(FLAGS_sample_range_end);
  if (sample_range_begin >= sample_range_end ||
      sample_range_end > sample_paths.size()) {
    std::cerr << "Error: invalid sample range specified." << std::endl;
    return 1;
  }

  const size_t num_samples = sample_range_end - sample_range_begin;
  const auto sample_paths_span =
      absl::MakeSpan(sample_paths).subspan(sample_range_begin, num_samples);
  ThreadPool thread_pool(absl::GetFlag(FLAGS_num_reader_threads));
  ReadSamplesResult read_samples_result(sample_paths_span.size());
  if (!ReadSamples(sample_paths_span, &thread_pool, &read_samples_result)) {
    std::cerr << "Error: failed to read samples." << std::endl;
    return 1;
  }

  std::cout << "Read " << num_samples << " samples." << std::endl;

  const auto &samples = read_samples_result.samples;

  if (absl::GetFlag(FLAGS_run_on_cpu)) {
    for (size_t i = 0; i < num_samples - 1; ++i) {
      for (size_t j = i + 1; j < num_samples; ++j) {
        const absl::Time time_before = absl::Now();
        const float king_coeff = ComputeKing(samples.data[i], samples.data[j]);
        const absl::Time time_after = absl::Now();
        std::cout << "KING coefficient between " << i << " and " << j << ": "
                  << king_coeff << " (took " << (time_after - time_before)
                  << ")";
        // Cut off at third degree
        // (https://www.kingrelatedness.com/manual.shtml).
        constexpr float kKingCutoff = 0.0884f;
        if (king_coeff < kKingCutoff) {
          std::cout << " (unrelated)";
        }
        std::cout << std::endl;
      }
    }
  } else {
    CudaArray<float> result(num_samples * num_samples);
    for (size_t i = 0; i < result.num_elements; ++i) {
      result.data[i] = 0.f;
    }

    ComputeKingKernel<<<num_samples, 64>>>(samples.data, result.data);

    // Wait for GPU to finish before accessing on host.
    cudaDeviceSynchronize();

    for (size_t i = 0; i < num_samples - 1; ++i) {
      for (size_t j = i + 1; j < num_samples; ++j) {
        std::cout << result.data[i * num_samples + j] << std::endl;
      }
    }
  }

  return 0;
}
