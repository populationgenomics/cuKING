#include <absl/base/thread_annotations.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <absl/synchronization/blocking_counter.h>
#include <absl/synchronization/mutex.h>
#include <zstd.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <thread>
#include <vector>

ABSL_FLAG(std::string, sample_list, "",
          "A text file listing one .cuking.zst input file path per line.");
ABSL_FLAG(int, num_reader_threads, 16,
          "How many threads to use for parallel file reading.");

__global__ void add_kernel(int n, float *x, float *y) {
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  const int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride) {
    y[i] = x[i] + y[i];
  }
}

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

using Sample = std::vector<uint16_t>;

// Reads and decompresses sample buffers from `paths`, until either all paths
// are read or the `max_total_size` has been reached. Callers should check the
// length of the `result` vector to determine how many samples have been read.
// Returns false if any failures occurred.
bool ReadSamples(const std::vector<std::string> &paths,
                 const size_t max_total_size, ThreadPool *const thread_pool,
                 std::vector<Sample> *const result) {
  result->clear();
  size_t total_size = 0;
  absl::Mutex mutex;
  absl::BlockingCounter blocking_counter(paths.size());
  bool errors = false;
  for (const std::string &path : paths) {
    thread_pool->Schedule([&] {
      // Early-out if we've already read too much.
      {
        absl::MutexLock lock(&mutex);
        if (total_size >= max_total_size) {
          blocking_counter.DecrementCount();
          return;
        }
      }

      // Determine the file size.
      std::error_code error_code;
      const size_t file_size = std::filesystem::file_size(path, error_code);
      if (error_code) {
        std::cerr << "Error: failed to determine size of \"" << path
                  << "\": " << error_code << std::endl;
        {
          absl::MutexLock lock(&mutex);
          errors = true;
        }
        blocking_counter.DecrementCount();
        return;
      }

      std::ifstream in(path, std::ifstream::binary);
      if (!in) {
        std::cerr << "Error: failed to open \"" << path << "\"." << std::endl;
        {
          absl::MutexLock lock(&mutex);
          errors = true;
        }
        blocking_counter.DecrementCount();
        return;
      }

      // Read the entire file.
      std::vector<uint8_t> buffer(file_size);
      in.read(reinterpret_cast<char *>(buffer.data()), file_size);
      if (!in) {
        std::cerr << "Error: failed to read \"" << path << "\"." << std::endl;
        {
          absl::MutexLock lock(&mutex);
          errors = true;
        }
        blocking_counter.DecrementCount();
        return;
      }
      in.close();

      // Validate the file header.
      static constexpr char magic[] = "CUK1";
      bool valid_header = true;
      if (buffer.size() < sizeof(magic) + sizeof(size_t)) {
        valid_header = false;
      } else {
        for (size_t i = 0; i < sizeof(magic); ++i) {
          if (buffer[i] != magic[i]) {
            valid_header = false;
            break;
          }
        }
      }

      if (!valid_header) {
        std::cerr << "Error: failed to validate header for \"" << path << "\"."
                  << std::endl;
        {
          absl::MutexLock lock(&mutex);
          errors = true;
        }
        blocking_counter.DecrementCount();
        return;
      }

      // Decompress the contents.
      size_t decompressed_size = 0;
      memcpy(&decompressed_size, buffer.data() + sizeof(magic),
             sizeof(decompressed_size));

      Sample sample(decompressed_size);
      const size_t zstd_result =
          ZSTD_decompress(sample.data(), decompressed_size,
                          buffer.data() + sizeof(magic) + sizeof(size_t),
                          buffer.size() - sizeof(magic) - sizeof(size_t));
      if (ZSTD_isError(zstd_result)) {
        std::cerr << "Error: failed to decompress \"" << path << "\"."
                  << std::endl;
        {
          absl::MutexLock lock(&mutex);
          errors = true;
        }
        blocking_counter.DecrementCount();
        return;
      }

      // Commit the result.
      {
        absl::MutexLock lock(&mutex);
        if (total_size + decompressed_size <= max_total_size) {
          total_size += decompressed_size;
          result->push_back(std::move(sample));
        }
        blocking_counter.DecrementCount();
      }
    });
  }

  blocking_counter.Wait();
  return errors;
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

  ThreadPool thread_pool(absl::GetFlag(FLAGS_num_reader_threads));
  std::vector<Sample> samples;
  if (!ReadSamples(sample_paths, static_cast<size_t>(-1), &thread_pool,
                   &samples)) {
    std::cerr << "Error: failed to read samples." << std::endl;
    return 1;
  }

  std::cout << "Read " << samples.size() << " samples." << std::endl;

  constexpr int N = 1 << 20;
  float *x, *y;

  // Allocate Unified Memory – accessible from CPU or GPU.
  cudaMallocManaged(&x, N * sizeof(float));
  cudaMallocManaged(&y, N * sizeof(float));

  // Initialize x and y arrays on the host.
  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }

  // Run kernel on 1M elements on the GPU.
  constexpr int blockSize = 256;
  constexpr int numBlocks = (N + blockSize - 1) / blockSize;
  add_kernel<<<numBlocks, blockSize>>>(N, x, y);

  // Wait for GPU to finish before accessing on host.
  cudaDeviceSynchronize();

  // Check for errors (all values should be 3.0f).
  float maxError = 0.0f;
  for (int i = 0; i < N; i++) {
    maxError = std::fmax(maxError, std::fabs(y[i] - 3.0f));
  }
  std::cout << "Max error: " << maxError << std::endl;

  cudaFree(x);
  cudaFree(y);

  return 0;
}
