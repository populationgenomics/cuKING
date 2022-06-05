#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

ABSL_FLAG(std::string, sample_list, "",
          "A text file listing one .cuking.zst input file path per line.");

__global__ void add_kernel(int n, float *x, float *y) {
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  const int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride) {
    y[i] = x[i] + y[i];
  }
}

namespace {

using Sample = std::vector<uint16_t>;

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
  std::vector<std::string> sample_files;
  while (std::getline(sample_list, line)) {
    if (line.empty()) {
      continue;
    }
    sample_files.push_back(line);
  }

  for (const std::string &sample_file : sample_files) {
    std::cout << sample_file << ": " << std::filesystem::file_size(sample_file)
              << std::endl;
  }

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
