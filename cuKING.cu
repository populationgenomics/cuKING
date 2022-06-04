#include <iostream>
#include <string>
#include <vector>

#include "absl/strings/str_join.h"
#include "vcflib/Variant.h"

__global__ void add_kernel(int n, float *x, float *y) {
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  const int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride) {
    y[i] = x[i] + y[i];
  }
}

int main(int argc, char **argv) {
  std::vector<std::string> v = {"foo", "bar", "baz"};
  std::string s = absl::StrJoin(v, "-");
  std::cout << "Joined string: " << s << "\n";

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
