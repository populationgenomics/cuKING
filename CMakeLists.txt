cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(cuKING LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(ABSL_PROPAGATE_CXX_STD ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES 60)  # P100

enable_language(CUDA)

add_subdirectory(abseil-cpp)

add_executable(main main.cu)

target_link_libraries(main absl::strings)