cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(cuking LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(ABSL_PROPAGATE_CXX_STD ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES 80)  # A100

add_compile_options(-Wall)

enable_language(CUDA)

find_package(absl REQUIRED)

add_executable(cuking cuking.cu)
add_executable(gvcf2cuking gvcf2cuking.cc)

add_library(url_reader
    url_reader.cc
)

target_link_libraries(cuking PRIVATE
    absl::flags
    absl::flags_parse
    zstd
)

target_link_libraries(gvcf2cuking PRIVATE
    absl::flags
    absl::flags_parse
    absl::flat_hash_map
    hts
    url_reader
    zstd
)