cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(cuKING LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(ABSL_PROPAGATE_CXX_STD ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES 60)  # P100

add_compile_options(-Wall)

enable_language(CUDA)

add_subdirectory(abseil-cpp)

add_executable(cuKING cuKING.cu)
add_executable(gvcf2cuking gvcf2cuking.cc)

target_link_libraries(cuKING
    absl::flags
    absl::flags_parse
)

target_link_libraries(gvcf2cuking
    absl::flags
    absl::flags_parse
    absl::flat_hash_map
    hts
)