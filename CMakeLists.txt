cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(cuking LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(ABSL_PROPAGATE_CXX_STD ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES 80) # A100

add_compile_options(-Wall)

enable_language(CUDA)

find_package(absl REQUIRED)
find_package(Arrow REQUIRED)
find_package(nlohmann_json REQUIRED)

add_executable(packed_parquet_to_bitset packed_parquet_to_bitset.cc)

add_executable(cuking cuking.cu)

target_link_libraries(packed_parquet_to_bitset PRIVATE
    absl::flags
    absl::flags_parse
    absl::status
    absl::statusor
    arrow_shared
    nlohmann_json::nlohmann_json
    parquet
)

target_link_libraries(cuking PRIVATE
    absl::flags
    absl::flags_parse
    absl::status
    arrow_shared
    nlohmann_json::nlohmann_json
)