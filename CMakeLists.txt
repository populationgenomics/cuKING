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
find_package(google_cloud_cpp_storage REQUIRED)
find_package(nlohmann_json REQUIRED)

add_executable(packed_parquet_to_bitset packed_parquet_to_bitset.cc)

add_executable(cuking cuking.cu)

add_library(gcs_client
    gcs_client.cc
)

target_link_libraries(gcs_client PRIVATE
    absl::status
    absl::statusor
    google-cloud-cpp::storage
)

target_link_libraries(packed_parquet_to_bitset PRIVATE
    absl::flags
    absl::flags_parse
    absl::status
    arrow_shared
    parquet
)

target_link_libraries(cuking PRIVATE
    absl::flags
    absl::flags_parse
    gcs_client
    nlohmann_json::nlohmann_json
)