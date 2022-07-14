cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(cuking LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES 80) # A100

add_compile_options(-Wall)

enable_language(CUDA)

find_package(absl REQUIRED)
find_package(Arrow REQUIRED)
find_package(google_cloud_cpp_storage REQUIRED)
find_package(nlohmann_json REQUIRED)

add_executable(cuking cuking.cu)

target_link_libraries(cuking PRIVATE
    absl::flags
    absl::flags_parse
    absl::status
    absl::statusor
    arrow_shared
    google-cloud-cpp::storage
    nlohmann_json::nlohmann_json
    parquet
)