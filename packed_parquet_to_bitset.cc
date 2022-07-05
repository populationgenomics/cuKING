#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <absl/status/status.h>
#include <absl/synchronization/blocking_counter.h>
#include <arrow/filesystem/filesystem.h>
#include <arrow/result.h>
#include <parquet/file_reader.h>
#include <parquet/metadata.h>

#include <iostream>
#include <vector>

#include "status_macros.h"
#include "thread_pool.h"

ABSL_FLAG(std::string, input_uri, "",
          "URI containing Parquet table partitions, supports file:// as well "
          "as gs://, e.g. gs://some/bucket/my_table.parquet");
ABSL_FLAG(size_t, num_threads, 100,
          "How many threads to use for processing of Parquet partitions. This "
          "influences the amount of memory required.");

namespace status_macros {

absl::Status ToAbslStatus(const arrow::Status& status) {
  return absl::UnknownError(status.ToString());
}

}  // namespace status_macros

namespace {

// Returns ceil(a / b) for integers a, b.
template <typename T>
inline T CeilIntDiv(const T a, const T b) {
  return (a + b - 1) / b;
}

// Atomically sets a bit in a bit set.
inline void AtomicSetBit(uint64_t* const bit_set, uint64_t index) {
  uint64_t* const ptr = bit_set + (index >> 6);
  // C++20 adds std::atomic_ref, but we're compiling with C++17. We can use
  // relaxed memory ordering as bit set values don't depend on one another.
  __atomic_or_fetch(ptr, uint64_t(1) << (index & size_t(63)), __ATOMIC_RELAXED);
}

// Atomically clears a bit in a bit set.
inline void AtomicClearBit(uint64_t* const bit_set, uint64_t index) {
  uint64_t* const ptr = bit_set + (index >> 6);
  // C++20 adds std::atomic_ref, but we're compiling with C++17. We can use
  // relaxed memory ordering as bit set values don't depend on one another.
  __atomic_and_fetch(ptr, ~(uint64_t(1) << (index & size_t(63))),
                     __ATOMIC_RELAXED);
}

// Keeps track of time intervals.
class StopWatch {
 public:
  absl::Duration GetElapsedAndReset() {
    const auto now = absl::Now();
    const auto result = now - last_time_;
    last_time_ = now;
    return result;
  }

 private:
  absl::Time last_time_ = absl::Now();
};

absl::Status Run() {
  // Initialize the dataset.
  std::cout << "Listing files...";
  std::cout.flush();
  StopWatch stop_watch;
  const auto input_uri = absl::GetFlag(FLAGS_input_uri);
  std::string input_path;
  ASSIGN_OR_RETURN(const auto input_fs,
                   arrow::fs::FileSystemFromUri(input_uri, &input_path));
  arrow::fs::FileSelector file_selector;
  file_selector.base_dir = input_path;
  ASSIGN_OR_RETURN(auto file_infos, input_fs->GetFileInfo(file_selector));
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  // Only keep Parquet files.
  file_infos.erase(std::remove_if(file_infos.begin(), file_infos.end(),
                                  [](const auto& file_info) {
                                    return file_info.extension() != "parquet";
                                  }));

  if (file_infos.empty()) {
    return absl::InvalidArgumentError("No Parquet files found");
  }

  std::cout << "Found " << file_infos.size() << " files." << std::endl;

  std::cout << "Reading metadata...";
  std::cout.flush();
  std::vector<std::shared_ptr<parquet::FileMetaData>> file_metadata(
      file_infos.size());
  cuking::ThreadPool thread_pool(absl::GetFlag(FLAGS_num_threads));
  RETURN_IF_ERROR(cuking::ParallelFor(
      0, file_infos.size(),
      [&](const size_t i) {
        ASSIGN_OR_RETURN(auto input_file,
                         input_fs->OpenInputFile(file_infos[i]));
        file_metadata[i] = parquet::ReadMetaData(input_file);
        return absl::OkStatus();
      },
      &thread_pool));
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  std::vector<size_t> row_offsets;
  row_offsets.reserve(file_infos.size());
  size_t num_rows = 0;
  for (const auto& metadata : file_metadata) {
    row_offsets.push_back(num_rows);
    num_rows += metadata->num_rows();
  }
  const size_t num_cols = file_metadata.front()->num_columns();
  std::cout << "Dataset contains " << num_rows << " rows and " << num_cols
            << " columns." << std::endl;

  // Allocate the memory required for the bit sets. To tranpose the Hail
  // MatrixTable, we we concatenate all fragment columns, storing the results
  // sequentially in a bit set. Each sample gets its own set of words in the bit
  // set. Each genotype is split into two bit sets, one for is_het, another one
  // for is_hom_var. We initialize all bits to 1 as that indicates a missing
  // value (is_het and is_hom_var set simultaneously).
  const size_t words_per_sample = CeilIntDiv(num_rows, size_t(64));
  std::vector<uint64_t> is_het(words_per_sample * num_cols, uint64_t(-1));
  std::vector<uint64_t> is_hom_var(words_per_sample * num_cols, uint64_t(-1));

  // Read all the fragments fully, updating the bit sets as we go.
  // TODO: reuse metadata?

  // Write the bit sets to a single output file.
  // TODO: define header and actually write.

  return absl::OkStatus();
}

}  // namespace

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  if (const auto status = Run(); !status.ok()) {
    std::cerr << "Error: " << status << std::endl;
    return 1;
  }

  return 0;
}