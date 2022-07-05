#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <absl/status/status.h>
#include <absl/synchronization/blocking_counter.h>
#include <arrow/filesystem/filesystem.h>
#include <arrow/result.h>
#include <parquet/column_reader.h>
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
  // Initialize the input dataset.
  std::cout << "Listing input files...";
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
    return absl::FailedPreconditionError("No Parquet files found");
  }
  std::cout << "Found " << file_infos.size() << " input files." << std::endl;

  // Read the metadata only to determine the row count per partition.
  std::cout << "Reading metadata...";
  std::cout.flush();
  std::vector<std::shared_ptr<parquet::FileMetaData>> file_metadata(
      file_infos.size());
  cuking::ThreadPool thread_pool(absl::GetFlag(FLAGS_num_threads));
  std::atomic<size_t> num_processed(0);
  RETURN_IF_ERROR(cuking::ParallelFor(
      0, file_infos.size(),
      [&](const size_t i) {
        ASSIGN_OR_RETURN(auto input_file,
                         input_fs->OpenInputFile(file_infos[i]));
        file_metadata[i] = parquet::ReadMetaData(input_file);
        if ((++num_processed & ((size_t(1) << 10) - 1)) == 0) {
          std::cout << ".";  // Progress indicator.
          std::cout.flush();
        }
        return absl::OkStatus();
      },
      &thread_pool));
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  // Compute the overall row offset per partition (a reduction).
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
  // MatrixTable, we concatenate each sample's columns across all partitions,
  // storing the results sequentially in a bit set. Each sample gets its own set
  // of words in the bit set, to facilitate comparing any two samples. Each
  // genotype is split into two bits (is_het / is_hom_var), in distinct bit
  // sets. We initialize all bits to 1 as that indicates a missing value (i.e.
  // is_het and is_hom_var are both set). That's why below we only ever have to
  // clear bits (AtomicClearBit).
  const size_t words_per_sample = CeilIntDiv(num_rows, size_t(64));
  const size_t bit_set_size = words_per_sample * num_cols;
  std::cout << "Allocating "
            << CeilIntDiv(bit_set_size * sizeof(uint64_t), size_t(1) << 30)
            << " GiB of memory for bit sets...";
  std::cout.flush();
  std::vector<uint64_t> is_het(words_per_sample * num_cols, uint64_t(-1));
  std::vector<uint64_t> is_hom_var(words_per_sample * num_cols, uint64_t(-1));
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  // Read all the partitions fully, updating the bit sets as we go.
  std::cout << "Processing partitions...";
  std::cout.flush();
  num_processed = 0;
  RETURN_IF_ERROR(cuking::ParallelFor(
      0, file_infos.size(),
      [&](const size_t i) {
        const auto& metadata = file_metadata[i];
        if (size_t(metadata->num_columns()) != num_cols) {
          return absl::FailedPreconditionError(absl::StrCat(
              "Expected ", num_cols, " columns, found ",
              metadata->num_columns(), " in ", file_infos[i].path()));
        }

        ASSIGN_OR_RETURN(auto input_file,
                         input_fs->OpenInputFile(file_infos[i]));
        auto file_reader = parquet::ParquetFileReader::Open(
            std::move(input_file), parquet::default_reader_properties(),
            metadata);

        static constexpr size_t kBatchBufferSize = 1024;
        std::vector<int32_t> batch_buffer(
            kBatchBufferSize);  // Used to read a column.
        const size_t num_row_groups = metadata->num_row_groups();
        for (size_t row_group = 0; row_group < num_row_groups; ++row_group) {
          auto row_group_reader = file_reader->RowGroup(row_group);
          for (size_t column = 0; column < num_cols; ++column) {
            auto column_reader = row_group_reader->Column(column);
            if (column_reader->type() != parquet::Type::INT32) {
              return absl::FailedPreconditionError(
                  absl::StrCat("Expected INT32 type, found ",
                               parquet::TypeToString(column_reader->type()),
                               " in ", file_infos[i].path()));
            }
            parquet::Int32Reader* const int32_reader =
                static_cast<parquet::Int32Reader*>(column_reader.get());

            // Pointers to the beginning of the bit set for this sample.
            uint64_t* const is_het_ptr =
                is_het.data() + column * words_per_sample;
            uint64_t* const is_hom_var_ptr =
                is_hom_var.data() + column * words_per_sample;

            // Read all values in this column.
            size_t row = row_offsets[i];
            while (int32_reader->HasNext()) {
              // We can set def_levels and rep_levels to nullptr, as there are
              // no undefined values.
              int64_t values_read = 0;
              int32_reader->ReadBatch(
                  kBatchBufferSize, /* def_levels */ nullptr,
                  /* rep_levels */ nullptr, batch_buffer.data(), &values_read);
              for (int64_t k = 0; k < values_read; ++k, ++row) {
                switch (batch_buffer[k]) {
                  case 0:  // hom-ref
                    AtomicClearBit(is_het_ptr, row);
                    AtomicClearBit(is_hom_var_ptr, row);
                    break;
                  case 1:  // het
                    AtomicClearBit(is_hom_var_ptr, row);
                    break;
                  case 2:  // hom-var
                    AtomicClearBit(is_het_ptr, row);
                    break;
                  case 3:  // missing
                    break;
                  default:
                    return absl::FailedPreconditionError(
                        absl::StrCat("Invalid value ", batch_buffer[k],
                                     " encountered in ", file_infos[i].path()));
                }
              }
            }
          }
        }

        if ((++num_processed & ((size_t(1) << 10) - 1)) == 0) {
          std::cout << ".";  // Progress indicator.
          std::cout.flush();
        }
        return absl::OkStatus();
      },
      &thread_pool));
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  // Write the bit sets to a single output file, with a JSON header.
  std::cout << "Writing output...";
  std::cout.flush();
  // TODO: define header and actually write.
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  return absl::OkStatus();
}

}  // namespace

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  if (const auto status = Run(); !status.ok()) {
    std::cerr << std::endl << "Error: " << status << std::endl;
    return 1;
  }

  return 0;
}