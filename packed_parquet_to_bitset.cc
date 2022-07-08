#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <absl/status/status.h>
#include <absl/status/statusor.h>
#include <absl/strings/str_cat.h>
#include <absl/synchronization/blocking_counter.h>
#include <arrow/filesystem/filesystem.h>
#include <arrow/io/memory.h>
#include <arrow/result.h>
#include <arrow/status.h>
#include <parquet/column_reader.h>
#include <parquet/file_reader.h>
#include <parquet/metadata.h>

#include <algorithm>
#include <iostream>
#include <nlohmann/json.hpp>
#include <queue>
#include <thread>
#include <vector>

#include "utils.h"

ABSL_FLAG(std::string, input_uri, "",
          "URI containing Parquet table partitions. Supports file:// as well "
          "as gs://, e.g. gs://some/bucket/my_table.parquet");
ABSL_FLAG(std::string, output_uri, "",
          "URI of a directory to write outputs to. Supports file:// as well "
          "as gs://, e.g. gs://some/bucket/my_table.cuking");
ABSL_FLAG(size_t, num_threads, 100,
          "How many threads to use for processing of Parquet partitions. This "
          "influences the amount of memory required.");

namespace cuking {

absl::Status ToAbslStatus(const arrow::Status& status) {
  return absl::UnknownError(status.ToString());
}

}  // namespace cuking

namespace {

// Atomically clears a bit in a bit set.
inline void AtomicClearBit(uint64_t* const bit_set, uint64_t index) {
  uint64_t* const ptr = bit_set + (index >> 6);
  // C++20 adds std::atomic_ref, but we're compiling with C++17. We can use
  // relaxed memory ordering as bit set values don't depend on one another.
  __atomic_and_fetch(ptr, ~(uint64_t(1) << (index & size_t(63))),
                     __ATOMIC_RELAXED);
}

// Adapted from the Abseil thread pool.
class ThreadPool {
 public:
  explicit ThreadPool(const size_t num_threads) {
    for (size_t i = 0; i < num_threads; ++i) {
      threads_.push_back(std::thread(&ThreadPool::WorkLoop, this));
    }
  }

  ThreadPool(const ThreadPool&) = delete;
  ThreadPool& operator=(const ThreadPool&) = delete;

  ~ThreadPool() {
    {
      absl::MutexLock l(&mu_);
      for (size_t i = 0; i < threads_.size(); i++) {
        queue_.push(nullptr);  // Shutdown signal.
      }
    }
    for (auto& t : threads_) {
      t.join();
    }
  }

  // Schedule a function to be run on a ThreadPool thread immediately.
  void Schedule(std::function<void()> func) {
    assert(func != nullptr);
    absl::MutexLock l(&mu_);
    queue_.push(std::move(func));
  }

 private:
  bool WorkAvailable() const ABSL_EXCLUSIVE_LOCKS_REQUIRED(mu_) {
    return !queue_.empty();
  }

  void WorkLoop() {
    while (true) {
      std::function<void()> func;
      {
        absl::MutexLock l(&mu_);
        mu_.Await(absl::Condition(this, &ThreadPool::WorkAvailable));
        func = std::move(queue_.front());
        queue_.pop();
      }
      if (func == nullptr) {  // Shutdown signal.
        break;
      }
      func();
    }
  }

  absl::Mutex mu_;
  std::queue<std::function<void()>> queue_ ABSL_GUARDED_BY(mu_);
  std::vector<std::thread> threads_;
};

// Invokes the given function in parallel over a range. If any invocations
// returns an error, only one such error will be returned
// (non-deterministically).
absl::Status ParallelFor(ThreadPool* const thread_pool, const size_t begin,
                         const size_t end,
                         std::function<absl::Status(size_t index)> func) {
  absl::BlockingCounter blocking_counter(end - begin);
  absl::Mutex mu;
  absl::Status result ABSL_GUARDED_BY(mu);
  for (size_t i = begin; i < end; ++i) {
    thread_pool->Schedule([&, i] {
      auto status = func(i);
      if (!status.ok()) {
        const absl::MutexLock lock(&mu);
        result.Update(std::move(status));
      }
      blocking_counter.DecrementCount();
    });
  }
  blocking_counter.Wait();
  return result;
}

absl::StatusOr<std::vector<uint8_t>> ReadFile(
    arrow::fs::FileSystem* const file_system,
    const arrow::fs::FileInfo& file_info) {
  ASSIGN_OR_RETURN(auto file, file_system->OpenInputFile(file_info));
  const auto size = file_info.size();
  std::vector<uint8_t> buffer(size);
  ASSIGN_OR_RETURN(const auto bytes_read, file->ReadAt(0, size, buffer.data()));
  if (bytes_read != size) {
    return absl::FailedPreconditionError(
        absl::StrCat("Expected to read ", size, " bytes, but got ", bytes_read,
                     " bytes instead"));
  }
  return std::move(buffer);
}

absl::Status Run() {
  // Validate flags.
  const auto input_uri = absl::GetFlag(FLAGS_input_uri);
  if (input_uri.empty()) {
    return absl::InvalidArgumentError("No input URI specified");
  }
  const auto output_uri = absl::GetFlag(FLAGS_output_uri);
  if (output_uri.empty()) {
    return absl::InvalidArgumentError("No output URI specified");
  }
  const auto num_threads = absl::GetFlag(FLAGS_num_threads);
  if (num_threads == 0) {
    return absl::InvalidArgumentError("Invalid number of threads requested");
  }

  // Initialize the input dataset.
  std::cout << "Listing input files...";
  std::cout.flush();
  cuking::StopWatch stop_watch;
  std::string input_path;
  ASSIGN_OR_RETURN(const auto input_fs,
                   arrow::fs::FileSystemFromUri(input_uri, &input_path));
  arrow::fs::FileSelector file_selector;
  file_selector.base_dir = input_path;
  ASSIGN_OR_RETURN(auto file_infos, input_fs->GetFileInfo(file_selector));

  // Only keep Parquet files.
  file_infos.erase(std::remove_if(file_infos.begin(), file_infos.end(),
                                  [](const auto& file_info) {
                                    return file_info.extension() != "parquet";
                                  }));
  if (file_infos.empty()) {
    return absl::FailedPreconditionError("No Parquet files found");
  }

  // Sorting isn't strictly necessary, but will make outputs deterministic.
  std::sort(
      file_infos.begin(), file_infos.end(),
      [](const auto& lhs, const auto& rhs) { return lhs.path() < rhs.path(); });

  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;
  std::cout << "Found " << file_infos.size() << " input files." << std::endl;

  // Read the metadata only to determine the row count per partition.
  std::cout << "Reading metadata...";
  std::cout.flush();
  std::vector<std::shared_ptr<parquet::FileMetaData>> file_metadata(
      file_infos.size());
  ThreadPool thread_pool(num_threads);
  std::atomic<size_t> num_processed(0);
  RETURN_IF_ERROR(
      ParallelFor(&thread_pool, 0, file_infos.size(), [&](const size_t i) {
        // Reading the entire file is significantly faster than using a standard
        // RandomAccessFile on GCS, probably due to less roundtrips.
        ASSIGN_OR_RETURN(auto file_buffer,
                         ReadFile(input_fs.get(), file_infos[i]));
        auto buffer_reader = std::make_shared<arrow::io::BufferReader>(
            file_buffer.data(), file_buffer.size());
        auto file_reader =
            parquet::ParquetFileReader::Open(std::move(buffer_reader));
        file_metadata[i] = file_reader->metadata();
        if ((++num_processed & ((size_t(1) << 10) - 1)) == 0) {
          std::cout << ".";  // Progress indicator.
          std::cout.flush();
        }
        return absl::OkStatus();
      }));
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
  // genotype is split into two bits (is_het followed by is_hom_var for each
  // sample), in distinct bit sets. We initialize all bits to 1 as that
  // indicates a missing value (i.e. is_het and is_hom_var are both set). That's
  // why below we only ever have to clear bits (AtomicClearBit).
  const size_t words_per_sample = 2 * cuking::CeilIntDiv(num_rows, size_t(64));
  const size_t bit_set_size = words_per_sample * num_cols;
  std::cout << "Allocating "
            << cuking::CeilIntDiv(bit_set_size * sizeof(uint64_t), size_t(1)
                                                                       << 20)
            << " MiB of memory for bit sets...";
  std::cout.flush();
  std::vector<uint64_t> bit_set(words_per_sample * num_cols, uint64_t(-1));
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  // Read all the partitions fully, updating the bit sets as we go.
  std::cout << "Processing partitions...";
  std::cout.flush();
  num_processed = 0;
  RETURN_IF_ERROR(
      ParallelFor(&thread_pool, 0, file_infos.size(), [&](const size_t i) {
        const auto& metadata = file_metadata[i];
        if (size_t(metadata->num_columns()) != num_cols) {
          return absl::FailedPreconditionError(absl::StrCat(
              "Expected ", num_cols, " columns, found ",
              metadata->num_columns(), " in ", file_infos[i].path()));
        }

        // Reading the entire file is significantly faster than using a standard
        // RandomAccessFile on GCS, probably due to less roundtrips.
        ASSIGN_OR_RETURN(auto file_buffer,
                         ReadFile(input_fs.get(), file_infos[i]));
        auto buffer_reader = std::make_shared<arrow::io::BufferReader>(
            file_buffer.data(), file_buffer.size());
        auto file_reader = parquet::ParquetFileReader::Open(
            std::move(buffer_reader), parquet::default_reader_properties(),
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
                bit_set.data() + column * words_per_sample;
            uint64_t* const is_hom_var_ptr = is_het_ptr + words_per_sample / 2;

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
      }));
  std::cout << " (" << stop_watch.GetElapsedAndReset() << ")" << std::endl;

  // Write both a JSON metadata descriptor as well as the packed bit sets to the
  // output directory.
  std::cout << "Writing output...";
  std::cout.flush();
  std::string output_path;
  ASSIGN_OR_RETURN(const auto output_fs,
                   arrow::fs::FileSystemFromUri(output_uri, &output_path));
  RETURN_IF_ERROR(output_fs->CreateDir(output_path));
  while (output_path.back() == '/') {
    output_path.pop_back();  // Remove trailing slashes.
  }

  auto* const schema = file_metadata.front()->schema();
  std::vector<std::string_view> sample_names;
  sample_names.reserve(num_cols);
  for (size_t i = 0; i < num_cols; ++i) {
    sample_names.push_back(schema->Column(i)->name());
  }

  nlohmann::json metadata;
  metadata["samples"] = std::move(sample_names);
  metadata["num_sites"] = num_rows;
  metadata["words_per_sample"] = words_per_sample;
  metadata["creation_time"] = absl::FormatTime(absl::Now());

  ASSIGN_OR_RETURN(
      auto metadata_output_stream,
      output_fs->OpenOutputStream(absl::StrCat(output_path, "/metadata.json")));
  RETURN_IF_ERROR(metadata_output_stream->Write(metadata.dump()));
  RETURN_IF_ERROR(metadata_output_stream->Close());

  ASSIGN_OR_RETURN(
      auto bit_set_output_stream,
      output_fs->OpenOutputStream(absl::StrCat(output_path, "/bit_set.bin")));
  RETURN_IF_ERROR(bit_set_output_stream->Write(
      bit_set.data(), bit_set.size() * sizeof(uint64_t)));
  RETURN_IF_ERROR(bit_set_output_stream->Close());

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