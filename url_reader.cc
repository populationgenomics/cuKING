#include "url_reader.h"

#include <absl/strings/strip.h>
#include <google/cloud/storage/client.h>

#include <optional>

namespace seqr {

namespace gcs = google::cloud::storage;

namespace {

class GcsReader : public UrlReader {
 public:
  explicit GcsReader(const size_t max_connection_pool_size)
      : shared_gcs_client_(
            google::cloud::Options{}.set<gcs::ConnectionPoolSizeOption>(
                max_connection_pool_size)) {}

  absl::StatusOr<std::vector<char>> Read(std::string_view url) const override {
    if (!absl::ConsumePrefix(&url, "gs://")) {
      return absl::InvalidArgumentError(absl::StrCat("Unsupported URL: ", url));
    }

    const size_t slash_pos = url.find_first_of('/');
    if (slash_pos == std::string_view::npos) {
      return absl::InvalidArgumentError(
          absl::StrCat("Incomplete blob URL ", url));
    }

    const auto bucket = url.substr(0, slash_pos);
    const auto blob = url.substr(slash_pos + 1);

    // Make a copy of the GCS client for thread-safety.
    gcs::Client gcs_client = shared_gcs_client_;

    try {
      auto reader =
          gcs_client.ReadObject(std::string(bucket), std::string(blob));
      if (reader.bad()) {
        return absl::InvalidArgumentError(
            absl::StrCat("Failed to read blob: ", reader.status().message()));
      }

      std::optional<int> content_length;
      for (const auto& header : reader.headers()) {
        if (header.first == "content-length") {
          int value = 0;
          if (!absl::SimpleAtoi(header.second, &value)) {
            return absl::NotFoundError(
                "Couldn't parse content-length header value");
          }
          content_length = value;
        }
      }
      if (!content_length) {
        return absl::NotFoundError("Couldn't find content-length header");
      }

      std::vector<char> result(*content_length, '\0');
      reader.read(result.data(), result.size());
      if (reader.bad()) {
        return absl::InvalidArgumentError(
            absl::StrCat("Failed to read blob: ", reader.status().message()));
      }
      return result;
    } catch (const std::exception& e) {
      // Unfortunately the googe-cloud-storage library throws exceptions.
      return absl::InternalError(
          absl::StrCat("Exception during reading of ", url, ": ", e.what()));
    }
  }

 private:
  // Share connection pool, but need to make copies for thread-safety.
  gcs::Client shared_gcs_client_;
};

}  // namespace

absl::StatusOr<std::unique_ptr<UrlReader>> MakeGcsReader(
    const size_t max_connection_pool_size) {
  return std::make_unique<GcsReader>(max_connection_pool_size);
}

}  // namespace seqr
