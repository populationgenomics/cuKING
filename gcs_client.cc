#include "gcs_client.h"

#include <absl/strings/strip.h>
#include <google/cloud/storage/client.h>

#include <optional>

namespace gcs_client {

namespace gcs = google::cloud::storage;

namespace {

// Returns the bucket + blob name from a full gs:// URL.
absl::StatusOr<std::pair<std::string_view, std::string_view>> SplitBlobPath(
    std::string_view url) {
  if (!absl::ConsumePrefix(&url, "gs://")) {
    return absl::InvalidArgumentError(absl::StrCat("Unsupported URL: ", url));
  }

  const size_t slash_pos = url.find_first_of('/');
  if (slash_pos == std::string_view::npos) {
    return absl::InvalidArgumentError(
        absl::StrCat("Incomplete blob URL ", url));
  }

  return std::make_pair(url.substr(0, slash_pos), url.substr(slash_pos + 1));
}

class GcsClientImpl : public GcsClient {
 public:
  explicit GcsClientImpl(const size_t max_connection_pool_size)
      : shared_gcs_client_(
            google::cloud::Options{}.set<gcs::ConnectionPoolSizeOption>(
                max_connection_pool_size)) {}

  absl::StatusOr<std::string> Read(std::string_view url) const override {
    auto bucket_and_object = SplitBlobPath(url);
    if (!bucket_and_object.ok()) {
      return bucket_and_object.status();
    }
    auto [bucket, blob] = *bucket_and_object;

    // Make a copy of the GCS client for thread-safety.
    gcs::Client gcs_client = shared_gcs_client_;

    try {
      auto reader =
          gcs_client.ReadObject(std::string(bucket), std::string(blob));
      if (reader.bad()) {
        return absl::InvalidArgumentError(absl::StrCat(
            "Failed to read blob ", url, ": ", reader.status().message()));
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

      std::ostringstream os;
      os << reader.rdbuf();
      if (reader.bad()) {
        return absl::InvalidArgumentError(absl::StrCat(
            "Failed to read blob ", url, ": ", reader.status().message()));
      }
      return os.str();
    } catch (const std::exception& e) {
      // Unfortunately the googe-cloud-storage library throws exceptions.
      return absl::InternalError(
          absl::StrCat("Exception during reading of ", url, ": ", e.what()));
    }
  }

  absl::Status Write(std::string_view url,
                     std::string contents) const override {
    auto bucket_and_object = SplitBlobPath(url);
    if (!bucket_and_object.ok()) {
      return bucket_and_object.status();
    }
    // Make a copy of the GCS client for thread-safety.
    gcs::Client gcs_client = shared_gcs_client_;

    auto status_or_metadata = gcs_client.InsertObject(
        std::string(bucket_and_object->first),
        std::string(bucket_and_object->second), std::move(contents));

    if (!status_or_metadata.ok()) {
      return absl::InvalidArgumentError(
          absl::StrCat("Failed to write blob ", url, ": ",
                       status_or_metadata.status().message()));
    }

    return absl::OkStatus();
  }

 private:
  // Share connection pool, but need to make copies for thread-safety.
  gcs::Client shared_gcs_client_;
};

}  // namespace

std::unique_ptr<GcsClient> NewGcsClient(const size_t max_connection_pool_size) {
  return std::make_unique<GcsClientImpl>(max_connection_pool_size);
}

}  // namespace gcs_client
