#pragma once

#include <string_view>

namespace cuking {

enum class VariantCategory {
  kHet = 0,
  kHomAlt = 1,
};

inline constexpr std::string_view kExpectedMagic = "CUK1";

struct FileHeader {
  uint8_t magic[kExpectedMagic.size()];
  uint64_t decompressed_size;
  uint32_t num_hets;
};

}  // namespace cuking