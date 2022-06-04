// Converts a single-sample GVCF to a compact binary format used by cuKING.
//
// Note: assumes alignment with GRCh38 to convert to global positions.
//
// Example: ./gvcf2cuking --input=NA12878.g.vcf.gz --output=NA12878.cuking

#include <absl/container/flat_hash_map.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <vcflib/Variant.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

ABSL_FLAG(std::string, input, "",
          "The GVCF input filename, e.g. NA12878.g.vcf.gz");
ABSL_FLAG(std::string, output, "",
          "The cuking output filename, e.g. NA12878.cuking");

enum class VariantCategory {
  kHet = 0,
  kHomAlt = 1,
  kSkip = 2,
};

constexpr uint16_t kMaxLocusDelta = (1 << 14) - 1;

inline uint16_t encode(const int64_t locus_delta,
                       const VariantCategory variant_category) {
  assert(locus_delta <= kMaxLocusDelta);
  return (locus_delta << 2) | static_cast<uint16_t>(variant_category);
}

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  const auto& input_file = absl::GetFlag(FLAGS_input);
  if (input_file.empty()) {
    std::cerr << "Error: no input file specified." << std::endl;
    return 1;
  }

  const auto& output_file = absl::GetFlag(FLAGS_output);
  if (output_file.empty()) {
    std::cerr << "Error: no output file specified." << std::endl;
    return 1;
  }

  vcflib::VariantCallFile vcf;
  std::string input_file_copy = input_file;  // open() below requires a copy.
  vcf.open(input_file_copy);
  if (!vcf.is_open()) {
    std::cerr << "Error: failed to open \"" << input_file << "\"." << std::endl;
    return 1;
  }

  // Used to convert to global position (for GRCh38).
  const absl::flat_hash_map<string_view, int64_t> chr_offsets = {
      {"chr1", 0},           {"chr2", 248956422},   {"chr3", 491149951},
      {"chr4", 689445510},   {"chr5", 879660065},   {"chr6", 1061198324},
      {"chr7", 1232004303},  {"chr8", 1391350276},  {"chr9", 1536488912},
      {"chr10", 1674883629}, {"chr11", 1808681051}, {"chr12", 1943767673},
      {"chr13", 2077042982}, {"chr14", 2191407310}, {"chr15", 2298451028},
      {"chr16", 2400442217}, {"chr17", 2490780562}, {"chr18", 2574038003},
      {"chr19", 2654411288}, {"chr20", 2713028904}, {"chr21", 2777473071},
      {"chr22", 282418305},  {"chrX", 287500152},   {"chrY", 303104241},
  };

  vcflib::Variant variant(vcf);
  int64_t last_position = -1, num_skips = 0, processed = 0;
  std::vector<uint16_t> encoded;
  while (vcf.getNextVariant(variant)) {
    if ((++processed & ((1 << 20) - 1)) == 0) {
      std::cout << "Processed " << (processed >> 20) << " Mi records..."
                << std::endl;
    }

    if (variant.sampleNames.size() != 1) {
      std::cerr << "Error: unexpected number of samples: "
                << variant.sampleNames.size()
                << "; only single-sample GVCFs are supported." << std::endl
                << variant << std::endl;
      return 1;
    }

    const auto offset = chr_offsets.find(variant.sequenceName);
    if (offset == chr_offsets.end()) {
      continue;  // Ignore this contig, e.g. alts.
    }
    const int64_t global_position =
        offset->second + variant.zeroBasedPosition();

    if (global_position <= last_position) {
      std::cerr << "Error: variants not sorted by global position." << std::endl
                << variant << std::endl;
      return 1;
    }

    int64_t locus_delta = global_position - last_position;
    // If the delta doesn't fit, add skips.
    while (locus_delta > kMaxLocusDelta) {
      encoded.push_back(encode(kMaxLocusDelta - 1, VariantCategory::kSkip));
      locus_delta -= kMaxLocusDelta - 1;
      ++num_skips;
    }

    const std::string genotype =
        variant.getGenotype(variant.sampleNames.front());
    const auto decomposed_gt = vcflib::decomposeGenotype(genotype);

    if (vcflib::isHet(decomposed_gt)) {
      encoded.push_back(encode(locus_delta, VariantCategory::kHet));
    } else if (vcflib::isHomNonRef(decomposed_gt)) {
      encoded.push_back(encode(locus_delta, VariantCategory::kHomAlt));
    }

    last_position = global_position;
  }

  std::ofstream out(output_file, std::ios::out | std::ios::binary);
  out.write(reinterpret_cast<const char*>(encoded.data()),
            encoded.size() * sizeof(uint16_t));
  out.close();

  std::cout << "Stats:" << std::endl
            << "  entries: " << encoded.size() << std::endl
            << "  skips: " << num_skips << std::endl;

  return 0;
}