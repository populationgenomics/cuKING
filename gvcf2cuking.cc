// Converts a single sample GVCF to a compact binary format used by cuKING.
//
// Note: assumes alignment with GRCh38 to convert to global positions.
//
// Example: ./gvcf2cuking --input=NA12878.g.vcf.gz --output=NA12878.cuking

#include <absl/container/flat_hash_map.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <vcflib/Variant.h>

#include <iostream>
#include <string>
#include <string_view>

ABSL_FLAG(std::string, input, "",
          "The GVCF input filename, e.g. NA12878.g.vcf.gz");
ABSL_FLAG(std::string, output, "",
          "The cuking output filename, e.g. NA12878.cuking");

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
  std::string input_file_copy = input_file;  // open() below requires a copy
  vcf.open(input_file_copy);
  if (!vcf.is_open()) {
    std::cerr << "Error: failed to open " << input_file << std::endl;
    return 1;
  }

  // Used to convert to global position (for GRCh38).
  const absl::flat_hash_map<string_view, uint64_t> chr_offsets = {
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
  int count = 0;
  while (vcf.getNextVariant(variant)) {
    const auto offset_iter = chr_offsets.find(variant.sequenceName);
    if (offset_iter == chr_offsets.end()) {
      continue;  // Ignore this contig, e.g. alts.
    }
    // VCF positions are 1-based, so subtract one.
    const uint64_t global_position = offset_iter->second + variant.position - 1;

    std::cout << variant << std::endl;
    std::cout << variant.sequenceName << std::endl;
    std::cout << variant.position << std::endl;
    std::cout << variant.id << std::endl;
    std::cout << variant.ref << std::endl;
    // TODO: Get GT from `samples`.
    if (++count >= 100) {
      break;
    }
  }

  return 0;
}