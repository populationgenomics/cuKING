// Converts a single sample GVCF to a compact binary format used by cuKING.

#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <vcflib/Variant.h>

#include <iostream>
#include <string>

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

  vcflib::Variant variant(vcf);
  int count = 0;
  while (vcf.getNextVariant(variant)) {
    std::cout << variant << std::endl;
    if (++count >= 100) {
      break;
    }
  }

  return 0;
}