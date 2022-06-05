// Converts a single-sample GVCF to a compact binary format used by cuKING.
//
// Note: assumes alignment with GRCh38 to convert to global positions.
//
// Example: ./gvcf2cuking --input=NA12878.g.vcf.gz --output=NA12878.cuking

#include <absl/cleanup/cleanup.h>
#include <absl/container/flat_hash_map.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <zstd.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

#include "cuking.h"

ABSL_FLAG(std::string, input, "",
          "The GVCF input filename, e.g. NA12878.g.vcf.gz");
ABSL_FLAG(std::string, output, "",
          "The cuking output filename, e.g. NA12878.cuking");

namespace {

constexpr uint16_t kMaxLocusDelta = (1 << 14) - 1;

inline uint16_t Encode(const int64_t locus_delta,
                       const cuking::VariantCategory variant_category) {
  assert(locus_delta <= kMaxLocusDelta);
  return (locus_delta << 2) | static_cast<uint16_t>(variant_category);
}

}  // namespace

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

  htsFile* const hts_file = bcf_open(input_file.c_str(), "r");
  if (hts_file == nullptr) {
    std::cerr << "Error: could not open \"" << input_file << "\"." << std::endl;
    return 1;
  }
  const absl::Cleanup hts_file_closer = [hts_file] { bcf_close(hts_file); };

  bcf_hdr_t* const hts_header = bcf_hdr_read(hts_file);
  if (hts_header == nullptr) {
    std::cerr << "Error: failed to read header." << std::endl;
    return 1;
  }
  const absl::Cleanup hts_header_destroyer = [hts_header] {
    bcf_hdr_destroy(hts_header);
  };

  const int num_samples = bcf_hdr_nsamples(hts_header);
  if (num_samples != 1) {
    std::cerr << "Error: unexpected number of samples: " << num_samples
              << "; only single-sample GVCFs are supported." << std::endl;
    return 1;
  }

  int num_seqs = 0;
  const auto seq_names = bcf_hdr_seqnames(hts_header, &num_seqs);
  if (seq_names == nullptr) {
    std::cerr << "Error: failed to read sequence names." << std::endl;
    return 1;
  }
  const absl::Cleanup seq_names_free = [seq_names] { free(seq_names); };

  // Used to convert to global position (for GRCh38).
  const absl::flat_hash_map<std::string_view, int64_t> chr_offsets = {
      {"chr1", 0},           {"chr2", 248956422},   {"chr3", 491149951},
      {"chr4", 689445510},   {"chr5", 879660065},   {"chr6", 1061198324},
      {"chr7", 1232004303},  {"chr8", 1391350276},  {"chr9", 1536488912},
      {"chr10", 1674883629}, {"chr11", 1808681051}, {"chr12", 1943767673},
      {"chr13", 2077042982}, {"chr14", 2191407310}, {"chr15", 2298451028},
      {"chr16", 2400442217}, {"chr17", 2490780562}, {"chr18", 2574038003},
      {"chr19", 2654411288}, {"chr20", 2713028904}, {"chr21", 2777473071},
      {"chr22", 2824183054}, {"chrX", 2875001522},  {"chrY", 3031042417},
  };

  bcf1_t* const record = bcf_init();
  if (record == nullptr) {
    std::cerr << "Error: failed to init record." << std::endl;
    return 1;
  }
  const absl::Cleanup record_destroyer = [record] { bcf_destroy(record); };

  int64_t last_position = -1, num_skips = 0, processed = 0;
  uint32_t num_hets = 0;
  int* genotypes = nullptr;
  const absl::Cleanup genotypes_free = [genotypes] { free(genotypes); };
  std::vector<uint16_t> encoded;
  while (bcf_read(hts_file, hts_header, record) == 0) {
    if ((++processed & ((1 << 20) - 1)) == 0) {
      std::cout << "Processed " << (processed >> 20) << " Mi records, encoded "
                << encoded.size() << " entries..." << std::endl;
    }

    const auto offset = chr_offsets.find(seq_names[record->rid]);
    if (offset == chr_offsets.end()) {
      continue;  // Ignore this contig, e.g. alts.
    }

    const int64_t global_position = offset->second + record->pos;
    if (global_position <= last_position) {
      std::cerr << "Error: variants not sorted by global position ("
                << seq_names[record->rid] << ", " << record->pos << ")."
                << std::endl;
      return 1;
    }

    int num_genotypes = 0;
    if (bcf_get_format_int32(hts_header, record, "GT", &genotypes,
                             &num_genotypes) <= 0 ||
        num_genotypes != 2) {
      continue;  // GT not present or unexpected number of GT entries.
    }

    const int num_ref_alleles =
        (bcf_gt_allele(genotypes[0]) == 0) + (bcf_gt_allele(genotypes[1]) == 0);

    if (num_ref_alleles == 2) {  // HomRef
      continue;
    }

    int64_t locus_delta = global_position - last_position;
    // If the delta doesn't fit, add skips.
    while (locus_delta > kMaxLocusDelta) {
      encoded.push_back(
          Encode(kMaxLocusDelta - 1, cuking::VariantCategory::kSkip));
      locus_delta -= kMaxLocusDelta - 1;
      ++num_skips;
    }

    if (num_ref_alleles == 1) {
      encoded.push_back(Encode(locus_delta, cuking::VariantCategory::kHet));
      ++num_hets;
    } else if (num_ref_alleles == 0) {
      encoded.push_back(Encode(locus_delta, cuking::VariantCategory::kHomAlt));
    }

    last_position = global_position;
  }

  std::cout << "Stats:" << std::endl
            << "  entries: " << encoded.size() << std::endl
            << "  hets: " << num_hets << std::endl
            << "  skips: " << num_skips << std::endl;

  // zstd-compress the result.
  const size_t encoded_byte_size = encoded.size() * sizeof(uint16_t);
  const size_t compress_bound = ZSTD_compressBound(encoded_byte_size);
  std::vector<uint8_t> compressed(compress_bound);
  constexpr int kZstdCompressionLevel = 3;
  const size_t zstd_result =
      ZSTD_compress(compressed.data(), compressed.size(), encoded.data(),
                    encoded_byte_size, kZstdCompressionLevel);
  if (ZSTD_isError(zstd_result)) {
    std::cerr << "Error: failed to zstd-compress the result." << std::endl;
    return 1;
  }
  compressed.resize(zstd_result);

  // Write the output file.
  std::ofstream out(output_file, std::ios::out | std::ios::binary);
  cuking::FileHeader file_header;
  memcpy(file_header.magic, cuking::kExpectedMagic.data(),
         sizeof(file_header.magic));
  file_header.decompressed_size = encoded_byte_size;
  file_header.num_hets = num_hets;
  out.write(reinterpret_cast<const char*>(&file_header),
            sizeof(cuking::FileHeader));  // Decompressed size.
  out.write(reinterpret_cast<const char*>(compressed.data()),
            compressed.size());  // Compressed data.
  out.close();

  return 0;
}