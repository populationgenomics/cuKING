// Converts a single-sample GVCF to a compact binary format used by cuKING.
//
// Note: assumes alignment with GRCh38 to convert to global positions.
//
// To use GCS input paths, set the GCS_OAUTH_TOKEN beforehand
// (see https://github.com/samtools/htslib/pull/446):
// export GOOGLE_APPLICATION_CREDENTIALS="/path/to/key.json"
// export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
//
// Example:
//   ./gvcf2cuking --input=gs://some/bucket/NA12878.g.vcf.gz
//                 --output=gs://another/bucket/NA12878.cuking

#include <absl/cleanup/cleanup.h>
#include <absl/container/flat_hash_map.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <cassert>
#include <iostream>
#include <string>
#include <string_view>

#include "gcs_client.h"

ABSL_FLAG(std::string, input, "",
          "The GVCF input path, e.g. gs://some/bucket/NA12878.g.vcf.gz");
ABSL_FLAG(std::string, output, "",
          "The cuking output path, e.g. gs://another/bucket/NA12878.cuking");
ABSL_FLAG(std::string, loci_table, "",
          "The set of loci to retain; everything else will be filtered, e.g. "
          "gs://some/bucket/loci_popmax_af_gt_0.05.bin");

namespace {

void SetBit(uint64_t* const bit_set, uint64_t index) {
  bit_set[index >> 6] |= 1ull << (index & 63ull);
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

  const auto& loci_table = absl::GetFlag(FLAGS_loci_table);
  if (loci_table.empty()) {
    std::cerr << "Error: no locus table." << std::endl;
    return 1;
  }

  // Read the loci table.
  auto gcs_client = cuking::NewGcsClient(1);
  auto loci_str = gcs_client->Read(loci_table);
  if (!loci_str.ok()) {
    std::cerr << "Error: " << loci_str.status() << std::endl;
    return 1;
  }
  const uint64_t* loci_array =
      reinterpret_cast<const uint64_t*>(loci_str->data());
  const size_t num_loci = loci_str->size() / sizeof(uint64_t);
  std::cout << "Read " << num_loci << " loci." << std::endl;

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
  const absl::flat_hash_map<std::string_view, uint64_t> chr_offsets = {
      {"chr1", 0ull},           {"chr2", 248956422ull},
      {"chr3", 491149951ull},   {"chr4", 689445510ull},
      {"chr5", 879660065ull},   {"chr6", 1061198324ull},
      {"chr7", 1232004303ull},  {"chr8", 1391350276ull},
      {"chr9", 1536488912ull},  {"chr10", 1674883629ull},
      {"chr11", 1808681051ull}, {"chr12", 1943767673ull},
      {"chr13", 2077042982ull}, {"chr14", 2191407310ull},
      {"chr15", 2298451028ull}, {"chr16", 2400442217ull},
      {"chr17", 2490780562ull}, {"chr18", 2574038003ull},
      {"chr19", 2654411288ull}, {"chr20", 2713028904ull},
      {"chr21", 2777473071ull}, {"chr22", 2824183054ull},
      {"chrX", 2875001522ull},  {"chrY", 3031042417ull},
  };

  bcf1_t* const record = bcf_init();
  if (record == nullptr) {
    std::cerr << "Error: failed to init record." << std::endl;
    return 1;
  }
  const absl::Cleanup record_destroyer = [record] { bcf_destroy(record); };

  // Allocate space for two bit sets.
  std::vector<uint64_t> bit_sets((num_loci + 64 - 1) / 64 * 2);
  uint64_t* const het_bit_set = bit_sets.data();
  uint64_t* const hom_alt_bit_set = bit_sets.data() + bit_sets.size() / 2;
  uint64_t last_position = 0, locus_index = 0;
  uint32_t processed = 0, num_het = 0, num_hom_alt = 0;
  int* genotypes = nullptr;
  const absl::Cleanup genotypes_free = [genotypes] { free(genotypes); };
  while (bcf_read(hts_file, hts_header, record) == 0) {
    if ((++processed & ((1 << 20) - 1)) == 0) {
      std::cout << "Processed " << (processed >> 20) << " Mi records..."
                << std::endl;
    }

    const auto offset = chr_offsets.find(seq_names[record->rid]);
    if (offset == chr_offsets.end()) {
      continue;  // Ignore this contig, e.g. alts.
    }

    const uint64_t global_position = offset->second + record->pos;
    if (global_position < last_position) {
      std::cerr << "Error: variants not sorted by global position ("
                << seq_names[record->rid] << ", " << record->pos << ")."
                << std::endl;
      return 1;
    }
    last_position = global_position;

    // Walk the loci table to check if it contains the current position.
    while (locus_index < num_loci &&
           loci_array[locus_index] < global_position) {
      ++locus_index;
    }
    if (locus_index >= num_loci) {
      break;  // Reached the end of the table.
    }
    if (loci_array[locus_index] > global_position) {
      continue;  // Locus not contained.
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

    if (num_ref_alleles == 1) {
      SetBit(het_bit_set, locus_index);
      ++num_het;
    } else if (num_ref_alleles == 0) {
      SetBit(hom_alt_bit_set, locus_index);
      ++num_hom_alt;
    }
  }

  std::cout << "Stats:" << std::endl
            << "  het: " << num_het << std::endl
            << "  hom_alt: " << num_hom_alt << std::endl;

  // Write the output file.
  if (auto status = gcs_client->Write(
          output_file, std::string(reinterpret_cast<char*>(bit_sets.data()),
                                   bit_sets.size() * sizeof(uint64_t)));
      !status.ok()) {
    std::cerr << "Error: " << status << std::endl;
    return 1;
  }

  return 0;
}