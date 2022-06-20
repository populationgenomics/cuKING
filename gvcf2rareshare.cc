#include <absl/cleanup/cleanup.h>
#include <absl/container/flat_hash_map.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <queue>
#include <string>
#include <string_view>

#include "gcs_client.h"

ABSL_FLAG(std::string, input, "",
          "The GVCF input path, e.g. gs://some/bucket/NA12878.g.vcf.gz");
ABSL_FLAG(std::string, output, "",
          "The cuking output path, e.g. gs://another/bucket/NA12878.cuking");
ABSL_FLAG(std::string, af_table, "",
          "The allele frequency table to use, e.g. "
          "gs://some/bucket/gnomad_v3_popmax_af.bin");
ABSL_FLAG(size_t, num_output_variants, 10000,
          "How many rare variants to collect.");

namespace {

#pragma pack(push, 1)
struct AfTableEntry {
  uint64_t global_position;
  char allele;
  float allele_frequency;
};
#pragma pack(pop)

std::optional<float> FindMinAf(const uint64_t global_position,
                               const char allele0, const char allele1,
                               const AfTableEntry* const af_table,
                               const size_t af_table_size) {
  // First use binary search to find the applicable range.
  AfTableEntry search_entry;
  search_entry.global_position = global_position;
  const auto iter_pair =
      std::equal_range(af_table, af_table + af_table_size, search_entry,
                       [](const AfTableEntry& lhs, const AfTableEntry& rhs) {
                         return lhs.global_position < rhs.global_position;
                       });

  // Then iterate over the range to find the minimum allele frequency, given the
  // alleles.
  bool found = false;
  float min_af = 1.f;
  for (auto iter = iter_pair.first; iter != iter_pair.second; ++iter) {
    if (iter->allele == allele0 || iter->allele == allele1) {
      min_af = std::min(min_af, iter->allele_frequency);
      found = true;
    }
  }
  return found ? min_af : std::optional<float>();
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

  const auto& af_table = absl::GetFlag(FLAGS_af_table);
  if (af_table.empty()) {
    std::cerr << "Error: no allele frequency table specified." << std::endl;
    return 1;
  }

  // Read the AF table.
  auto gcs_client = cuking::NewGcsClient(1);
  auto af_table_str = gcs_client->Read(af_table);
  if (!af_table_str.ok()) {
    std::cerr << "Error: " << af_table_str.status() << std::endl;
    return 1;
  }
  const AfTableEntry* af_table_array =
      reinterpret_cast<const AfTableEntry*>(af_table_str->data());
  const size_t af_table_size = af_table_str->size() / sizeof(AfTableEntry);
  std::cout << "Read " << af_table_size << " allele frequency table entries."
            << std::endl;

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

  struct RareVariant {
    uint64_t global_position;
    float allele_frequency;

    bool operator<(const RareVariant& other) const {
      return std::tie(allele_frequency, global_position) <
             std::tie(other.allele_frequency, other.global_position);
    }
  };

  // Use a priority queue to keep the `num_output_variants` with the smallest
  // allele frequencies, while scanning through all variants in the GVCF file.
  const size_t num_output_variants = absl::GetFlag(FLAGS_num_output_variants);
  std::priority_queue<RareVariant> variant_queue;

  uint32_t processed = 0;
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

    if (bcf_unpack(record, BCF_UN_ALL) != 0) {
      std::cerr << "Error: failed to unpack record." << std::endl;
      return 1;
    }

    int num_genotypes = 0;
    if (bcf_get_format_int32(hts_header, record, "GT", &genotypes,
                             &num_genotypes) <= 0 ||
        num_genotypes != 2) {
      continue;  // GT not present or unexpected number of GT entries.
    }

    const int gt0 = bcf_gt_allele(genotypes[0]);
    const int gt1 = bcf_gt_allele(genotypes[1]);

    if (gt0 == 0 && gt1 == 0) {
      continue;  // hom-ref
    }

    char allele0 = 0;
    if (gt0 != 0) {
      const auto allele = record->d.allele[gt0];
      if (strlen(allele) == 1) {
        allele0 = allele[0];
      }
    }

    char allele1 = 0;
    if (gt1 != 0) {
      const auto allele = record->d.allele[gt1];
      if (strlen(allele) == 1) {
        allele1 = allele[0];
      }
    }

    // Find the minimum allele frequency for this locus + genotypes.
    const auto min_af = FindMinAf(global_position, allele0, allele1,
                                  af_table_array, af_table_size);

    if (min_af) {
      if (variant_queue.size() < num_output_variants) {
        variant_queue.push({global_position, *min_af});
      } else if (variant_queue.top().allele_frequency > *min_af) {
        variant_queue.pop();
        variant_queue.push({global_position, *min_af});
      }
    }
  }

  // Encode the output variant loci using deltas, so we only need 32 bits per
  // locus.
  std::vector<uint64_t> global_positions;
  global_positions.reserve(variant_queue.size());
  while (!variant_queue.empty()) {
    global_positions.push_back(variant_queue.top().global_position);
    variant_queue.pop();
  }

  std::sort(global_positions.begin(), global_positions.end());

  std::vector<uint32_t> deltas;
  deltas.reserve(global_positions.size());
  uint64_t last_position = 0;
  for (const uint64_t position : global_positions) {
    const uint64_t delta = position - last_position;
    if (delta >= 1ULL << 32) {
      std::cerr << "Error: can't encode delta at " << position
                << " using 32 bits." << std::endl;
      return 1;
    }
    deltas.push_back(delta);
    last_position = position;
  }

  // Write the output file.
  if (auto status = gcs_client->Write(
          output_file, std::string(reinterpret_cast<char*>(deltas.data()),
                                   deltas.size() * sizeof(uint32_t)));
      !status.ok()) {
    std::cerr << "Error: " << status << std::endl;
    return 1;
  }

  return 0;
}