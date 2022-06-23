// Converts a single-sample GVCF to a compact binary format used by cuKING.
//
// Note: assumes alignment with GRCh38 to convert to global positions.
//
// To use GCS input paths, set the GCS_OAUTH_TOKEN beforehand
// (see https://github.com/samtools/htslib/pull/446):
// export GOOGLE_APPLICATION_CREDENTIALS="/path/to/key.json"
// export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

#include <absl/cleanup/cleanup.h>
#include <absl/container/flat_hash_map.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <cassert>
#include <functional>
#include <iostream>
#include <string>
#include <string_view>

#include "gcs_client.h"

ABSL_FLAG(std::string, input, "",
          "The GVCF input path, e.g. gs://some/bucket/NA12878.g.vcf.gz");
ABSL_FLAG(std::string, output, "",
          "The cuking output path, e.g. gs://another/bucket/NA12878.cuking");
ABSL_FLAG(std::string, sites_table, "",
          "The set of sites (locus + var) to retain; everything else will be "
          "filtered out, e.g. "
          "gs://some/bucket/sites_gnomad_ld_pruned_combined_variants_v2.bin");

namespace {

// Based on boost::hash_combine. Hashes `v` before combining.
template <typename T>
inline size_t HashCombine(const size_t seed, const T& v) {
  return seed ^ (std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
}

// The sites table on disk doesn't use any padding.
#pragma pack(push, 1)
struct SitesTableEntry {
  uint64_t global_position;
  char allele;
};
#pragma pack(pop)

struct SitesTableEntryHash {
  size_t operator()(const SitesTableEntry& e) const {
    return HashCombine(std::hash<uint64_t>()(e.global_position), e.allele);
  }
};

struct SitesTableEntryEqual {
  size_t operator()(const SitesTableEntry& lhs,
                    const SitesTableEntry& rhs) const {
    return lhs.global_position == rhs.global_position &&
           lhs.allele == rhs.allele;
  }
};

inline void SetBit(uint64_t* const bit_set, uint64_t index) {
  bit_set[index >> 6] |= 1ull << (index & 63ull);
}

// Returns ceil(a / b) for integers a, b.
template <typename T>
inline T CeilIntDiv(const T a, const T b) {
  return (a + b - 1) / b;
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

  const auto& sites_table = absl::GetFlag(FLAGS_sites_table);
  if (sites_table.empty()) {
    std::cerr << "Error: no sites table specified." << std::endl;
    return 1;
  }

  // Read the sites table.
  auto gcs_client = cuking::NewGcsClient(1);
  auto sites_table_str = gcs_client->Read(sites_table);
  if (!sites_table_str.ok()) {
    std::cerr << "Error reading \"" << sites_table
              << "\": " << sites_table_str.status() << std::endl;
    return 1;
  }

  const SitesTableEntry* sites_table_array =
      reinterpret_cast<const SitesTableEntry*>(sites_table_str->data());
  const size_t num_sites = sites_table_str->size() / sizeof(SitesTableEntry);
  std::cout << "Read " << num_sites << " sites." << std::endl;

  // Hash each site to its index, which will be used in the bit set.
  absl::flat_hash_map<SitesTableEntry, uint64_t, SitesTableEntryHash,
                      SitesTableEntryEqual>
      index_by_site;
  for (uint64_t i = 0; i < num_sites; ++i) {
    index_by_site[sites_table_array[i]] = i;
  }

  sites_table_str->clear();  // Free memory.

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
  std::vector<uint64_t> bit_sets(CeilIntDiv<size_t>(num_sites, 64) * 2);
  uint64_t* const het_bit_set = bit_sets.data();
  uint64_t* const hom_alt_bit_set = bit_sets.data() + bit_sets.size() / 2;
  uint32_t processed = 0, num_het = 0, num_hom_alt = 0;
  int *genotypes = nullptr, *gq = nullptr, *dp = nullptr;
  const absl::Cleanup free_buffers = [genotypes, gq, dp] {
    free(genotypes);
    free(gq);
    free(dp);
  };

  // Read the GVCF one record at a time.
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

    const int gt[2] = {bcf_gt_allele(genotypes[0]),
                       bcf_gt_allele(genotypes[1])};

    if (gt[0] == 0 && gt[1] == 0) {
      continue;  // Ignore hom-ref.
    }

    if (strlen(record->d.allele[0]) != 1) {
      continue;  // Ignore indels.
    }

    int num_gq = 0, num_dp = 0;
    if (bcf_get_format_int32(hts_header, record, "GQ", &gq, &num_gq) <= 0 ||
        num_gq != 1 ||
        bcf_get_format_int32(hts_header, record, "DP", &dp, &num_dp) <= 0 ||
        num_dp != 1) {
      continue;
    }

    // Filter bad quality calls, inspired by gnomAD criteria:
    // https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#the-gnomad-hgdp-and-1000-genomes-callset
    if (gq[0] < 20 || dp[0] < 10) {
      continue;
    }

    const auto find_index = [&](const uint64_t global_position,
                                const int genotype) {
      const char* const allele = record->d.allele[genotype];
      if (strlen(allele) != 1) {
        return std::optional<uint64_t>();  // Ignore indels.
      }

      const auto index = index_by_site.find({global_position, allele[0]});
      if (index == index_by_site.end()) {
        return std::optional<uint64_t>();  // Filtered out.
      }

      return std::optional<uint64_t>(index->second);
    };

    if (gt[0] == gt[1]) {  // Hom-alt.
      if (const auto index = find_index(global_position, gt[0]); index) {
        SetBit(hom_alt_bit_set, *index);
        ++num_hom_alt;
      }
    } else {  // Heterozygous, potentially with two non-ref alleles.
      for (int k = 0; k < 2; ++k) {
        if (gt[k] == 0) {
          continue;
        }
        if (const auto index = find_index(global_position, gt[k]); index) {
          SetBit(het_bit_set, *index);
          ++num_het;
        }
      }
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
