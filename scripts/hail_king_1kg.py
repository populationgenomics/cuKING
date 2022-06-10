#!/usr/bin/env python3

import hail as hl
from cpg_utils.hail_batch import init_batch, output_path

DENSE_HGDP_1KG_TABLE = 'gs://gcp-public-data--gnomad/release/3.1.2/mt/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_dense.mt'
AF_THRESHOLD = 0.05
TARGET_COUNT = 250000


def main():
    init_batch()

    ht = hl.read_table(DENSE_HGDP_1KG_TABLE)
    ht = ht.filter_cols(ht.hgdp_tgp_meta.project == '1000 Genomes')
    ht = ht.filter(ht.gnomad_freq.AF < AF_THRESHOLD)
    ht = ht.sample(TARGET_COUNT / ht.count())
    kinship = hl.king(ht.GT)
    kinship.write(output_path('king_1kg.mt'))


if __name__ == '__main__':
    main()
