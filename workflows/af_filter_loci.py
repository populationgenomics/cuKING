#!/usr/bin/env python3

from array import array
from cloudpathlib import AnyPath
import hail as hl
from cpg_utils.hail_batch import init_batch, output_path

GNOMAD_V3_TABLE = 'gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht'
AF_THRESHOLD = 0.05


def main():
    init_batch(driver_cores=8)

    ht = hl.read_table(GNOMAD_V3_TABLE)

    # Filter rare variants.
    ht = ht.filter(ht.popmax.AF > AF_THRESHOLD)

    # Convert to global position and collect list.
    ht = ht.select(pos=ht.locus.global_position())
    global_positions = ht.pos.collect()
    print(f'Filtered to {len(global_positions)} loci.')

    # Write to binary output file.
    path = output_path(f'loci_popmax_af_gt_{str(AF_THRESHOLD)}.bin')
    with AnyPath(path).open('wb') as f:
        array('L', global_positions).tofile(f)  # 'L' == uint64_t


if __name__ == '__main__':
    main()
