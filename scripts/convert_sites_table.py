#!/usr/bin/env python3

"""Converts a Hail sites table to packed binary format."""

import click
import struct
import hail as hl
from cloudpathlib import AnyPath


@click.command()
@click.option(
    '--input', help='Input path for Hail sites table', required=True
)
@click.option('--output', help='Output path for packed binary file', required=True)
def main(input, output):
    """Script entry point."""

    hl.init(default_reference='GRCh38')

    ht = hl.read_table(input)

    # Filter to biallelic SNVs, assuming that multi-allelics have been split beforehand.
    ht = ht.filter((hl.len(ht.alleles) == 2) & (hl.len(ht.alleles[0]) == 1) & (hl.len(ht.alleles[1]) == 1))

    # Only retain global positon and the alt allele.
    ht = ht.select(global_position=ht.locus.global_position(), alt=ht.alleles[1])
    ht = ht.key_by()
    ht = ht.select(ht.global_position, ht.alt)

    # Write to output file.
    with AnyPath(output).open('wb') as f:
        for val in ht.collect():
            f.write(struct.pack('=Qc', val.global_position, str.encode(val.alt)))


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
