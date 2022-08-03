#!/usr/bin/env python3

"""
Converts a dense MatrixTable to the Parquet format suitable for cuKING.
It is assumed that the table only contains biallelics, which is the case
for gnomAD QC MTs.
"""

import click
import hail as hl
from hail.utils import hadoop_open
import json
import os


@click.command()
@click.option("--input", help="Input path for the dense MT", required=True)
@click.option("--output", help="Output path for Parquet files", required=True)
def main(input, output):
    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(input)

    # Compute the field required for KING.
    mt = mt.select_entries(n_alt_alleles=mt.GT.n_alt_alleles())

    # Create a table based on entries alone. By adding a row and column index, we can
    # drop the locus, alleles, and sample field, as well as skip writing missing
    # entries.
    mt = mt.select_globals()
    mt = mt.select_rows()
    mt = mt.select_cols()
    mt = mt.add_row_index()
    mt = mt.add_col_index()
    entries = mt.entries()
    entries = entries.key_by()
    entries = entries.select(entries.row_idx, entries.col_idx, entries.n_alt_alleles)

    # Export one compressed Parquet file per partition.
    entries.to_spark().write.option('compression', 'zstd').parquet(output)

    # Write metadata that's useful for postprocessing. Map `col_idx` to `s` explicitly
    # so we don't need to rely on a particular order returned by `collect`.
    col_idx_mapping = hl.struct(col_idx=mt.col_idx, s=mt.s).collect()
    col_idx_mapping.sort(key=lambda e: e.col_idx)
    metadata = {
        'num_sites': mt.count_rows(),
        'samples': [e.s for e in col_idx_mapping],
    }
    with hadoop_open(os.path.join(output, 'metadata.json'), 'w') as f:
        json.dump(metadata, f)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
