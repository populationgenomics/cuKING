#!/usr/bin/env python3

"""
Converts a Hail MatrixTable to the Parquet format suitable for cuKING.
It is assumed that the table only contains biallelics, which is the case
for gnomAD QC MTs.
"""

import argparse
import hail as hl
import json


def mt_to_cuking_inputs(mt: hl.MatrixTable, parquet_uri: str, overwrite: bool) -> None:
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
    entry_write = entries.to_spark().write.option('compression', 'zstd')
    if overwrite:
        entry_write = entry_write.mode('overwrite')
    entry_write.parquet(parquet_uri)

    # The Parquet files only contain `col_idx` instead of the sample name `s`, to reduce
    # storage requirements. Save a dict from `col_idx` to `s`, so cuKING can map back to
    # sample names for its results. By collecting both `col_idx` and `s` simultaneously,
    # we don't rely on a particular order being returned by `collect`.
    col_idx_mapping = hl.struct(col_idx=mt.col_idx, s=mt.s).collect()
    col_idx_mapping.sort(key=lambda e: e.col_idx)
    metadata = {
        'num_sites': mt.count_rows(),
        'samples': [e.s for e in col_idx_mapping],
    }
    with hl.utils.hadoop_open(f'{parquet_uri}/metadata.json', 'w') as f:
        json.dump(metadata, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--mt-uri',
        help='Input URI for the Hail MT',
        required=True,
    )
    parser.add_argument(
        '--parquet-uri',
        help='Output URI for the Parquet files',
        required=True,
    )
    parser.add_argument(
        '--overwrite',
        help='Overwrite output files',
        action='store_true',
    )
    args = parser.parse_args()

    hl.init(default_reference='GRCh38')
    mt = hl.read_matrix_table(args.mt_uri)
    mt_to_cuking_inputs(mt=mt, parquet_uri=args.parquet_uri, overwrite=args.overwrite)
