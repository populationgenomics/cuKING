#!/usr/bin/env python3

"""
Converts cuKING outputs in Parquet format to a Hail table.
"""

import argparse
import hail as hl


def cuking_outputs_to_ht(parquet_uri: str) -> hl.Table:
    spark = hl.utils.java.Env.spark_session()
    df = spark.read.parquet(parquet_uri)
    ht = hl.Table.from_spark(df)
    ht = ht.key_by(ht.i, ht.j)
    return ht


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--parquet-uri',
        help='Input URI for the Parquet files',
        required=True,
    )
    parser.add_argument(
        '--ht-uri',
        help='Output URI for the Hail table',
        required=True,
    )
    parser.add_argument(
        '--overwrite',
        help='Overwrite output files',
        action='store_true',
    )
    args = parser.parse_args()

    hl.init(default_reference='GRCh38')
    ht = cuking_outputs_to_ht(args.parquet_uri)
    ht.write(args.ht_uri, overwrite=args.overwrite)
