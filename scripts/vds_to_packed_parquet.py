#!/usr/bin/env python3

"""Filters a VDS using a given sites table, precomputes fields necessary for KING, and exports the result as Parquet partitions.

hailctl dataproc start --service-account=dataproc-service-account@cpg-gnomad-v4-relatedness.iam.gserviceaccount.com \
    --max-idle=5m --num-secondary-workers=20 --packages gnomad \
    --requester-pays-allow-buckets=sites-for-relatedness-transfer-au-tmp vds-to-packed-parquet && \
hailctl dataproc submit vds-to-packed-parquet vds_to_packed_parquet.py \
    --input=gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test.vds \
    --sites=gs://sites-for-relatedness-transfer-au-tmp/ukb_sites.ht \
    --output=gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_packed.parquet
"""

import click
from gnomad.utils.annotations import annotate_adj
import hail as hl

# str.removesuffix is only available in Python >= 3.9.
def remove_suffix(s: str, suffix: str) -> str:
    return s[:-len(suffix)] if s.endswith(suffix) else s

@click.command()
@click.option("--input", help="Input path for VDS", required=True)
@click.option("--sites", help="Sites table path for row filtering", required=True)
@click.option("--output", help="Output path for Parquet files", required=True)
def main(input, sites, output):
    hl.init(default_reference='GRCh38')

    vds = hl.vds.read_vds(input)

    # Filter samples.
    # TODO: remove this once the duplicate samples have been removed from the input.
    import collections

    duplicate_samples = [
        item
        for item, count in collections.Counter(vds.variant_data.s.collect()).items()
        if count > 1
    ]
    vds.variant_data = vds.variant_data.filter_cols(
        hl.set(duplicate_samples).contains(vds.variant_data.s), keep=False
    )

    # Convert to purely biallelic sites prior to site filtering. This also takes care of LGT -> GT conversion.
    vds = hl.vds.split_multi(vds)

    # Filter to sites table.
    sites_table = hl.read_table(sites)
    vds = hl.vds.filter_variants(vds, sites_table)

    # Fill in hom-ref calls before applying variant filtering.
    mt = hl.vds.to_dense_mt(vds)

    # Apply basic variant QC.
    mt = annotate_adj(mt)
    mt = mt.filter_entries(mt.adj)

    # Compute two-bit representation of fields necessary for KING: hom-ref (0), het (1), hom-var (2), missing (3).
    mt = mt.unfilter_entries()
    mt = mt.select_entries(
        packed_gt=hl.if_else(hl.is_defined(mt.GT), mt.GT.n_alt_alleles(), 3)
    )

    # Create a table based on entries alone.
    mt = mt.select_rows()
    ht = mt.make_table()

    # Remove locus and alleles key, as we only need a row index.
    ht = ht.key_by()
    ht = ht.drop(ht.locus, ht.alleles)

    # Remove the ".packed_gt" suffix from the column names.
    ht = ht.rename({s : remove_suffix(s, '.packed_gt') for s in ht.row})

    # Export to one Parquet file per partition.
    ht.to_spark().write.parquet(output)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
