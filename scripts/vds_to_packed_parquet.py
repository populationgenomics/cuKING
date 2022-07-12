#!/usr/bin/env python3

"""Filters a VDS using a given sites table, precomputes fields necessary for KING, and exports the result as Parquet partitions.

gcloud config set project cpg-gnomad-production-27bb && \
gcloud dataproc autoscaling-policies import vds-to-packed-parquet \
    --source=vds_to_packed_parquet_autoscaling_policy.yaml \
    --region=us-central1 && \
hailctl dataproc start \
    --region=us-central1 \
    --autoscaling-policy=vds-to-packed-parquet \
    --service-account=dataproc@cpg-gnomad-production-27bb.iam.gserviceaccount.com \
    --requester-pays-allow-buckets=gnomad,cpg-gnomad-production \
    --num-secondary-workers=1000 \
    --packages gnomad \
    --max-idle=5m \
    vds-to-packed-parquet && \
hailctl dataproc submit --region=us-central1 vds-to-packed-parquet vds_to_packed_parquet.py \
    --input=gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds \
    --sites=gs://cpg-gnomad-production/relatedness/ukb_sites.ht \
    --output=gs://cpg-gnomad-production/relatedness/gnomad_v4.0_king_packed.parquet
"""

import click
from gnomad.utils.annotations import annotate_adj
import hail as hl
import json
import os


@click.command()
@click.option("--input", help="Input path for VDS", required=True)
@click.option("--sites", help="Sites table path for row filtering", required=True)
@click.option("--output", help="Output path for Parquet files", required=True)
@click.option(
    "--remove_duplicate_sample_ids",
    help="Whether to remove duplicate sample IDs. Should not be necessary on production datasets.",
    is_flag=True,
    default=False,
)
@click.option(
    "--head_fraction",
    type=click.FLOAT,
    help="Only take this fraction of the beginning of the dataset, useful for benchmarking only.",
    default=None,
)
def main(input, sites, output, remove_duplicate_sample_ids, head_fraction):
    hl.init(default_reference='GRCh38')

    vds = hl.vds.read_vds(input)

    if remove_duplicate_sample_ids:
        import collections

        duplicate_samples = [
            item
            for item, count in collections.Counter(vds.variant_data.s.collect()).items()
            if count > 1
        ]
        vds.variant_data = vds.variant_data.filter_cols(
            hl.set(duplicate_samples).contains(vds.variant_data.s), keep=False
        )

    if head_fraction:
        vds.variant_data = vds.variant_data.head(
            int(vds.variant_data.count_rows() * head_fraction)
        )
        vds.reference_data = vds.reference_data.head(
            int(vds.reference_data.count_rows() * head_fraction)
        )

    # Convert to purely biallelic sites prior to site filtering. This also takes care of
    # LGT -> GT conversion.
    vds = hl.vds.split_multi(vds)

    # Filter to sites table.
    sites_table = hl.read_table(sites)
    vds = hl.vds.filter_variants(vds, sites_table)

    # Fill in hom-ref calls before applying variant filtering.
    mt = hl.vds.to_dense_mt(vds)

    # Apply basic variant QC.
    mt = annotate_adj(mt)
    mt = mt.filter_entries(mt.adj)

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
    mt = mt.annotate_cols(col_idx_mapping=hl.struct(col_idx=mt.col_idx, s=mt.s))
    col_idx_mapping = mt.col_idx_mapping.collect()
    col_idx_mapping.sort(key=lambda e: e.col_idx)
    metadata = {
        'num_sites': sites_table.count(),
        'samples': [e.s for e in col_idx_mapping],
    }
    with hl.utils.hadoop_open(os.path.join(output, 'metadata.json'), 'w') as f:
        json.dump(metadata, f)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
