#!/usr/bin/env python3

"""Runs KING and PC-Relate on the gnomAD v4 test dataset.

hailctl dataproc start --service-account=dataproc-service-account@cpg-gnomad-v4-relatedness.iam.gserviceaccount.com \
    --max-idle=5m --num-secondary-workers=20 --packages gnomad \
    --requester-pays-allow-buckets=sites-for-relatedness-transfer-au-tmp relatedness-comparison && \
hailctl dataproc submit relatedness-comparison relatedness_comparison.py
"""

from gnomad.utils.annotations import annotate_adj
import hail as hl


def write_if_not_exists(table, path):
    if not hl.utils.hadoop_exists(path):
        table.write(path)
    if isinstance(table, hl.table.Table):
        return hl.read_table(path)
    if isinstance(table, hl.matrixtable.MatrixTable):
        return hl.read_matrix_table(path)
    raise ValueError(f'Unexpected table type {type(table)}')


def main():
    hl.init(default_reference='GRCh38')

    vds = hl.vds.read_vds(
        'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test.vds'
    )

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
    sites_table = hl.read_table(
        'gs://sites-for-relatedness-transfer-au-tmp/ukb_sites.ht'
    )
    vds = hl.vds.filter_variants(vds, sites_table)

    # Fill in hom-ref calls before applying variant filtering.
    mt = hl.vds.to_dense_mt(vds)

    # Apply basic variant QC.
    mt = annotate_adj(mt)
    mt = mt.filter_entries(mt.adj)

    # Store GT-only table before running relatedness methods.
    mt = mt.select_entries(mt.GT)
    mt = write_if_not_exists(
        mt, 'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites.mt'
    )

    # Compute KING.
    king = hl.king(mt.GT)
    king = write_if_not_exists(
        king,
        'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites_king.mt',
    )

    # Compute PCA for PC-Relate.
    _, scores, _ = hl.hwe_normalized_pca(mt.GT, k=10, compute_loadings=False)

    scores = write_if_not_exists(
        scores,
        'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites_pca_scores.ht',
    )

    # Compute PC-Relate.
    for min_individual_maf in (0.001, 0.01, 0.05):
        # Typically we'd pass min_kinship=0.05 here, but in order to compare with KING
        # results, we don't filter.
        pc_relate = hl.pc_relate(
            mt.GT,
            min_individual_maf=min_individual_maf,
            scores_expr=scores[mt.col_key].scores,
            block_size=4096,
            statistics="all",
        )
        pc_relate = write_if_not_exists(
            pc_relate,
            f'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites_pc_relate_{min_individual_maf}_based_on_scores.ht',
        )


if __name__ == '__main__':
    main()
