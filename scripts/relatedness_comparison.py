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

    # Store GT-only table.
    mt = mt.select_entries(mt.GT)
    mt = write_if_not_exists(
        mt, 'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites.mt'
    )

    # LD-prune.
    mt = mt.unfilter_entries()
    pruned_ht = hl.ld_prune(mt.GT, r2=0.1)
    mt = mt.filter_rows(hl.is_defined(pruned_ht[mt.row_key]))
    mt = write_if_not_exists(
        mt, 'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites_pruned.mt'
    )

    # Compute KING.
    king = hl.king(mt.GT)
    king = write_if_not_exists(
        king,
        'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites_pruned_king.mt',
    )

    # Compute PCA for PC-Relate.
    _, scores, loadings = hl.hwe_normalized_pca(mt.GT, k=10, compute_loadings=True)
    scores = write_if_not_exists(
        scores,
        'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites_pruned_pca_scores.ht',
    )
    loadings = write_if_not_exists(
        loadings,
        'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites_pruned_pca_loadings.ht',
    )

    def pc_relate(mt, scores_expr, file_name_suffix):
        """Runs PC-Relate for a few different min_individual_maf values."""
        for min_individual_maf in (0.001, 0.01, 0.05):
            # Typically we'd pass min_kinship=0.05 here, but in order to compare with KING
            # results, we don't filter.
            ht = hl.pc_relate(
                mt.GT,
                min_individual_maf=min_individual_maf,
                scores_expr=scores_expr,
                block_size=4096,
                statistics="all",
            )
            write_if_not_exists(
                ht,
                f'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites_pruned_pc_relate_{min_individual_maf}_{file_name_suffix}.ht',
            )

    # Compute PC-Relate for the new scores.
    pc_relate(mt, scores[mt.col_key].scores, 'new_scores')

    # The previous scores only contain UKB samples, therefore subset samples first.
    previous_scores = hl.read_table('gs://sites-for-relatedness-transfer-au-tmp/pruned.ukb_pca_scores.ht')
    mt = mt.filter_cols(hl.is_defined(previous_scores[mt.s]))
    mt = write_if_not_exists(
        mt, 'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites_pruned_ukb_only.mt'
    )
    pc_relate(mt, previous_scores[mt.col_key].scores, 'previous_scores')

    # Compute new scores only for the UKB sample subset and recompute PC-Relate.
    _, scores, loadings = hl.hwe_normalized_pca(mt.GT, k=10, compute_loadings=True)
    scores = write_if_not_exists(
        scores,
        'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites_pruned_pca_scores_ukb_only.ht',
    )
    loadings = write_if_not_exists(
        loadings,
        'gs://sites-for-relatedness-transfer-au-tmp/gnomad_v4.0_test_ukb_sites_pruned_pca_loadings_ukb_only.ht',
    )
    pc_relate(mt, scores[mt.col_key].scores, 'new_scores_ukb_only')


if __name__ == '__main__':
    main()
