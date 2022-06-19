#!/usr/bin/env python3

import os
import hailtop.batch as hb
from cpg_utils.hail_batch import get_config, remote_tmpdir, output_path

DOCKER_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/cuking:1eaea078804541c43b7adb1db7fd5124fc3380cf'
LOCI_TABLE_PATH = 'gs://cpg-thousand-genomes-main/cuking/loci_gnomad_ld_pruned_combined_variants.bin'


def main():
    config = get_config()

    service_backend = hb.ServiceBackend(
        billing_project=config['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )

    batch = hb.Batch(backend=service_backend)

    with open('1kg_gvcfs.txt') as f:
        paths = [line.strip() for line in f.readlines()]

    # Process 20 files per job.
    for chunk in [paths[i : i + 20] for i in range(0, len(paths), 20)]:
        job = batch.new_job(chunk[0])
        job.image(DOCKER_IMAGE)
        job.memory('lowmem')
        job.command('set -x')
        for gvcf_path in chunk:
            # Need to refresh the token regularly.
            job.command(
                'export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)'
            )
            basename = os.path.basename(gvcf_path)
            cuking_path = output_path(basename.replace('.g.vcf.gz', '.cuking'))
            job.command(
                f'gvcf2cuking --input={gvcf_path} --output={cuking_path} --loci_table={LOCI_TABLE_PATH}'
            )

    batch.run(wait=False)


if __name__ == '__main__':
    main()
