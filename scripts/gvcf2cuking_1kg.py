#!/usr/bin/env python3

import click
import os
import hailtop.batch as hb
from cpg_utils.hail_batch import get_config, remote_tmpdir, output_path

DOCKER_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/cuking'
SITES_TABLE_PATH = 'gs://cpg-thousand-genomes-main/cuking/sites_gnomad_ld_pruned_combined_variants_v2.bin'


@click.command()
@click.option('--image-version', help='Docker image version to use', required=True)
def main(image_version):
    config = get_config()

    service_backend = hb.ServiceBackend(
        billing_project=config['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )

    batch = hb.Batch('gvcf2cuking', backend=service_backend)

    with open('1kg_gvcfs.txt') as f:
        paths = [line.strip() for line in f.readlines()]

    # Process 20 files per job.
    for chunk in [paths[i : i + 20] for i in range(0, len(paths), 20)]:
        job = batch.new_job(chunk[0])
        job.image(f'{DOCKER_IMAGE}:{image_version}')
        job.memory('lowmem')
        for gvcf_path in chunk:
            # Need to refresh the token regularly.
            job.command(
                'export GCS_OAUTH_TOKEN=$(print_google_service_account_access_token.py)'
            )
            basename = os.path.basename(gvcf_path)
            cuking_path = output_path(basename.replace('.g.vcf.gz', '.cuking'))
            # Wrap in "set -x" to print command, but not put access token in logs.
            job.command(
                f'set -x; gvcf2cuking --input={gvcf_path} --output={cuking_path} --sites_table={SITES_TABLE_PATH}; set +x'
            )

    batch.run(wait=False)


if __name__ == '__main__':
    main()
