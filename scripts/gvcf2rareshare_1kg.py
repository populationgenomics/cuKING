#!/usr/bin/env python3

import os
import hailtop.batch as hb
from cpg_utils.hail_batch import get_config, remote_tmpdir, output_path

# TODO: update SHA after fixing table.
DOCKER_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/cuking:e2e5c746cc5791769ffec482150c2e468ee6a5eb'
AF_TABLE_PATH = 'gs://cpg-thousand-genomes-main/cuking/gnomad_v3_popmax_af.bin'


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
        job.memory('16Gi')
        job.command('set -x')
        for gvcf_path in chunk:
            # Need to refresh the token regularly.
            job.command(
                'export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)'
            )
            basename = os.path.basename(gvcf_path)
            rareshare_path = output_path(basename.replace('.g.vcf.gz', '.rareshare'))
            # TODO: fix command after fixing Dockerfile
            job.command(
                f'/app/build/gvcf2rareshare --input={gvcf_path} --output={rareshare_path} --af_table={AF_TABLE_PATH}'
            )

    batch.run(wait=False)


if __name__ == '__main__':
    main()
