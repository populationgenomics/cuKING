#!/usr/bin/env python3

import os
import hailtop.batch as hb
from cpg_utils.hail_batch import get_config, remote_tmpdir, output_path

DOCKER_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/cuking:d4877256d0b6e3a46465390e249106832e2c90af'


def main():
    config = get_config()

    service_backend = hb.ServiceBackend(
        billing_project=config['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )

    batch = hb.Batch(backend=service_backend)

    with open('1kg_gvcfs.txt') as f:
        for gvcf_path in f.readlines():
            basename = os.path.basename(gvcf_path)
            job = batch.new_job(basename)
            job.image(DOCKER_IMAGE)
            job.memory('lowmem')
            job.command(
                'export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)'
            )
            cuking_path = output_path(basename.replace('.g.vcf.gz', '.cuking'))
            job.command(f'gvcf2cuking --input={gvcf_path} --output={cuking_path}')

    batch.run(wait=False)


if __name__ == '__main__':
    main()
