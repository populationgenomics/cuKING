# cuKING

`cuKING` is a CUDA-based [KING relatedness](https://www.chen.kingrelatedness.com/publications/pdf/BI26_2867.pdf) estimator that scales to large sample numbers. The resulting coefficients are identical to [Hail's `hl.king` implementation](hail.is/docs/0.2/methods/relatedness.html#hail.methods.king).

Hail's MatrixTable data layout is not well suited for the KING algorithm: KING compares genotypes pairwise, but a MatrixTable stores the genotypes of *all samples* at a locus in a block of memory. Effectively, `cuKING` therefore transposes the MatrixTable memory representation, allowing an efficient bit set based approach similar to the [C++ KING implementation](https://www.kingrelatedness.com/) or [`somalier`](https://github.com/brentp/somalier), but on the GPU.

## Usage

`cuKING` reads Parquet tables that only contain the minimum necessary information to compute KING coefficients:

```text
+-----------+-----------+-----------------+
|   row_idx |   col_idx |   n_alt_alleles |
|-----------+-----------+-----------------|
|     65360 |       288 |               2 |
|     65360 |       290 |               1 |
|     65360 |       291 |               0 |
|     65360 |       293 |               1 |
|     65360 |       294 |               0 |
|     65360 |       303 |               0 |
|     65360 |       304 |               0 |
|     65360 |       305 |               1 |
|     65360 |       306 |               0 |
|     65360 |       307 |               1 |
|     65360 |       308 |               2 |
|     65360 |       309 |               2 |
```

`row_idx` corresponds to the index of a genomic site, `col_idx` corresponds to a sample index, and `n_alt_alleles` is the number of non-ref alleles at that site.

These tables can be efficiently generated by Hail. When starting from a [Hail VDS](https://hail.is/docs/0.2/vds/index.html#the-data-model-of-variantdataset), use the [`vds_to_packed_parquet.py`](vds_to_packed_parquet.py) script to generate these tables on an autoscaling Dataproc cluster:

```sh
gcloud config set project cpg-gnomad-production-27bb && \
gcloud dataproc autoscaling-policies import vds-to-packed-parquet \
    --source=vds_to_packed_parquet_autoscaling_policy.yaml \
    --region=us-central1 && \
hailctl dataproc start \
    --region=us-central1 \
    --autoscaling-policy=vds-to-packed-parquet \
    --service-account=dataproc@cpg-gnomad-production-27bb.iam.gserviceaccount.com \
    --requester-pays-allow-buckets=gnomad,cpg-gnomad-production \
    --num-secondary-workers=100 \
    --packages gnomad \
    --max-idle=5m \
    vds-to-packed-parquet && \
hailctl dataproc submit --region=us-central1 vds-to-packed-parquet vds_to_packed_parquet.py \
    --input=gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds \
    --sites=gs://cpg-gnomad-production/relatedness/ukb_sites.ht \
    --output=gs://cpg-gnomad-production/relatedness/gnomad_v4.0_king_packed.parquet
```

In the meantime, build the Docker image for `cuKING` using Cloud Build:

```sh
gcloud builds submit --config cloudbuild.yaml .
```

Once the Parquet tables have been created, `cuKING` needs to run on a VM with an NVIDIA A100 GPU, which can be scheduled using Cloud Batch. A Container-Optimized OS instance would be ideal for this, but at the time of writing the preview version of Cloud Batch doesn't support COS yet. Instead, create an Ubuntu-based virtual machine instance template, which includes a [startup script](instance_startup_script.sh) to install recent NVIDIA drivers:

```sh
gcloud compute instance-templates create cuking-instance-template \
    --project=cpg-gnomad-production-27bb \
    --machine-type=a2-highgpu-1g \
    --network-interface=network=default,network-tier=PREMIUM,address="" \
    --metadata-from-file=startup-script=instance_startup_script.sh \
    --no-restart-on-failure \
    --maintenance-policy=TERMINATE \
    --provisioning-model=SPOT \
    --instance-termination-action=DELETE \
    --service-account=cuking@cpg-gnomad-production-27bb.iam.gserviceaccount.com \
    --scopes=https://www.googleapis.com/auth/cloud-platform \
    --accelerator=count=1,type=nvidia-tesla-a100 \
    --create-disk=auto-delete=yes,boot=yes,device-name=cuking-instance-template,image=projects/ubuntu-os-cloud/global/images/ubuntu-1804-bionic-v20220712,mode=rw,size=10,type=pd-balanced \
    --no-shielded-secure-boot \
    --shielded-vtpm \
    --shielded-integrity-monitoring \
    --reservation-affinity=any
```

To launch `cuKING`, submit a Cloud Batch job, e.g.:

```sh
gcloud beta batch jobs submit cuking-gnomad-v4 \
    --location=us-central1 \
    --config=batch_job.json \
    --script-text="sudo docker run --name cuking --gpus all
    us-central1-docker.pkg.dev/cpg-gnomad-production-27bb/images/cuking:latest cuking
    --input_uri=gs://cpg-gnomad-production/relatedness/gnomad_v4.0_king_packed.parquet
    --output_uri=gs://cpg-gnomad-production/relatedness/gnomad_v4.0_king_relatedness.parquet
    --king_coeff_threshold=0.05
    --requester_pays_project=cpg-gnomad-production-27bb"
```

Run the following to check the status of the job:

```sh
gcloud beta batch jobs describe cuking-gnomad-v4
```

If there are any errors, they'll show up in Cloud Logging for that particular Cloud Batch job.

## Sharding

If the number of samples and sites is so large that they won't fit into the memory of a single GPU (40 GB for `a2-highgpu-1g` machines), the computation can be sharded. Sharding works by splitting the full relatedness matrix into submatrices that are computed independently, so the results can be easily combined afterwards.

For example, to halve memory requirements, the full matrix can be split into $4 \cdot 4 = 16$ equally sized submatrices (i.e. a "split factor" of 4). Only the "upper triangular" submatrices need to be evaluated due to symmetry of relatedness, leading to 10 shards.

Sharding is implemented through two parameters, `--split_factor` ($k$) and `--shard_index` ($i$), with $0 \leq i < \frac{k(k + 1)}{2}$.

Each shard corresponds to a separate Parquet output partition, so results can easily be combined afterwards.

## Downstream analysis

To import the results in Hail and prune samples based on relatedness, run the following:

```python
import hail as hl
from hail.utils.java import Env

hl.init(default_reference='GRCh38')

spark = Env.spark_session()
df = spark.read.parquet('gs://cpg-gnomad-production/relatedness/gnomad_v4.0_king_relatedness.parquet')
ht = hl.Table.from_spark(df)

pairs = ht.filter(ht.phi > 0.1)
related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)
```