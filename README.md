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

Once the Parquet tables have been created, `cuKING` needs to run on a VM with an NVIDIA A100 GPU. Cloud Batch with a Container-Optimized OS instance would be ideal for this, but at the time of writing the supported NVIDIA drivers aren't recent enough. We therefore schedule an instance manually and use a startup script to install drivers manually and launch the Docker container. This is wrapped in [`run_cuking.sh`](run_cuking.sh):

```sh
PROJECT=cpg-gnomad-production-27bb \
INPUT_URI=gs://cpg-gnomad-production/relatedness/gnomad_v3.1_qc_mt_v2_sites_dense_king_packed.parquet \
OUTPUT_URI=gs://cpg-gnomad-production/relatedness/gnomad_v4.0_king_relatedness.json KING_COEFF_THRESHOLD=0.05 \
LOGGING_OUTPUT=gs://cpg-gnomad-production/relatedness/gnomad_v4.0_king_relatedness.json \
./run_cuking.sh
```

The resulting JSON file contains a sparse dictionary that only contains sample pairs with KING coefficients larger than the specified threshold.
