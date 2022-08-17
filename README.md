# cuKING

`cuKING` is a CUDA-based [KING relatedness](https://www.kingrelatedness.com) estimator that scales to large sample numbers. The resulting kinship coefficients are identical to [Hail's `hl.king` implementation](hail.is/docs/0.2/methods/relatedness.html#hail.methods.king).

Hail's MatrixTable data layout is not well suited for the KING algorithm: KING compares genotypes pairwise, but a MatrixTable stores the genotypes of *all samples* at a locus in a block of memory. Effectively, `cuKING` therefore transposes the MatrixTable memory representation, enabling an efficient bit set based approach similar to [`somalier`](https://github.com/brentp/somalier), but on the GPU.

Also inspired by `somalier`, IBS0, IBS1, and IBS2 values are computed as well. This can be helpful for pedigree inference. Note that `somalier` computes the "within-family" estimator (scaled by a factor of two), while `cuKING` matches Hail's "between-family" estimator.

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

These Parquet tables can be efficiently generated by Hail. Use the [`mt_to_cuking_inputs.py`](mt_to_cuking_inputs.py) script to annotate and convert a Hail MatrixTable.

Below, `$PROJECT_ID` refers to your GCP project ID.

```sh
gcloud config set project PROJECT_ID && \
hailctl dataproc start \
    --region=us-central1 \
    --service-account=dataproc@PROJECT_ID.iam.gserviceaccount.com \
    --requester-pays-allow-buckets=gnomad,cpg-gnomad-production \
    --num-secondary-workers=20 \
    --max-idle=5m \
    mt-to-cuking-inputs && \
hailctl dataproc submit --region=us-central1 mt-to-cuking-inputs mt_to_cuking_inputs.py \
    --mt-uri=gs://cpg-gnomad-production/relatedness/gnomad.exomes.v4.0.pre_ld_prune_qc_sites.bases_over_dp_20_filtered.mt \
    --parquet-uri=gs://cpg-gnomad-production/relatedness/gnomad.exomes.v4.0.pre_ld_prune_qc_sites.parquet
```

In the meantime, build the Docker image for `cuKING` using Cloud Build, which assumes that you have an Artifact Registry Docker repository called `images`:

```sh
gcloud builds submit --region=us-central1 --config cloudbuild.yaml --substitutions=TAG_NAME=$(git describe --tags) .
```

Once the Parquet tables have been created, `cuKING` needs to run on a VM with an NVIDIA A100 GPU, which can be scheduled using Cloud Batch. A Container-Optimized OS instance would be ideal for this, but at the time of writing the preview version of Cloud Batch doesn't support COS yet. Instead, create an Ubuntu-based virtual machine instance template, which includes a [startup script](instance_startup_script.sh) to install recent NVIDIA drivers:

```sh
gcloud compute instance-templates create cuking-instance-template \
    --project=PROJECT_ID \
    --machine-type=a2-highgpu-1g \
    --network-interface=network=default,network-tier=PREMIUM,address="" \
    --metadata-from-file=startup-script=instance_startup_script.sh \
    --maintenance-policy=TERMINATE \
    --provisioning-model=STANDARD \
    --service-account=cuking@PROJECT_ID.iam.gserviceaccount.com \
    --scopes=https://www.googleapis.com/auth/cloud-platform \
    --accelerator=count=1,type=nvidia-tesla-a100 \
    --create-disk=auto-delete=yes,boot=yes,device-name=cuking-instance-template,image=projects/ubuntu-os-cloud/global/images/ubuntu-1804-bionic-v20220712,mode=rw,size=10,type=pd-balanced \
    --no-shielded-secure-boot \
    --shielded-vtpm \
    --shielded-integrity-monitoring \
    --reservation-affinity=any
```

You can use Cloud Batch to run multiple instances of `cuKING` currently. Run [`cloud_batch_submit.py`](cloud_batch_submit.py) to generate and launch a Cloud Batch job:

```sh
./cloud_batch_submit.py \
    --location=us-central1 \
    --project-id=$PROJECT_ID \
    --tag-name=$(git describe --tags) \
    --input-uri=gs://cpg-gnomad-production/relatedness/gnomad.exomes.v4.0.pre_ld_prune_qc_sites.parquet \
    --output-uri=gs://cpg-gnomad-production/relatedness/gnomad.exomes.v4.0.pre_ld_prune_qc_sites_king.parquet \
    --requester-pays-project=$PROJECT_ID \
    --kin-threshold=0.05 \
    --split-factor=4
```

If there are any errors, they'll show up in Cloud Logging for the Cloud Batch job that the above script prints.

## Sharding

If the number of samples and sites is so large that they won't fit into the memory of a single GPU (40 GB for `a2-highgpu-1g` machines), the computation can be sharded. Sharding works by splitting the full relatedness matrix into submatrices that are computed independently, so the results can be easily combined afterwards.

For example, to halve memory requirements, the full matrix can be split into $4 \cdot 4 = 16$ equally sized submatrices (i.e. a "split factor" of 4). Only the "upper triangular" submatrices need to be evaluated due to symmetry of relatedness, leading to 10 shards.

Sharding is implemented through two parameters, `--split_factor` ($k$) and `--shard_index` ($i$), with $0 \leq i < \frac{k(k + 1)}{2}$. [`cloud_batch_submit.py`](cloud_batch_submit.py) sets these parameters accordingly.

Each shard corresponds to a separate Parquet output partition, so results can easily be combined afterwards.

## Downstream analysis

To convert the resulting Parquet files to a Hail table for downstream analysis, use the [`cuking_outputs_to_ht.py`](cuking_outputs_to_ht.py) script. Hail's [`maximal_independent_set`](https://hail.is/docs/0.2/methods/misc.html#hail.methods.maximal_independent_set) method is helpful to prune related samples.
