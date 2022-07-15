#!/usr/bin/env bash

set -ex

ZONE=us-central1-a
INSTANCE=cuking-instance

# Replace parameters for invocation based on environment variables.
envsubst < instance_startup_script.sh > instance_startup_script.sh.out

# Start the instance.
gcloud compute instances create ${INSTANCE} \
    --project=cpg-gnomad-production-27bb \
    --zone=us-central1-a \
    --machine-type=a2-highgpu-1g \
    --metadata-from-file=startup-script=instance_startup_script.sh.out \
    --network-interface=network-tier=PREMIUM,subnet=default \
    --maintenance-policy=TERMINATE \
    --provisioning-model=STANDARD \
    --service-account=dataproc@cpg-gnomad-production-27bb.iam.gserviceaccount.com \
    --scopes=https://www.googleapis.com/auth/cloud-platform \
    --accelerator=count=1,type=nvidia-tesla-a100 \
    --create-disk=auto-delete=yes,boot=yes,device-name=cuking-instance,image=projects/ubuntu-os-cloud/global/images/ubuntu-2204-jammy-v20220712a,mode=rw,size=10,type=projects/cpg-gnomad-production-27bb/zones/us-central1-a/diskTypes/pd-balanced \
    --no-shielded-secure-boot \
    --shielded-vtpm \
    --shielded-integrity-monitoring \
    --reservation-affinity=any

# Wait until the instance is shut down.
while [ "$(gcloud compute instances describe ${INSTANCE} --project ${PROJECT} --zone ${ZONE} --format='value(status)')" == "RUNNING" ];
do
    sleep 60
done

# Delete the instance.
gcloud -q compute instances delete ${INSTANCE} --project ${PROJECT} --zone=${ZONE}