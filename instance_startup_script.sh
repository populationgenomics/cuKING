#!/usr/bin/env bash

set -ex

# Enable logging and monitoring metrics.
export DEBIAN_FRONTEND=noninteractive
curl -sSO https://dl.google.com/cloudagents/add-google-cloud-ops-agent-repo.sh
sudo bash add-google-cloud-ops-agent-repo.sh --also-install

# Install Docker and NVIDIA drivers for CUDA 11.7.
apt-get update
apt-get install -y docker.io gcc make pkg-config
curl -O https://us.download.nvidia.com/XFree86/Linux-x86_64/515.57/NVIDIA-Linux-x86_64-515.57.run
sh NVIDIA-Linux-x86_64-515.57.run -s

# Install the NVIDIA Container Toolkit.
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
curl -s -L https://nvidia.github.io/libnvidia-container/ubuntu22.04/libnvidia-container.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
apt-get update
apt-get install -y nvidia-docker2
systemctl restart docker

# Register the Artifact Registry Docker repository.
gcloud -q auth configure-docker us-central1-docker.pkg.dev

# Run cuKING.
docker run --name cuking --gpus all \
    us-central1-docker.pkg.dev/${PROJECT}/images/cuking:latest cuking \
    --input_uri=${INPUT_URI} \
    --output_uri=${OUTPUT_URI} \
    --requester_pays_project=${PROJECT} \
    --king_coeff_threshold=${KING_COEFF_THRESHOLD}

# Copy the logs.
docker logs cuking | gsutil -u ${PROJECT} cp - ${LOGGING_OUTPUT}

# Shut down the instance.
shutdown -h now
