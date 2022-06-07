FROM nvidia/cuda:11.7.0-devel-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive

# Remove nvidia repositories to work around https://github.com/NVIDIA/nvidia-docker/issues/1402
RUN rm /etc/apt/sources.list.d/cuda-ubuntu2204-x86_64.list && \
    apt update && apt install --no-install-recommends -y \
        apt-transport-https \
        apt-utils \
        automake \
        build-essential \
        ca-certificates \
        ccache \
        cmake \
        cmake \
        curl \
        g++ \
        gcc \
        git \
        libc-ares-dev \
        libc-ares2 \
        libhts-dev \
        libre2-dev \
        libssl-dev \
        m4 \
        make \
        pkg-config \
        tar \
        wget \
        zlib1g-dev \
        libzstd-dev

RUN mkdir -p /deps/abseil-cpp && cd /deps/abseil-cpp && \
    curl -sSL https://github.com/abseil/abseil-cpp/archive/20211102.0.tar.gz | tar -xzf - --strip-components=1 && \
    cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_STANDARD=17 \
      -DBUILD_TESTING=OFF \
      -DBUILD_SHARED_LIBS=yes \
      -S . -B cmake-out && \
    cmake --build cmake-out --target install -- -j 16 && \
    ldconfig

RUN mkdir -p /deps/protobuf && cd /deps/protobuf && \
    curl -sSL https://github.com/protocolbuffers/protobuf/archive/v21.1.tar.gz | tar -xzf - --strip-components=1 && \
    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DBUILD_SHARED_LIBS=yes \
        -Dprotobuf_BUILD_TESTS=OFF \
        -Dprotobuf_ABSL_PROVIDER=package \
        -S . -B cmake-out && \
    cmake --build cmake-out --target install -- -j 16 && \
    ldconfig

RUN mkdir -p /deps/grpc && cd /deps/grpc && \
    curl -sSL https://github.com/grpc/grpc/archive/v1.46.3.tar.gz | tar -xzf - --strip-components=1 && \
    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DBUILD_SHARED_LIBS=yes \
        -DgRPC_INSTALL=ON \
        -DgRPC_BUILD_TESTS=OFF \
        -DgRPC_ABSL_PROVIDER=package \
        -DgRPC_CARES_PROVIDER=package \
        -DgRPC_PROTOBUF_PROVIDER=package \
        -DgRPC_RE2_PROVIDER=package \
        -DgRPC_SSL_PROVIDER=package \
        -DgRPC_ZLIB_PROVIDER=package \
        -S . -B cmake-out && \
    cmake --build cmake-out --target install -- -j 16 && \
    ldconfig

RUN mkdir -p /deps/crc32c && cd /deps/crc32c && \
    curl -sSL https://github.com/google/crc32c/archive/1.1.2.tar.gz | tar -xzf - --strip-components=1 && \
    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DBUILD_SHARED_LIBS=yes \
        -DCRC32C_BUILD_TESTS=OFF \
        -DCRC32C_BUILD_BENCHMARKS=OFF \
        -DCRC32C_USE_GLOG=OFF \
        -S . -B cmake-out && \
    cmake --build cmake-out --target install -- -j 16 && \
    ldconfig

RUN mkdir -p /deps/json && cd /deps/json && \
    curl -sSL https://github.com/nlohmann/json/archive/v3.10.5.tar.gz | tar -xzf - --strip-components=1 && \
    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DBUILD_SHARED_LIBS=yes \
        -DBUILD_TESTING=OFF \
        -DJSON_BuildTests=OFF \
        -S . -B cmake-out && \
    cmake --build cmake-out --target install -- -j 16 && \
    ldconfig

RUN mkdir -p /deps/google-cloud-cpp && cd /deps/google-cloud-cpp && \
    curl -sSL https://github.com/googleapis/google-cloud-cpp/archive/refs/tags/v1.41.0.tar.gz | tar -xzf - --strip-components=1 && \
    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DCMAKE_INSTALL_MESSAGE=NEVER \
        -DBUILD_TESTING=OFF \
        -DGOOGLE_CLOUD_CPP_ENABLE_EXAMPLES=OFF \
        -H. -B cmake-out && \
    cmake --build cmake-out --target install -- -j 16 && \
    ldconfig

# Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt update && apt install -y google-cloud-sdk

COPY . /app/
WORKDIR /app

RUN rm -rf build && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    cmake --build . -j 16 && \
    mv cuking gvcf2cuking /usr/local/bin