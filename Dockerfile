FROM nvidia/cuda:11.7.0-devel-ubuntu22.04 AS dev

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
        apt-transport-https \
        ca-certificates \
        cmake \
        curl \
        g++ \
        git \
        gnupg \
        libc-ares-dev \
        libcurl4-gnutls-dev \
        libre2-dev \
        libssl-dev \
        make \
        ninja-build \
        pkg-config \
        zlib1g-dev

# Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update && apt-get install -y google-cloud-sdk

RUN mkdir -p /deps/abseil-cpp && cd /deps/abseil-cpp && \
    curl -sSL https://github.com/abseil/abseil-cpp/archive/refs/tags/20220623.0.tar.gz | tar -xzf - --strip-components=1 && \
    cmake \
        -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DABSL_PROPAGATE_CXX_STD=ON \
        -DBUILD_TESTING=OFF \
        -DBUILD_SHARED_LIBS=yes \
        -S . -B cmake-out && \
    cmake --build cmake-out --target install && \
    ldconfig

RUN mkdir -p /deps/protobuf && cd /deps/protobuf && \
    curl -sSL https://github.com/protocolbuffers/protobuf/archive/v21.1.tar.gz | tar -xzf - --strip-components=1 && \
    cmake \
        -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DBUILD_SHARED_LIBS=yes \
        -Dprotobuf_BUILD_TESTS=OFF \
        -Dprotobuf_ABSL_PROVIDER=package \
        -S . -B cmake-out && \
    cmake --build cmake-out --target install && \
    ldconfig

RUN mkdir -p /deps/grpc && cd /deps/grpc && \
    curl -sSL https://github.com/grpc/grpc/archive/v1.46.3.tar.gz | tar -xzf - --strip-components=1 && \
    cmake \
        -G Ninja \
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
    cmake --build cmake-out --target install && \
    ldconfig

RUN mkdir -p /deps/crc32c && cd /deps/crc32c && \
    curl -sSL https://github.com/google/crc32c/archive/1.1.2.tar.gz | tar -xzf - --strip-components=1 && \
    cmake \
        -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DBUILD_SHARED_LIBS=yes \
        -DCRC32C_BUILD_TESTS=OFF \
        -DCRC32C_BUILD_BENCHMARKS=OFF \
        -DCRC32C_USE_GLOG=OFF \
        -S . -B cmake-out && \
    cmake --build cmake-out --target install && \
    ldconfig

RUN mkdir -p /deps/json && cd /deps/json && \
    curl -sSL https://github.com/nlohmann/json/archive/v3.10.5.tar.gz | tar -xzf - --strip-components=1 && \
    cmake \
        -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DBUILD_SHARED_LIBS=yes \
        -DBUILD_TESTING=OFF \
        -DJSON_BuildTests=OFF \
        -S . -B cmake-out && \
    cmake --build cmake-out --target install && \
    ldconfig

RUN mkdir -p /deps/google-cloud-cpp && cd /deps/google-cloud-cpp && \
    curl -sSL https://github.com/googleapis/google-cloud-cpp/archive/refs/tags/v1.41.0.tar.gz | tar -xzf - --strip-components=1 && \
    cmake \
        -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DBUILD_SHARED_LIBS=yes \
        -DGOOGLE_CLOUD_CPP_ENABLE=storage \
        -DBUILD_TESTING=OFF \
        -DGOOGLE_CLOUD_CPP_ENABLE_EXAMPLES=OFF \
        -DCMAKE_INSTALL_MESSAGE=NEVER \
        -H. -B cmake-out && \
    cmake --build cmake-out --target install && \
    ldconfig

RUN mkdir -p /deps/arrow && cd /deps/arrow && \
    curl -sSL https://github.com/apache/arrow/archive/refs/tags/apache-arrow-8.0.0.tar.gz | tar -xzf - --strip-components=1 && \
    cmake cpp \
        -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DARROW_BUILD_STATIC=OFF \
        -DARROW_PARQUET=ON \
        -DARROW_WITH_ZSTD=ON \
        -B cmake-out && \
    cmake --build cmake-out --target install && \
    ldconfig

# extract-elf-so tars .so files to create small Docker images.
RUN curl -sSL -o /deps/extract-elf-so https://github.com/William-Yeh/extract-elf-so/releases/download/v0.6/extract-elf-so_static_linux-amd64 && \
    chmod +x /deps/extract-elf-so

FROM dev as extract

COPY . /app/
WORKDIR /app

RUN cmake \
        -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -B cmake-out && \
    cmake --build cmake-out

RUN /deps/extract-elf-so --cert /app/cmake-out/cuking

FROM nvidia/cuda:11.7.0-base-ubuntu22.04 AS minimal

RUN --mount=type=bind,from=extract,source=/app/rootfs.tar,target=/rootfs.tar \
    tar xf /rootfs.tar && \
    ldconfig
