FROM ubuntu:26.04 AS build

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install --yes --no-install-recommends \
    build-essential \
    ca-certificates \
    cmake \
    ninja-build \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /source
COPY . .

RUN cmake -S . -B /build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \
        -DCHARGEFW_BUILD_CLI=ON \
        -DCHARGEFW_BUILD_PYTHON=OFF \
        -DCHARGEFW_BUILD_TESTS=OFF \
        -DCHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS=OFF \
    && cmake --build /build --target chargefw_cli --parallel \
    && cmake --install /build --prefix /opt/chargefw --strip

FROM ubuntu:26.04

RUN apt-get update && apt-get install --yes --no-install-recommends \
    libstdc++6 \
    zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/chargefw /opt/chargefw

ENTRYPOINT ["/opt/chargefw/bin/chargefw"]
CMD ["--help"]
