# Multi-architecture Dockerfile for Geant4 Stopping Simulation
# Supports both linux/arm64 (Apple Silicon M1/M2/M3) and linux/amd64 (Intel/AMD)
#
# NOTE: versions are pinned on purpose. The base image is pinned by digest and
# geant4 by exact version, so the prebuilt image published to GHCR is
# reproducible. To upgrade, bump these deliberately and re-verify a full
# `./run_docker.sh` pipeline before re-tagging a release.

# condaforge/miniforge3:latest as of 2026-10 (linux/amd64 manifest)
FROM condaforge/miniforge3@sha256:f752860f77bd417aa4db35be3d2e906fdb740fb80f40755ed651fff9e8873aad

LABEL maintainer="Razvan Lica <razvan.lica@cern.ch>"
LABEL description="Geant4 Stopping Simulation environment with Geant4, CMake, C++ compilers, and Python plotting dependencies"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Install system utilities and build tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Geant4, compilation toolchain, and Python analysis packages via conda-forge
RUN mamba install -y -c conda-forge \
    geant4=11.4.3 \
    cmake \
    make \
    cxx-compiler \
    c-compiler \
    python>=3.10 \
    numpy \
    matplotlib \
    pandas \
    && mamba clean -afy

# Install expat/zlib development files in a SEPARATE transaction.
# Geant4 CMake config does find_dependency(EXPAT) and find_dependency(ZLIB),
# which needs headers (zlib.h) and dev symlinks (libz.so); conda-forge splits
# these into runtime-only packages (libexpat, libzlib) that geant4 pulls in.
# IMPORTANT: when requested in the same transaction as the packages above, the
# solver silently substitutes the zlib request with zlib-ng (which provides the
# zlib virtual but ships no headers), so the dev files never land in the image
# and find_package(Geant4) fails. Installing them in their own transaction
# reliably links the real zlib/expat packages.
RUN mamba install -y -c conda-forge \
    expat \
    zlib \
    && mamba clean -afy

# Make CMake search the conda prefix for find_package (Geant4, EXPAT, ...)
ENV CMAKE_PREFIX_PATH=/opt/conda

# Set working directory
WORKDIR /workspace

# Copy entrypoint script
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]
