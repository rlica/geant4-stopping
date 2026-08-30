# Multi-architecture Dockerfile for Geant4 Stopping Simulation
# Supports both linux/arm64 (Apple Silicon M1/M2/M3) and linux/amd64 (Intel/AMD)

FROM condaforge/miniforge3:latest

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
    geant4 \
    cmake \
    make \
    cxx-compiler \
    c-compiler \
    python>=3.10 \
    numpy \
    matplotlib \
    pandas \
    && mamba clean -afy

# Set working directory
WORKDIR /workspace

# Copy entrypoint script
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]
