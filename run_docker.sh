#!/bin/bash
set -e

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCAL_IMAGE="geant4-stopping:latest"
REGISTRY_IMAGE="ghcr.io/rlica/geant4-stopping:latest"
IMAGE_NAME="$LOCAL_IMAGE"

# Ensure docker is available
if ! command -v docker &> /dev/null; then
    echo "Error: 'docker' command not found. Please install and start Docker."
    echo "Tip: run ./install_docker.sh to set Docker up automatically."
    exit 1
fi

# Function to build Docker image
build_image() {
    echo "==> Building Docker image '${LOCAL_IMAGE}'..."
    docker build -t "${LOCAL_IMAGE}" -f "${SCRIPT_DIR}/docker/Dockerfile" "${SCRIPT_DIR}"
    echo "==> Docker image built successfully."
}

# Handle build-image command explicitly (maintainer path, skips pull)
if [ "$1" == "build-image" ]; then
    build_image
    exit 0
fi

# Resolve which image to use, in order of preference:
#   1. a local image already present        -> use it (fast, works offline)
#   2. the prebuilt image published on GHCR -> pull it (normal user path)
#   3. fallback                             -> build locally from ./docker
if [[ "$(docker images -q ${LOCAL_IMAGE} 2> /dev/null)" != "" ]]; then
    echo "==> Using existing local image '${LOCAL_IMAGE}'."
elif docker pull "${REGISTRY_IMAGE}"; then
    echo "==> Using prebuilt image '${REGISTRY_IMAGE}'."
    IMAGE_NAME="${REGISTRY_IMAGE}"
else
    echo "==> Could not pull prebuilt image '${REGISTRY_IMAGE}'."
    echo "    Building locally from ${SCRIPT_DIR}/docker/Dockerfile instead..."
    build_image
fi

# If no arguments provided, show usage and default to full run
if [ $# -eq 0 ]; then
    echo "Usage: ./run_docker.sh [command] [options]"
    echo ""
    echo "Commands:"
    echo "  all [macro]       Build, run simulation, and plot (default macro: ../run.mac)"
    echo "  compile           Compile the C++ simulation"
    echo "  run [macro]       Run the simulation with specified macro"
    echo "  plot              Generate Python histograms/plots"
    echo "  bash              Open interactive shell inside container"
    echo "  build-image       Rebuild the Docker image"
    echo "  <cmd>             Run any custom command inside the container"
    echo ""
    echo "==> Executing default workflow: all ../run.mac"
    ACTION="all"
else
    ACTION="$1"
    shift
fi

# Set interactive flags if terminal is attached
DOCKER_FLAGS=""
if [ -t 0 ] && [ -t 1 ]; then
    DOCKER_FLAGS="-it"
elif [ -t 0 ]; then
    DOCKER_FLAGS="-i"
fi

# Run docker container mounting current directory
docker run --rm ${DOCKER_FLAGS} \
    -v "${SCRIPT_DIR}:/workspace" \
    -w /workspace \
    "${IMAGE_NAME}" "${ACTION}" "$@"
