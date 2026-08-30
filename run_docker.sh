#!/bin/bash
set -e

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE_NAME="geant4-stopping:latest"

# Ensure docker is available
if ! command -v docker &> /dev/null; then
    echo "Error: 'docker' command not found. Please install and start Docker."
    exit 1
fi

# Function to build Docker image
build_image() {
    echo "==> Building Docker image '${IMAGE_NAME}'..."
    docker build -t "${IMAGE_NAME}" -f "${SCRIPT_DIR}/docker/Dockerfile" "${SCRIPT_DIR}"
    echo "==> Docker image built successfully."
}

# Check if image exists, build if missing
if [[ "$(docker images -q ${IMAGE_NAME} 2> /dev/null)" == "" ]]; then
    echo "==> Docker image '${IMAGE_NAME}' not found locally."
    build_image
fi

# Handle build-image command explicitly
if [ "$1" == "build-image" ]; then
    build_image
    exit 0
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
