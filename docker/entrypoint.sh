#!/bin/bash
set -e

# Activate conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate base

# Helper function to configure and build
build_project() {
    mkdir -p /workspace/build
    cd /workspace/build
    # If CMakeCache.txt exists from a different path/system, clear stale cache
    if [ -f "CMakeCache.txt" ]; then
        if grep -q -v "/workspace" CMakeCache.txt 2>/dev/null; then
            echo "==> Stale CMake cache detected from another environment. Resetting CMake cache..."
            rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake Makefile
        fi
    fi
    cmake /workspace
    make -j$(nproc 2>/dev/null || echo 4) "$@"
}

# If no arguments provided, open bash
if [ $# -eq 0 ]; then
    set -- bash
fi

COMMAND="$1"

case "$COMMAND" in
    compile|build)
        shift
        echo "==> Building Geant4 Stopping Simulation..."
        build_project "$@"
        echo "==> Build completed successfully: /workspace/build/stopping"
        ;;
        
    run)
        shift
        MACRO="${1:-../run.mac}"
        if [ ! -f "/workspace/build/stopping" ]; then
            echo "==> Executable not found. Compiling first..."
            build_project
        fi
        cd /workspace/build
        echo "==> Running simulation with macro: $MACRO"
        ./stopping "$MACRO"
        ;;
        
    plot)
        shift
        cd /workspace/build
        echo "==> Generating analysis histograms and plots..."
        if [ -f "plot_histogram.py" ]; then
            python3 plot_histogram.py "$@"
        elif [ -f "/workspace/plot_histogram.py" ]; then
            python3 /workspace/plot_histogram.py "$@"
        else
            echo "Error: plot_histogram.py not found."
            exit 1
        fi
        ;;
        
    all)
        shift
        MACRO="${1:-../run.mac}"
        echo "==> Running full workflow (compile -> simulate -> plot)..."
        build_project
        cd /workspace/build
        echo "==> Running simulation with macro: $MACRO"
        ./stopping "$MACRO"
        echo "==> Generating plots..."
        if [ -f "plot_histogram.py" ]; then
            python3 plot_histogram.py
        elif [ -f "/workspace/plot_histogram.py" ]; then
            python3 /workspace/plot_histogram.py
        fi
        echo "==> All steps completed! Check build/ for generated PDF/PNG plots and CSV data."
        ;;
        
    *)
        exec "$@"
        ;;
esac
