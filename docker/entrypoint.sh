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

# Helper: resolve the reaction case directory from a macro path like
# reactions/<case>/macros/<name>.mac (empty when not a case macro).
resolve_case_dir() {
    local MACRO="$1"
    MACRO="${MACRO#./}"
    MACRO="${MACRO#/workspace/}"
    case "$MACRO" in
        reactions/*/macros/*) echo "${MACRO%%/macros/*}" ;;
        *) echo "" ;;
    esac
}

# Helper: run the simulation. For case macros, cd into
# reactions/<case>/outputs/ so the reaction inputs are found via
# ../angular_distribution.txt and ../cross_sections.txt, and the CSV
# and plots land in that case's outputs/ instead of a shared build/.
run_simulation() {
    local MACRO="$1"
    local CASE_DIR RELMACRO OUT_DIR
    CASE_DIR="$(resolve_case_dir "$MACRO")"
    if [ ! -f "/workspace/build/stopping" ]; then
        echo "==> Executable not found. Compiling first..."
        build_project
    fi
    if [ -n "$CASE_DIR" ]; then
        OUT_DIR="/workspace/$CASE_DIR/outputs"
        RELMACRO="../macros/$(basename "$MACRO")"
        mkdir -p "$OUT_DIR"
        cd "$OUT_DIR"
        echo "==> Running simulation for case '$CASE_DIR' with macro: $MACRO"
        /workspace/build/stopping "$RELMACRO"
    else
        cd /workspace/build
        echo "==> Running simulation with macro: $MACRO"
        ./stopping "$MACRO"
    fi
}

# Helper: generate plots. For case runs the CSV lives in
# reactions/<case>/outputs/; abort with a clear message when the
# simulation has not produced it.
generate_plots() {
    local CASE_DIR="$1"
    local OUT_DIR LABEL
    if [ -n "$CASE_DIR" ]; then
        OUT_DIR="/workspace/$CASE_DIR/outputs"
        LABEL="$(basename "$CASE_DIR")"
    else
        OUT_DIR="/workspace/build"
    fi
    cd "$OUT_DIR"
    echo "==> Generating analysis histograms and plots..."
    if [ ! -f "stopping_data_nt_RecoilData.csv" ]; then
        echo "Error: stopping_data_nt_RecoilData.csv not found in $OUT_DIR."
        echo "The simulation must run first, e.g.:"
        echo "  ./run_docker.sh all <path-to-macro>"
        exit 1
    fi
    python3 /workspace/plot_histogram.py ${LABEL:+"$LABEL"}
}

case "$COMMAND" in
    compile|build)
        shift
        echo "==> Building Geant4 Stopping Simulation..."
        build_project "$@"
        echo "==> Build completed successfully: /workspace/build/stopping"
        ;;
        
    run)
        shift
        MACRO="${1:-$(ls /workspace/reactions/*/macros/*.mac 2>/dev/null | head -1)}"
        run_simulation "$MACRO"
        ;;
        
    plot)
        shift
        CASE_DIR="$(resolve_case_dir "$1")"
        if [ -n "$CASE_DIR" ]; then
            generate_plots "$CASE_DIR"
        elif [ -n "$1" ] && [ -d "/workspace/reactions/$1" ]; then
            generate_plots "reactions/$1"
        else
            generate_plots ""
        fi
        ;;
        
    all)
        shift
        MACRO="${1:-$(ls /workspace/reactions/*/macros/*.mac 2>/dev/null | head -1)}"
        echo "==> Running full workflow (compile -> simulate -> plot)..."
        build_project
        run_simulation "$MACRO"
        generate_plots "$(resolve_case_dir "$MACRO")"
        echo "==> All steps completed! Check the case outputs directory for generated PDF/PNG plots and CSV data."
        ;;
        
    *)
        exec "$@"
        ;;
esac
