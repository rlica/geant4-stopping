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

# Convenience: passing a macro path directly (e.g. when invoking the
# container manually) runs the full workflow with it instead of trying
# to execute the macro file as a command.
if [[ "$COMMAND" == *.mac ]]; then
    set -- "$COMMAND" "${@:2}"
    COMMAND="all"
fi

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

# Helper: outputs directory for a case + macro: outputs_<macro-stem> so
# different macros for the same reaction never overwrite each other.
outputs_dir_for() {
    local CASE_DIR="$1" MACRO="$2"
    echo "/workspace/$CASE_DIR/outputs_$(basename "$MACRO" .mac)"
}

# Helper: newest outputs directory of a case (any outputs_<macro> dir,
# falling back to the legacy unsuffixed outputs/). Used by 'plot <case>'.
newest_outputs_dir() {
    local CASE_DIR="$1" NEWEST
    NEWEST="$(ls -dt "/workspace/$CASE_DIR/"outputs_* 2>/dev/null | head -1)"
    if [ -z "$NEWEST" ] && [ -d "/workspace/$CASE_DIR/outputs" ]; then
        NEWEST="/workspace/$CASE_DIR/outputs"
    fi
    echo "$NEWEST"
}

# Helper: run the simulation. For case macros, cd into
# reactions/<case>/outputs_<macro-stem>/ so the reaction inputs are found
# via ../angular_distribution.txt and ../cross_sections.txt, and the CSV
# and plots land in that macro's own outputs directory instead of a
# shared build/.
run_simulation() {
    local MACRO="$1"
    local CASE_DIR RELMACRO OUT_DIR
    CASE_DIR="$(resolve_case_dir "$MACRO")"
    if [ ! -f "/workspace/build/stopping" ]; then
        echo "==> Executable not found. Compiling first..."
        build_project
    fi
    if [ -n "$CASE_DIR" ]; then
        OUT_DIR="$(outputs_dir_for "$CASE_DIR" "$MACRO")"
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

# Helper: generate plots. For case runs the CSV lives in the macro's
# outputs_<macro-stem>/ directory; abort with a clear message when the
# simulation has not produced it.
generate_plots() {
    local CASE_DIR="$1" MACRO="$2"
    local OUT_DIR LABEL
    if [ -n "$CASE_DIR" ]; then
        if [ -n "$MACRO" ]; then
            OUT_DIR="$(outputs_dir_for "$CASE_DIR" "$MACRO")"
        else
            OUT_DIR="$(newest_outputs_dir "$CASE_DIR")"
        fi
        if [ -z "$OUT_DIR" ]; then
            echo "Error: no outputs directory found in /workspace/$CASE_DIR."
            echo "The simulation must run first, e.g.:"
            echo "  ./run_docker.sh all <path-to-macro>"
            exit 1
        fi
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
        MACRO="${1:-}"
        CASE_DIR="$(resolve_case_dir "$MACRO")"
        if [ -n "$CASE_DIR" ]; then
            generate_plots "$CASE_DIR" "$MACRO"
        elif [ -n "$MACRO" ] && [ -d "/workspace/reactions/$MACRO" ]; then
            generate_plots "reactions/$MACRO" ""
        else
            generate_plots "" ""
        fi
        ;;

    all)
        shift
        MACRO="${1:-$(ls /workspace/reactions/*/macros/*.mac 2>/dev/null | head -1)}"
        echo "==> Running full workflow (compile -> simulate -> plot)..."
        build_project
        run_simulation "$MACRO"
        generate_plots "$(resolve_case_dir "$MACRO")" "$MACRO"
        echo "==> All steps completed! Check the macro's outputs directory for generated PDF/PNG plots and CSV data."
        ;;
        
    *)
        exec "$@"
        ;;
esac
