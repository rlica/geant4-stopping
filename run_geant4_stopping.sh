#!/bin/bash
set -euo pipefail

INVOCATION_DIR="$PWD"

# Source Geant4 environment. Override the install location with
# GEANT4_INSTALL=/path/to/geant4-install if yours lives elsewhere.
G4ENV="${GEANT4_INSTALL:-$HOME/geant4/geant4-v11.4.3-install}"
if [ -f "$G4ENV/bin/geant4.sh" ]; then
    source "$G4ENV/bin/geant4.sh"
else
    echo "ERROR: Geant4 environment not found at: $G4ENV/bin/geant4.sh" >&2
    echo "Install it with ./install_geant4.sh, or point GEANT4_INSTALL at your existing install." >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Resolve the macro argument BEFORE changing directory anywhere, so both
# relative paths (from the invocation directory) and absolute paths work.
MACRO="${1:-}"
if [ -n "$MACRO" ]; then
    if [ "${MACRO:0:1}" != "/" ]; then
        MACRO="$INVOCATION_DIR/$MACRO"
    fi
    if [ ! -f "$MACRO" ]; then
        echo "ERROR: macro file not found: $MACRO" >&2
        exit 1
    fi
else
    MACRO="$(ls "$SCRIPT_DIR"/reactions/*/macros/*.mac 2>/dev/null | head -1 || true)"
    if [ -z "$MACRO" ]; then
        echo "ERROR: no macro given and none found under reactions/*/macros/." >&2
        exit 1
    fi
fi
MACRO="$(cd "$(dirname "$MACRO")" && pwd)/$(basename "$MACRO")"

# Compile the simulation
mkdir -p "$SCRIPT_DIR/build"
cd "$SCRIPT_DIR/build"

# If CMakeCache.txt exists from a different environment (e.g. the Docker
# image, where the build lived at /workspace/build) or after the source
# tree moved, clear the stale cache so native cmake reconfigures cleanly.
if [ -f CMakeCache.txt ] && ! grep -q "CMAKE_HOME_DIRECTORY:INTERNAL=$SCRIPT_DIR" CMakeCache.txt; then
    echo "==> Stale CMake cache detected from another environment. Resetting CMake cache..."
    rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake Makefile
fi

# Embed the Geant4 library directory as rpath so build/stopping also runs
# in shells where the Geant4 environment was not sourced. Use legacy
# DT_RPATH (--disable-new-dtags): modern RUNPATH is not consulted for
# transitive dependencies (e.g. libG4zlib.so pulled in by other G4 libs).
cmake -DCMAKE_BUILD_RPATH="$G4ENV/lib" -DCMAKE_INSTALL_RPATH="$G4ENV/lib" \
    -DCMAKE_EXE_LINKER_FLAGS="-Wl,--disable-new-dtags" "$SCRIPT_DIR"
make -j"$(nproc 2>/dev/null || echo 4)"

# Run the simulation inside the case's outputs/ directory so inputs
# (../angular_distribution.txt, ../cross_sections.txt) and outputs (CSV,
# plots) stay per reaction case.
CASE_DIR="$(dirname "$(dirname "$MACRO")")"        # reactions/<case>
mkdir -p "$CASE_DIR/outputs"
cd "$CASE_DIR/outputs"
"$SCRIPT_DIR/build/stopping" "../macros/$(basename "$MACRO")"

# Generate histograms and plots
python3 "$SCRIPT_DIR/plot_histogram.py"

# Display PDF plots automatically (optional viewer)
if command -v xpdf >/dev/null 2>&1; then
    for pdf in stopping_times_histogram.pdf recoil_energies_histogram.pdf angular_distribution.pdf; do
        [ -f "$pdf" ] && xpdf "$pdf" &
    done
fi

echo "Simulation and plotting completed successfully!"

