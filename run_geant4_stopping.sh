#!/bin/bash

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

# Create build directory if it doesn't exist
mkdir -p "$SCRIPT_DIR/build"
cd "$SCRIPT_DIR/build"

# Compile the simulation
cmake "$SCRIPT_DIR"
make -j4

# Run the simulation: pass a reaction macro as $1, otherwise the first
# macro under reactions/*/macros/ is used. The simulation runs inside the
# case's outputs/ directory so inputs (../angular_distribution.txt,
# ../cross_sections.txt) and outputs (CSV, plots) stay per reaction case.
MACRO="${1:-$(ls "$SCRIPT_DIR"/reactions/*/macros/*.mac 2>/dev/null | head -1 || true)}"
if [ -n "$MACRO" ]; then
    MACRO="$(cd "$(dirname "$MACRO")" && pwd)/$(basename "$MACRO")"
    CASE_DIR="$(dirname "$(dirname "$MACRO")")"        # reactions/<case>
    mkdir -p "$CASE_DIR/outputs"
    cd "$CASE_DIR/outputs"
    "$SCRIPT_DIR/build/stopping" "../macros/$(basename "$MACRO")"
else
    echo "No reaction macros found under reactions/*/macros/."
    "$SCRIPT_DIR/build/stopping"
    cd "$SCRIPT_DIR/build"
fi

# Generate histograms and plots
python3 "$SCRIPT_DIR/plot_histogram.py"

# Display PDF plots automatically
xpdf stopping_times_histogram.pdf &
xpdf recoil_energies_histogram.pdf &
xpdf angular_distribution.pdf &

echo "Simulation and plotting completed successfully!"
