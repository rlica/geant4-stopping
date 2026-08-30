#!/bin/bash

# Source Geant4 environment
source /home/razvan/geant4-v11.4.1-install/bin/geant4.sh

# Create build directory if it doesn't exist
mkdir -p build
cd build

# Compile the simulation
cmake ..
make -j4

# Run the simulation
./stopping ../run.mac

# Generate histograms and plots
python3 plot_histogram.py

# Display PDF plots automatically
xpdf stopping_times_histogram.pdf &
xpdf recoil_energies_histogram.pdf &
xpdf angular_distribution.pdf &

echo "Simulation and plotting completed successfully!"
