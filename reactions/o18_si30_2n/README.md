# Reaction case: 18O + 30Si -> 16O + 32Si* (2n transfer), E* = 5.5-6.5 MeV

Beam: 18O at 20-26 MeV on a 30Si target. Fresco calculations for this
case live in `Frescox/si30-o18-2n/` (per-energy subfolders E22..E26).

## Contents

- `angular_distribution.txt` - CM dSigma/dOmega extracted from the fresco
  output (`angular_distribution_fresco.txt`), consumed by the Geant4
  simulation. Regenerate with:
  `python3 extract_fresco.py angular_distribution_fresco.txt angular_distribution.txt`
- `cross_sections.txt` - total 2n cross section (mb) vs beam energy (MeV),
  used for depth-weighted reaction probability in thick targets.
- `macros/` - one macro per target backing material:
  - `ta_backing.mac` - tantalum backing
  - `si30_backing.mac` - 30Si powder
  - `ti46_backing.mac` - 46Ti backing
- `outputs/` - CSV and plots produced by the simulation (gitignored).

## Run

From the repository root:

    ./run_docker.sh all reactions/o18_si30_2n/macros/si30_backing.mac

### Run natively (without Docker)

Requires a local Geant4 installation and the compiled `build/stopping`
executable (see the main README, "Option B: Native Build"). Source your
Geant4 environment first if needed (e.g.
`source <geant4-install>/bin/geant4.sh`).

The simulation must run **inside this case's `outputs/` directory**, so
the reaction inputs are found one level up (`../angular_distribution.txt`,
`../cross_sections.txt`) and the CSV/plots are written next to them:

    cd reactions/o18_si30_2n
    mkdir -p outputs
    cd outputs
    ../../../build/stopping ../macros/si30_backing.mac

Optionally regenerate the histograms/plots from the same directory:

    python3 ../../../plot_histogram.py o18_si30_2n

Or use the helper script from the repository root, which builds, resolves
all paths, and plots automatically:

    ./run_geant4_stopping.sh reactions/o18_si30_2n/macros/si30_backing.mac
