# Geant4 Stopping Simulation

## Introduction
The `geant4-stopping` project is a Geant4-based application designed to simulate nuclear transfer reactions (e.g., 18O + 30Si -> 16O + 32Si*) and the subsequent slowing-down process of heavy ions in a solid target medium. It employs a custom minimal physics list tuned for high-accuracy, low-energy electromagnetic interactions (including atomic deexcitation, fluorescence, and Auger electron emission) alongside step limiters to ensure precise $dE/dx$ integration. 

The primary objective is to simulate the primary beam's energy loss in a "thick" target, dynamically model the interaction point using energy-dependent cross sections, and track the kinematics, stopping time, and stopping depth of the recoiling ions. It incorporates angular distributions from FRESCOX calculations to sample realistic ejectile scattering angles in the Center-of-Mass (CM) frame, adjusting reaction Q-values for specific excitation energies. The energy-dependent transfer cross sections used in thick-target mode are produced with the GRAZING code (see `GRAZING_9/` and the *Reaction Cases* section below).

## Quickstart

### Option A: Using Docker (Recommended for any OS / Apple Silicon / Linux / Windows)

No local installation of Geant4 or compilation tools is required, and **no image build is needed**: on first use, `run_docker.sh` automatically pulls the prebuilt, multi-architecture image (`ghcr.io/rlica/geant4-stopping:latest`) which contains Geant4, the compilers, all datasets, and the Python analysis packages. The image is cached locally, so subsequent runs start instantly.

1. **One-time prerequisite: install Docker** (skip if `docker --version` already works):
   ```bash
   ./install_docker.sh
   ```
   Then open a **new terminal** (or run `newgrp docker` on Linux) so the group change takes effect.

2. **Run full simulation & analysis pipeline with one command:**
   ```bash
   ./run_docker.sh
   # or explicitly with a specific reaction macro:
   ./run_docker.sh all reactions/o18_si30_2n/macros/si30_backing.mac
   ```

3. **Other convenient Docker commands:**
   ```bash
   ./run_docker.sh compile          # Compile C++ simulation inside Docker
   ./run_docker.sh run reactions/s36_o18_2n/macros/ta_backing.mac # Run a specific macro
   ./run_docker.sh plot reactions/o18_si30_2n # Regenerate plots for a case
   ./run_docker.sh bash             # Open an interactive shell inside Docker
   ./run_docker.sh build-image      # Maintainers: rebuild the image locally
   ```

   If the prebuilt image cannot be pulled (e.g. offline), the script falls back to building the image locally from `docker/Dockerfile`.

#### Troubleshooting Docker runs

- **`./run_docker.sh compile` fails with `Could NOT find ZLIB` / `Could NOT find EXPAT`** (or any CMake configure error) right after cloning or updating the repo: your cached local `geant4-stopping:latest` image was built from an older Dockerfile and is stale. Rebuild it once and retry:
   ```bash
   ./run_docker.sh build-image
   ```
- **`docker` command not found:** run `./install_docker.sh` first. On Linux it adds your user to the `docker` group, which only takes effect in a new terminal (or after `newgrp docker`).
- **Apple Silicon / ARM Macs and Intel/AMD PCs are both supported:** the prebuilt GHCR image is multi-architecture, and `./run_docker.sh build-image` natively builds for the local architecture.

---

### Option B: Native Build (Requires local Geant4)

#### Prerequisites
- **Geant4** (version 11.x or similar) with appropriate physics datasets.
- **CMake** (version 3.16 or higher).
- **Python 3** with `matplotlib` and `pandas` (for plotting results).

#### Building the Project
To compile the simulation natively, you can use the provided bash script or build it manually using CMake.

**Using the provided script:**
```bash
./run_geant4_stopping.sh
```
This script sets up the environment, builds the executable, runs the default macro, and generates PDF plots automatically.

**Manual Build:**
```bash
mkdir -p build && cd build
cmake ..
make -j4
```

### Running the Simulation
Execute the simulation by providing a macro file. For reaction cases, run
from the case's `outputs/` directory so the inputs are found one level up
and the CSV/plots land next to it:
```bash
mkdir -p reactions/o18_si30_2n/outputs
cd reactions/o18_si30_2n/outputs
../../build/stopping ../macros/si30_backing.mac
```
If run without arguments, the simulation will execute 10 test events:
```bash
./stopping
```

## User Manual

### Output Data
The application leverages the `G4AnalysisManager` to export event-by-event data to a CSV file (`stopping_data.csv`).
Data recorded per event includes:
- `StoppingTime_ps`: The total time taken for the recoil to stop (in picoseconds).
- `RecoilEnergy_MeV`: The initial kinetic energy of the recoil (in MeV) right after the transfer reaction.

Additional scripts (like `plot_histogram.py`) can parse this CSV to produce histograms. Additionally, at the end of every run, the simulation prints a detailed summary to the standard output containing the mean stopping times, interaction depths, standard deviations, and standard errors.

### Macro Commands Reference
The physical target and reaction kinematics are highly configurable via custom UI commands, typically defined in your reaction's macro file (e.g. `reactions/o18_si30_2n/macros/si30_backing.mac`).

#### 1. Target Configuration
You can define a custom, homogeneous stopping medium without recompiling.
- `/detector/customZ <Z>`: Set the Atomic Number (e.g., `14` for Si).
- `/detector/customA <A>`: Set the Mass Number (e.g., `30`).
- `/detector/customDensity <density> <unit>`: Set the mass density (e.g., `1.25 g/cm3`).
- `/detector/customSymbol <Symbol>`: Set the chemical symbol (e.g., `Si`).

*Note: The target is built dynamically using these properties. Changes must be executed before `/run/initialize`.*

- `/detector/targetThickness <value> <unit>`: Set the full target thickness along the beam axis (e.g. `50 um`, `0.05 mm`). Default: 1 mm.
- `/detector/targetArealDensity <value> mg/cm2`: Set the target thickness as an areal (mass) thickness; it is converted to a physical thickness using `/detector/customDensity` and overrides `targetThickness`.

*Thin-target transmission mode:* if the thickness is small enough that recoils exit the target instead of stopping in it, the run summary additionally reports the number of transmitted recoils and their average energy and velocity (`v/c`) **before and after** the medium (forward exits only; backward exits are counted separately). The CSV gains the columns `RecoilInitBeta`, `ExitEnergy_MeV` and `ExitBeta`, where the exit energy/beta are **signed** (positive = forward exit, negative = backward exit), and `plot_histogram.py` produces an `exit_energy_histogram` showing both populations. For transmitted rows, `StoppingTime_ps` holds the time at which the recoil exits the target and `RecoilEnergy_MeV` its initial energy. A ready-made example is `thin_target.mac` (10 um layer). With the default (thick) target the output is unchanged.*

#### 2. Reaction Kinematics
- `/reaction/thinTarget <true/false>`: 
  - `false` (default): Simulates beam slowing down inside the target. The interaction depth is sampled using the beam's energy loss and external energy-dependent cross-sections.
  - `true`: Bypasses beam slowing down, forcing the reaction to occur directly at the surface of the target.
  - *Fallback:* if `false` but `cross_sections.txt` is missing or empty, the run automatically falls back to `true` (with a warning), since depth-weighted sampling is impossible without cross sections.
- `/reaction/targetZ <Z>` and `/reaction/targetA <A>`: Target nuclide involved in the reaction.
- `/reaction/recoilZ <Z>` and `/reaction/recoilA <A>`: Recoiling heavy fragment.
- `/reaction/ejectileZ <Z>` and `/reaction/ejectileA <A>`: Light ejectile fragment.
- `/reaction/excitationEnergy <energy> <unit>`: Excitation energy E* (MeV) of the recoil fragment used in the Q-value (default `6 MeV`, tuned for the 18O+30Si 2n case). Set this to match the E* window of your FRESCO calculation. *Note for inverse kinematics: the heavy, beam-like fragment must be given as the "recoil" (it is the one whose stopping is simulated) and the light fragment as the "ejectile".*
- `/reaction/thetaMin <angle> <unit>` and `/reaction/thetaMax <angle> <unit>`: Limits the Center-of-Mass scattering angle (e.g., `0 deg` to `180 deg`).

#### 3. Primary Beam Generation
The simulation uses the standard `G4ParticleGun` to define the incoming beam.
- `/gun/particle ion`: Set particle type to ion.
- `/gun/ion <Z> <A>`: Define the beam ion (e.g., `/gun/ion 8 18` for 18O).
- `/gun/energy <energy> <unit>`: Set the initial beam kinetic energy (e.g., `/gun/energy 26 MeV`).

## Reaction Cases

Each simulated reaction lives in its own folder under `reactions/<case>/`,
containing everything that reaction needs:

```
reactions/
├── o18_si30_2n/          # 18O + 30Si -> 16O + 32Si* (2n transfer)
│   ├── README.md         # description, provenance, how to regenerate inputs
│   ├── angular_distribution.txt       # CM dSigma/dOmega (from FRESCO)
│   ├── angular_distribution_fresco.txt# raw FRESCO output
│   ├── cross_sections.txt             # total xsec (mb) vs beam energy (MeV)
│   ├── macros/           # one macro per target backing (Ta, 30Si, 46Ti, ...)
│   └── outputs/          # CSV + plots per backing (gitignored, regenerated)
└── s36_o18_2n/           # 36S + 18O -> 38S + 16O (inverse kinematics)
    └── ...               # same structure
```

- To **switch reaction/backing**, just pass a different macro:
  `./run_docker.sh all reactions/s36_o18_2n/macros/ta_backing.mac`.
  The simulation runs inside `reactions/<case>/outputs/`, so each case's
  CSV (`stopping_data_nt_RecoilData.csv`) and plots never mix.
- To **add a new reaction**: create `reactions/<new-case>/` with
  `angular_distribution.txt`, `cross_sections.txt`, and one macro per
  target backing (see the existing cases and their READMEs). The physics
  inputs are produced with the two codes below.

### Physics input codes: FRESCOX and GRAZING

The reaction inputs for the simulation are produced by two nuclear-reaction
codes kept in this repository, one folder each:

- **`Frescox/`** — coupled-channels code (upstream: LLNL/Frescox, v7.3).
  Produces the **angular distributions** (`angular_distribution.txt`,
  dσ/dΩ vs θ_CM) that determine the ejectile emission angles.
  - Calculations live in one subfolder per reaction (e.g.
    `Frescox/si30-o18-2n/`, with one subfolder per beam energy:
    `E22/` … `E26/`).
  - Build the serial `frescox` binary with the top-level convenience
    makefile: `cd Frescox && make install` (installs to `~/binw/<arch>`;
    see `Frescox/README.md` for details and options).
  - Extract the distribution with
    `python3 extract_fresco.py <raw_fresco_output> angular_distribution.txt`.

- **`GRAZING_9/`** — semi-classical nucleon-transfer code (GRAZING v9).
  Produces the **energy-dependent total transfer cross sections**
  (`cross_sections.txt`, mb vs beam energy in MeV) used for the
  interaction-depth weighting in thick-target mode.
  - `grazing_src_jun2005/` — Fortran sources, with `README`/`HOWTO`
    describing the calculation LEGS (0–4) and the output FOR units.
  - `grazing_exe/` — prebuilt executables (`grazing_9`,
    `grazing_9_64bit`) and the required data tables (`mtablex_2004.dat`,
    `massexp_2004.dat`, `be23.dat`). The code is interactive; the input
    of each run is stored in the `fort.10` file it writes.
  - `32Si/` — GRAZING outputs for the `o18_si30_2n` case, one folder per
    beam energy (`E_20.0` … `E_25.6`). Each contains the LEG 3
    excitation-energy distributions (`fort.95`) that were integrated
    over the E* window 5.5–6.5 MeV to build the
    `reactions/o18_si30_2n/cross_sections.txt` file (3 columns:
    E_lab, total 2n cross section, cross section in the E* window).
  - Note: GRAZING cross sections are only needed for **thick-target**
    runs (`/reaction/thinTarget false`). In thin-target mode the reaction
    occurs at the target surface and the simulation falls back gracefully
    to thin-target behaviour when `cross_sections.txt` is absent
    (see the fallback note above).

**Work procedure for a new reaction:**
1. **FRESCOX**: create `Frescox/<reaction>/`, run the frescox
   calculations per beam energy, and extract the angular distribution
   into `reactions/<case>/angular_distribution.txt`.
2. **GRAZING**: run `GRAZING_9/grazing_exe/grazing_9` per beam energy
   (LEG 3 for the excitation-energy distributions), integrate over the
   E* window of interest and write
   `reactions/<case>/cross_sections.txt` (only needed for thick-target
   simulations).
3. Create one macro per target backing in `reactions/<case>/macros/`
   (mind the inverse-kinematics recoil/ejectile convention and
   `/reaction/excitationEnergy`), then run:
   `./run_docker.sh all reactions/<case>/macros/<backing>.mac`.

### Required Input Data Files
The primary generator relies on external ASCII data files to accurately sample reaction probabilities and kinematics. They are resolved **relative to the working directory of the simulation**, i.e. one level above `outputs/` inside the case folder (the `../` below):
- `angular_distribution.txt`: Contains Center-of-Mass differential cross-sections ($d\sigma/d\Omega$) versus the scattering angle $\theta$ (in degrees). Often imported from FRESCO/DWBA outputs.
- `cross_sections.txt`: Total cross-section data (in mb) versus beam energy (in MeV). This is vital for correctly weighting the interaction probability at different depths as the beam loses energy in a thick target.

## Acknowledgments & Disclaimer

This project was developed with the assistance of **Antigravity**, an agentic AI coding assistant developed by Google DeepMind.
