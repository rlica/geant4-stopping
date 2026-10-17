# Geant4 Stopping Simulation

## Introduction
The `geant4-stopping` project is a Geant4-based application designed to simulate nuclear transfer reactions (e.g., 18O + 30Si -> 16O + 32Si*) and the subsequent slowing-down process of heavy ions in a solid target medium. It employs a custom minimal physics list tuned for high-accuracy, low-energy electromagnetic interactions (including atomic deexcitation, fluorescence, and Auger electron emission) alongside step limiters to ensure precise $dE/dx$ integration. 

The primary objective is to simulate the primary beam's energy loss in a "thick" target, dynamically model the interaction point using energy-dependent cross sections, and track the kinematics, stopping time, and stopping depth of the recoiling ions. It incorporates angular distributions from FRESCO calculations to sample realistic ejectile scattering angles in the Center-of-Mass (CM) frame, adjusting reaction Q-values for specific excitation energies.

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
   # or explicitly with a specific macro:
   ./run_docker.sh all ../run.mac
   ```

3. **Other convenient Docker commands:**
   ```bash
   ./run_docker.sh compile          # Compile C++ simulation inside Docker
   ./run_docker.sh run ../30Si_run.mac # Run a specific macro
   ./run_docker.sh plot             # Generate/update histogram plots
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
Execute the simulation by providing a macro file:
```bash
./stopping ../run.mac
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
The physical target and reaction kinematics are highly configurable via custom UI commands, typically defined in your `run.mac` file.

#### 1. Target Configuration
You can define a custom, homogeneous stopping medium without recompiling.
- `/detector/customZ <Z>`: Set the Atomic Number (e.g., `14` for Si).
- `/detector/customA <A>`: Set the Mass Number (e.g., `30`).
- `/detector/customDensity <density> <unit>`: Set the mass density (e.g., `1.25 g/cm3`).
- `/detector/customSymbol <Symbol>`: Set the chemical symbol (e.g., `Si`).

*Note: The target is built dynamically using these properties. Changes must be executed before `/run/initialize`.*

#### 2. Reaction Kinematics
- `/reaction/thinTarget <true/false>`: 
  - `false` (default): Simulates beam slowing down inside the target. The interaction depth is sampled using the beam's energy loss and external energy-dependent cross-sections.
  - `true`: Bypasses beam slowing down, forcing the reaction to occur directly at the surface of the target.
- `/reaction/targetZ <Z>` and `/reaction/targetA <A>`: Target nuclide involved in the reaction.
- `/reaction/recoilZ <Z>` and `/reaction/recoilA <A>`: Recoiling heavy fragment.
- `/reaction/ejectileZ <Z>` and `/reaction/ejectileA <A>`: Light ejectile fragment.
- `/reaction/thetaMin <angle> <unit>` and `/reaction/thetaMax <angle> <unit>`: Limits the Center-of-Mass scattering angle (e.g., `0 deg` to `180 deg`).

#### 3. Primary Beam Generation
The simulation uses the standard `G4ParticleGun` to define the incoming beam.
- `/gun/particle ion`: Set particle type to ion.
- `/gun/ion <Z> <A>`: Define the beam ion (e.g., `/gun/ion 8 18` for 18O).
- `/gun/energy <energy> <unit>`: Set the initial beam kinetic energy (e.g., `/gun/energy 26 MeV`).

### Required Input Data Files
The primary generator relies on external ASCII data files to accurately sample reaction probabilities and kinematics. These should be present in the project's root folder (one level above the `build` directory):
- `angular_distribution.txt`: Contains Center-of-Mass differential cross-sections ($d\sigma/d\Omega$) versus the scattering angle $\theta$ (in degrees). Often imported from FRESCO/DWBA outputs.
- `cross_sections_18O_30Si_2n_Estar_5p5_6p5.txt` (or equivalent): Total cross-section data (in mb) versus beam energy (in MeV). This is vital for correctly weighting the interaction probability at different depths as the beam loses energy in a thick target.

## Acknowledgments & Disclaimer

This project was developed with the assistance of **Antigravity**, an agentic AI coding assistant developed by Google DeepMind.
