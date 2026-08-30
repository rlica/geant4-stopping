# Geant4 Stopping Simulation

## Introduction
The `geant4-stopping` project is a Geant4-based application designed to simulate nuclear transfer reactions (e.g., 18O + 30Si -> 16O + 32Si*) and the subsequent slowing-down process of heavy ions in a solid target medium. It employs a custom minimal physics list tuned for high-accuracy, low-energy electromagnetic interactions (including atomic deexcitation, fluorescence, and Auger electron emission) alongside step limiters to ensure precise $dE/dx$ integration. 

The primary objective is to simulate the primary beam's energy loss in a "thick" target, dynamically model the interaction point using energy-dependent cross sections, and track the kinematics, stopping time, and stopping depth of the recoiling ions. It incorporates angular distributions from FRESCO calculations to sample realistic ejectile scattering angles in the Center-of-Mass (CM) frame, adjusting reaction Q-values for specific excitation energies.

## Quickstart

### Prerequisites
- **Geant4** (version 11.x or similar) with appropriate physics datasets.
- **CMake** (version 3.16 or higher).
- **Python 3** with `matplotlib` and `pandas` (for plotting results).

### Building the Project
To compile the simulation, you can use the provided bash script or build it manually using CMake.

**Using the provided script:**
```bash
./run_geant4_stopping.sh
```
This script sets up the environment, builds the executable, runs the default macro, and generates PDF plots automatically.

**Manual Build:**
```bash
mkdir build && cd build
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
