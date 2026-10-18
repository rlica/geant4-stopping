# Reaction case: 36S + 18O -> 38S + 16O (inverse kinematics, 2n transfer)

Beam: 36S at ~64 MeV on an 18O target. Inverse kinematics: the heavy 36S
beam impinges on the light 18O target, and the 38S recoil is tracked.

## Status / TODO

- [ ] Fresco calculation for this reaction (create `Frescox/s36-o18-2n/`
      and run frescox per beam energy).
- [ ] Extract the angular distribution:
      `python3 extract_fresco.py <fresco_out> angular_distribution.txt`
- [ ] Produce `cross_sections.txt` (total 2n cross section in mb vs beam
      energy in MeV, same 3-column format as the o18_si30_2n case).

## Contents (target layout)

- `angular_distribution.txt` - CM dSigma/dOmega from the new fresco run.
- `cross_sections.txt` - energy-dependent total cross section.
- `macros/` - one macro per target backing, e.g. `ta_backing.mac`.
- `outputs/` - CSV and plots produced by the simulation (gitignored).

## Run

From the repository root:

    ./run_docker.sh all reactions/s36_o18_2n/macros/ta_backing.mac
