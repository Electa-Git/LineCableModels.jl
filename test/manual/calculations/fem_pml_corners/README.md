# PML corner qualification

The anisotropic production option was removed at the user's request on
2026-09-30: its measured speed gain did not justify the added complexity and
conductance degradation. Production now uses structured PML corners only.
The qualification scripts, [feature-results.md](feature-results.md) and saved
run data remain historical evidence, not supported feature instructions.

This is a serial Julia experiment, not production code. Original numerical
results are retained. Evidence and the append-only log are under
`.linecablemodels/fem/pml-corner-mesh/`.

The frozen meshes use localized physical-domain refinement. Both baseline and
candidate solves use the current corrected quasi-fw sources. Only the four
corner interiors change. Physical elements, noncorner PML elements and native
voltage-path nodes must have identical hashes before a candidate is solved.
Each corner is meshed in an isolated native Gmsh model with its frozen boundary
segments, then inserted into the reference mesh. This avoids remeshing any
physical material. It is qualification tooling, not a proposed runtime pipeline.

The proposed anisotropic metric uses the existing frequency/material/Gamma,
clearance, cubic stretch and coefficient-variation density in each direction.
The density in x is multiplied by the envelope of the finite-layer modal
amplitude in y to power 2/3, and conversely. Modes include the static solution;
the amplitude envelope is floored at 1/8, limiting this cross-coordinate
relaxation to a factor of four in element length. This is an interpolation
heuristic for qualification, not an error estimator for cable quantities.
Gmsh BAMG consumes a tensor background view. The two prescribed metric-length
multipliers are 1.0 and 1.35; they are fixed before any candidate solve.

Pilots: air/air 0.1 Hz Gamma=0, mixed 0.1 Hz and 1 MHz Gamma=0.99 gamma_earth.
All use soil rho=0.1 ohm m and wire radius=0.0425 m. Raw per-entry R/X/G/B
changes must stay within 2%, with no new signs. Exact zero references require
exact zero candidates; no clipping/floors enter comparisons. Analytical values
are reported separately. The runtime gate is 20% lower median total warmed
elapsed time over matched pilot workloads, including mesh construction.

Only a surviving candidate proceeds to the existing fifteen placement/Gamma/
frequency cases and low/high three-wire, screen, tube and sector fixtures.
Only after qualification may a feature option be added for managed Julia and
native detached ONELAB, including cache identity and maintenance tests.
Failure of both candidates ends this approach without production edits.

Run from the repository root:

```bash
JULIA_DEPOT_PATH=/tmp/lcm-onelab-depot:/home/amartins/.julia \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/home/amartins/.juliaup/bin/julia --compiled-modules=existing --startup-file=no \
  --project=test test/manual/calculations/fem_pml_corners/qualify.jl \
  >> .linecablemodels/fem/pml-corner-mesh/live.log 2>&1
```

```bash
tail -n 80 -F .linecablemodels/fem/pml-corner-mesh/live.log
```

`mesh-only` prepares candidates without invoking GetDP. Completed meshes and
solves are resumed after source/mesh hash verification. First-call Julia
compilation is recorded separately from native GetDP time.

The follow-up `diagonals_and_soil.jl` separates two questions. It reverses only
the structured corner diagonals on the saved aerial 0.1 Hz mesh, retaining all
nodes and element counts. Independently, it compares a fresh detached mixed
0.1 Hz export with and without soil localization, keeping air localization
enabled on both. The fresh bundle includes all ten manual-runner frequencies
for inspection, but this focused soil comparison solves only 0.1 Hz. It does
not overwrite `.linecablemodels/fem/onelab-two-bare-wires/`.

```bash
JULIA_DEPOT_PATH=/tmp/lcm-onelab-depot:/home/amartins/.julia \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/home/amartins/.juliaup/bin/julia --compiled-modules=existing --startup-file=no \
  --project=test test/manual/calculations/fem_pml_corners/diagonals_and_soil.jl
```

This script appends native progress to the same `live.log` itself. Run
`diagnostic_images.jl` with the same Julia environment to render its saved
soil meshes without meshing or solving again.

The feature-only `verify_feature.jl` and `benchmark_feature.jl` callers were
removed together with the production option. Their completed measurements
remain in the results notes and saved run data.

`verify_runner_meshes.jl` checks all meshes of the current two-wire manual
runner with structured corners, retaining its frequency order and conductor CAD reuse, without running
the solver grid or opening plots. It covers the fifth-frequency failure that
the earlier two-frequency check missed, and appends to the same live log.
