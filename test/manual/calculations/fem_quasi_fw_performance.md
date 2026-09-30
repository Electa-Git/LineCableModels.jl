# Quasi-fw performance diagnosis — 2026-09-29

The user reports that the current quasi-fw two-wire results now meet their
correctness requirement. Preserve that calculation as the performance baseline.
Retiring quasi-TEM is conditional on improving quasi-fw performance; neither the
physics selector nor production defaults have been changed by this investigation.
No new scientific acceptance or refinement policy belongs in the engine.

## Measured current cost

The latest seven retained scans are `run-0D8nGb`, `run-DoZCtM`, `run-vEE4tu`,
`run-XhkyiL`, `run-mNELmx`, `run-DBjPXG`, and `run-Il3Sjn`, under
`.linecablemodels/fem/runs/`. They cover 70 frequency solves / 140 source columns.
The last three are the radius branch: their run durations total 1613.85 seconds,
consistent with the reported 1620.98 seconds for the enclosing parametric study.
Only about 41.34 seconds precede the first GetDP process across those three runs.
Thus meshing and Julia setup are not the dominant elapsed cost here.

Accumulated worker time across all seven scans (not elapsed scan time):

| Stage | Seconds | Share of worker wall time |
|---|---:|---:|
| Assembly | 5893.92 | 83.4% |
| Solve | 785.21 | 11.1% |
| Constraint updates | 118.22 | 1.7% |
| Output | 26.61 | 0.4% |
| Total worker wall time | 7070.25 | 100% |

There are ten factorizations for twenty source columns in each scan. The backend
already uses `GenerateRHSGroup` and `SolveAgain` for the second source. There is
no missing factorization reuse to recover between the two conductors. Each
frequency currently has a different mesh/material matrix, so factorization reuse
across frequencies would require a different algorithm.

## PML cost

On the retained 100 ohm m, 4.25 cm mixed wire case (`run-vEE4tu`), the 0.1 Hz mesh
has 400128 PML triangles, 33298 physical air/earth triangles and 626 conductor
triangles. At 1 MHz it has the same 400128 PML triangles, 29412 physical
medium triangles and 6004 conductor triangles. Thus roughly 92% of its triangles
are in the PML.

The four tensor-product corners alone contain `4*2*192^2 = 294912` triangles.
Their count is fixed by `pml_layers`, not by the physical domain dimensions or
ordinary mesh-size factors. Halving `domain_skin_depths` therefore does not
remove that dominant fixed cost. Coarsening the ordinary physical mesh also
risks undoing local accuracy for comparatively little benefit. A reduced PML
layer count, or separately designed side/top/bottom grading, is a possible
second-stage experiment, not a presently qualified replacement.

Native structured constraints override ordinary mesh-size prescriptions:
[Gmsh transfinite meshing](https://gmsh.info/doc/texinfo/#Structured-grids).

## Same-mesh coefficient pilot

Evidence and the reproducible runner are in
`.linecablemodels/fem/performance-pilot-20260929/`. It copies the saved GetDP
inputs, references the same retained meshes and runs serially. Production files,
mesh sizes, quadrature, equations, boundary conditions and voltage paths stay
unchanged. The candidate evaluates the stretch once per coordinate in the PML
tensors using native registers, and selects the exact identity-stretch material
coefficients in physical regions. This uses native
[GetDP expressions and registers](https://getdp.info/doc/texinfo/getdp.html#Run_002dtime-variables-and-registers).

Both endpoint pairs completed, with both sources at each frequency:

| Frequency | Baseline wall [s] | Candidate wall [s] | Wall reduction | Baseline assembly [s] | Candidate assembly [s] |
|---|---:|---:|---:|---:|---:|
| 0.1 Hz | 243.261 | 159.187 | 34.56% | 202.679 | 118.514 |
| 1 MHz | 239.550 | 155.481 | 35.09% | 202.899 | 118.947 |

Every saved Z and P component is identical; inversion consequently produces
identical Y, including the individual real and imaginary parts and their signs.
At 0.1 Hz the smallest conductance entry, -8.635583649900808e-22 S/m, is
unchanged. This is 1.53–1.54x in the two native pairs. It is not a measured
full Julia-study speedup or a repeatability result. Timings on the current host
are slower than the earlier scan, so only the paired pilot is compared; no
cross-run extrapolation is used. These are existing native binaries with no
Julia compilation in the timed process.

All four solves and their assessment are complete. `comparison.csv` records
zero change in every compared Z/P/Y entry and no R/X/G/B sign changes.
`run.py` skips completed cases; `assess.py` reproduces the comparisons.
The first exploratory candidate had a GetDP piecewise-function
syntax error before assembly; its log is retained as `parse-failed.log`. The
corrected candidate uses `[All]` for the default piecewise coefficient.

## Order of further work

1. The endpoint pilot and full public-route comparison are complete. Production
   integration is recorded below.
2. The measured four-worker/one-thread setting is selected in the manual runner.
   Assembly dominates, and more BLAS threads alone cannot accelerate its
   interpreted coefficient evaluation. The package defaults remain two workers
   and one thread; the benchmark does not justify a host-independent increase.
3. Next compare PML grading/layer prescriptions
   while preserving the conductor and near-interface mesh. Use raw individual
   G and B entries and signs, not only relative complex Y, across the current
   air/buried/mixed fixtures and resistivities. Extending to 100 MHz needs its
   own upper-frequency check; the current retained grid ends at 1 MHz.

The next bounded selection completed in
`.linecablemodels/fem/performance-public-20260929/`: three serial public scans of
the retained ten-frequency, 100 ohm m mixed-wire problem. They compare baseline
2 workers/4 threads, candidate 2 workers/4 threads, and candidate 4 workers/1
thread. Each scan uses the same declared mesh controls and permits mesh reuse,
but executes fresh solves. Each result retains native timings, public matrices,
raw primitive columns and first-call compilation time. The final assessment
requires identical mesh hashes and records every R/X/G/B difference and sign.
The first public call and subsequent warmed calls are labelled separately;
first-use compilation must not be credited to the coefficient optimization.

Those scans used test-only executable wrappers to select frozen baseline/candidate
solver copies, before production integration. The launcher locks against duplicate
batches, skips completed scans and automatically assesses the outputs. It runs
the three scans in sequence; frequency concurrency is the measured engine
option, with only the main implementation agent involved.

```sh
tail -F .linecablemodels/fem/performance-public-20260929/progress.log
```

All three public scans finished on 2026-09-29:

| Variant | Public wall [s] | Compilation [s] | Accumulated assembly [s] |
|---|---:|---:|---:|
| Baseline, 2 workers / 4 threads | 1254.381 | 66.948 | 1960.866 |
| Coefficient candidate, 2 workers / 4 threads | 650.113 | 0 | 941.221 |
| Coefficient candidate, 4 workers / 1 thread | 383.842 | 0 | 883.861 |

The assessment verified identical mesh hashes across all three. Candidate
2/4 has identical Z and Y entries. Candidate 4/1 has maximum relative G-entry
change 2.4644e-6 (0.00024644%), maximum absolute Y-entry change 4.3548e-13 S/m,
and no R/X/G/B sign changes. The warmed 650.113-to-383.842 comparison reduces
wall time by 41.0%; the first baseline includes compilation and must not be
presented as an entirely warmed speedup. These measurements cover one
ten-frequency mixed-wire case. These results were obtained with the copied
prototype before production integration.

The user's separate read-only PML-resolution inquiry is recorded in
`fem_pml_resolution_diagnostics.md`. All ten saved 100 ohm m meshes contain
400128 PML triangles; its CSVs separate prescribed counts, measured geometry,
historical qualification evidence and untested count projections. No PML
resolution change or new solve was made for that inquiry.

## Production integration — 2026-09-29

`ext/LineCableModelsGmshExt/getdp/pml.pro` now contains exactly the qualified
candidate expressions, with explanatory comments. Native registers reuse sx/sy
within each tensor evaluation. Material coefficients in the physical air, earth,
conductors and passive regions use the exact identity-stretch expression. PML
coefficients, profiles, meshes, quadrature, equations, voltage extraction and
boundary conditions retain their benchmarked values.

`test/manual/calculations/run_two_bare_wires_fem.jl` now uses
`frequency_workers=4, solver_threads=1`. Its mesh controls are unchanged.
The ordinary backend/export source capture carries the optimization into newly
created runs and detached exports; there is no runtime prototype wrapper.
Existing detached bundles are saved copies and must be re-exported to receive
the new source. Start a fresh Julia session before rerunning the manual script:
the extension captures GetDP sources when it loads.

Evidence is in `.linecablemodels/fem/performance-integration-20260929/`.
`source-verification.json` confirms production coefficient expressions match
the benchmark candidate after excluding comments, and that `pml.pro` is the
only changed production file relative to the benchmark's src/ext manifest.
`checks.log` records the existing source/resume and detached export/parity
regressions: **179/179 assertions passed across three maintained test items**,
including native Julia/detached numerical parity at two frequencies for both
quasi-fw and quasi-TEM. The first fresh backend run's captured `pml.pro` also
matches the production file byte for byte. Manual-runner syntax and scoped
whitespace checks passed. No mesh-resolution experiment is included in this
integration.

Reproduction:

```sh
JULIA_DEPOT_PATH=/tmp/lcm-onelab-depot:/home/amartins/.julia \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
julia --compiled-modules=existing --startup-file=no --project=. test/runtests.jl \
  'resume requires effective inputs' \
  'detached native export and caller ownership' \
  'detached native numerical parity and relocation'
```

Removing the unused quasi-TEM branch would simplify maintenance; it would not
speed up a run that already selects quasi-fw. That retirement remains a separate
conditional decision.
