# Native FEM performance execution — 2026-09-30

**Historical screening record.** The withholding decision below is superseded
by the researcher's implementation directive. Harmonic consolidation,
prescribed physical quadrature and quadrilateral PML are now integrated as
described in [README.md](README.md). New feature-execution plots and timings
are kept separately in `.linecablemodels/fem/native-feature-review`; the original
measurements and numerical discrepancies below remain intact.

**Researcher review is pending.** The pass/fail labels below describe the
previously prescribed numerical checks, not a scientific acceptance decision.
The researcher decides whether the observed changes are acceptable. The earlier
decision to exclude candidates without first providing plots was premature.
Use [plot_results.jl](plot_results.jl) to inspect the retained raw results with
the public plotting API; this reads files and launches no FEM computations.

**Faster candidates exist; their scientific acceptance is pending your review.**
The earlier delivery added binary MSH 4.1 and optional native MUMPS ordering
and sparse-row preallocation in Julia and detached ONELAB. At that stage the
assembly and quadrangle candidates remained manual experiments. They have since
been integrated: the current implementation and runner comparison are described
in [README.md](README.md). The measurements below refer to the original frozen
comparison, not to the subsequent feature execution.

## Delivered changes

- Binary MSH 4.1 in managed meshing and exported geometry. Native round trips
  preserve every node coordinate, element connectivity and physical-group
  membership. All three retained pilot solves produce exactly identical raw
  matrices in ASCII and binary. Existing ASCII meshes remain readable and
  eligible for the existing mesh cache; there is no second cache.
- `mumps_ordering=nothing` and `petsc_prealloc=nothing` preserve GetDP defaults.
  Explicit values are forwarded to `-mat_mumps_icntl_7` and `-petsc_prealloc`.
  Export accepts `solver_options=(mumps_ordering=0, petsc_prealloc=256)` and
  presents both controls in the ONELAB Numerics group. The native runtime
  remains independent of Julia. Available ordering libraries depend on the
  GetDP build; production performs no search, acceptance or retry.
- Existing factorization reuse across source columns and mesh reuse are
  already implemented. The runner's explicit `mesh_policy=:remesh` remains
  untouched. Solver options participate in the existing run identity.

Illustrative explicit controls, **not a recommendation that AMD/256 is faster**:

```julia
compute(problem, formulation;
    options=(mumps_ordering=0, petsc_prealloc=256))
export_data(:onelab, problem, formulation; file_name="study/study.pro",
    solver_options=(mumps_ordering=0, petsc_prealloc=256))
```

Restart Julia to load changed extension code. Existing detached bundles retain
their own copied sources; re-export to obtain the new controls and mesh format.

## Measured costs

Five warmed Gmsh read/write repetitions, with the first call excluded:

| Retained mesh | ASCII write/read (s) | Binary write/read (s) | ASCII/binary bytes |
|---|---:|---:|---:|
| Aerial, 0.1 Hz, Γ=0 | 0.1422 / 0.1043 | 0.01645 / 0.01689 | 8,863,463 / 9,072,950 |
| Mixed, 0.1 Hz, Γ=0.99γearth | 0.1417 / 0.1042 | 0.01906 / 0.01842 | 8,807,532 / 9,030,186 |
| Mixed, 1 MHz, Γ=0.99γearth | 0.1225 / 0.1121 | 0.01846 / 0.01795 | 9,295,200 / 9,534,642 |

This saves approximately **0.20–0.21 s per Gmsh write/read pair**, not seconds
of factorization. Binary files are slightly larger here. No large end-to-end
speedup follows from this result.

The combined assembly candidate was measured in two warmed, serial repetitions
with execution order reversed. Each workload consists of all three pilot
solves, including both source columns:

| Native workload | Reference | Combined candidate |
|---|---:|---:|
| Median total GetDP elapsed | 94.300 s | 89.275 s |
| Mean summed assembly | 66.645 s | 61.763 s |
| Mean summed solve | 19.199 s | 19.185 s |
| Maximum observed peak RSS | 1,917,196 KiB | 1,916,352 KiB |

The **5.33% total improvement has not been integrated into production** because
expanded fixtures exceed the numerical screening thresholds. Acceptance remains
the researcher's decision. DOF counts and meshes are identical. Native process
timings exclude Julia compilation, mesh preparation and analytical comparisons;
they are not timings of the user's complete parametric study. These are two
matched repetitions, not a broad performance characterization.

Reference pilot meshes were reused, so meshing cost is zero in those timings.
The first quadrangle conversion recorded 1.048 s preparation including 0.575 s
Julia compilation, plus 0.123 s mesh writing, separately from GetDP. Fresh
fixture mesh preparation records `mesh_seconds` and `compile_seconds` in each
`mesh.toml`; it is not included in the pilot speed comparison.

Native AMD ordering was also screened on the same three meshes. Baseline/AMD
total times were 27.74/28.19, 33.04/33.03 and 33.80/34.09 s. All passed the pilot
preservation checks, but there is no measured speed advantage here and no
default change. Sparse preallocation is an exposed native control, not a claimed
optimization with a measured preferred value.

## Numerical candidates and screening flags

The fixed rule was at most 2% change in **every nonzero raw R/X/G/B entry**,
exact zero preservation and no sign reversals. No absolute floor, clipping or
analytical substitution was used.

**Structured PML quadrangles.** Pairing existing triangles along their grid
diagonals retained all coordinates, physical elements, shared boundaries and
voltage paths. Aerial 0.1 Hz unknowns fell from 359,568 to 275,018. Initial
GetDP times fell from 27.74 to 17.17 s with 9-point Gauss–Legendre integration;
16 points took 23.05 s. These are initial screening times, not warmed claims.
However G21 changed from −2.32662e−25 to −2.87486e−25 S/m (23.56%); the 16-point
rule still changed it by 23.47%. No signs reversed. The analytical value is
−5.95373e−25 S/m, but moving closer does not establish convergence. Both failed
the pilot gate, so neither entered expanded qualification or production.

**Combined harmonic assembly and physical three-point quadrature.** The former
combines σ and jωε on identical basis supports; the latter is polynomial-exact
in exact arithmetic on affine first-order physical triangles, retaining the
12-point PML rule and existing line rule. All 15 air/soil/mixed cases passed
(0.1 Hz, 1 kHz, 1 MHz at Γ=0; low/high endpoints at Γ=0.99γearth). The maximum
relative component change was 0.00371%, with no sign reversals. All warmed
pilot results also passed.

Expanded Γ=0 endpoints gave:

| Fixture | 0.1 Hz | 1 MHz |
|---|---|---|
| Three conductors in different layers | Pass | Pass |
| 49-wire screen with foil | Fail | Fail |
| Tubular sheath | Fail | Pass |
| Three sectors with wire neutral | Pass | Pass |

For the low-frequency screen, combined G32 changes from +3.54802e−28 to
−3.10004e−28 S/m; three conductance entries reverse sign. At 1 MHz screen G31
changes from −2.12288e−14 to −2.61125e−14 S/m (23.0%). Tube G21 at 0.1 Hz changes
from −3.22104e−31 to −2.00420e−30 S/m. These are raw preservation failures;
they do not by themselves establish a physical error of that relative size.

The two changes were then separated on the **same saved screen mesh/reference**.
Harmonic assembly alone produced four sign reversals; physical three-point
quadrature alone produced three. Both therefore stopped at that fixture.
This establishes that neither independently satisfies the locked raw-component
criterion. It does not identify whether the last-bit differences originate
principally in assembly, factorization or cancellation in extraction; no such
diagnosis is claimed. No physical-quadrature option or new weak terms were
promoted merely because the bare-wire pilots passed.

## Evidence and verification

All runs were serial under the main agent, with native output appended to
`.linecablemodels/fem/native-performance-20260930/live.log`. The same directory
retains the 24-file source snapshot, mesh/source hashes, native matrices,
analytical comparison columns for bare-wire cases, timings and completion files.
The original studies and their logs were not rewritten.

- `*/comparison-*.toml`, `components-*.csv`: pilot and spectrum comparisons.
- `fixtures/*/comparison.toml`, `components.csv`: original combined shape results.
- `fixtures/*/comparison-harmonic.toml`, `comparison-physical3.toml`: separation.
- `warmed-costs.csv`, `warmed-assessment.toml`, `warmed-preservation.toml`:
  matched native costs and retained-result checks.
- `mesh-io/results.toml`: lossless round trips and five-repetition I/O medians.

Focused regression: **557 passed, 2 failed, 0 errors**, across 12 maintained
test items. Export ownership, binary mesh loading, native solver controls,
managed/detached numerical parity and relocation, import, worker recovery and
option normalization passed. The two failures were both reproduced separately
with the frozen pre-change production sources in `/tmp/fem-native-before`:

1. `test/extensions/fem_export.jl:257`: exact string equality between tables from
   a frequency scan and a separately remeshed/repeated native solve. The
   pre-change test is at line 250. The assertion and its comparison policy were
   not modified.
2. `test/extensions/fem_resume.jl:88`: the test expects an unchanged mesh
   fingerprint after doubling evaluated earth permittivity. The same two
   differing fingerprints occur before and after these changes. The assertion
   was left intact.

Reproduce the focused suite (ONELAB's scan test requires local socket access):

```bash
env OPENBLAS_NUM_THREADS=1 OPENBLAS_DEFAULT_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  julia --project=test test/runtests.jl unit/engine/options.jl \
  fem_solver_controls fem_export fem_resume fem_import fem_workers
```

No scientific tolerances or test assertions were relaxed to obtain a pass.
`git diff --check` passes. The frozen source comparison confirms that the PDE,
integration rules, geometry construction and customary manual runner are
unchanged. All qualification processes have finished.

Production ownership: [options](../../../../src/engine/options.jl),
[managed command](../../../../ext/LineCableModelsGmshExt/workers.jl),
[mesh format](../../../../ext/LineCableModelsGmshExt/mesh.jl),
[export](../../../../ext/LineCableModelsGmshExt/export.jl), and
[detached native options](../../../../ext/LineCableModelsGmshExt/getdp/onelab.pro).

Native features used: [GetDP integration objects](https://getdp.info/doc/texinfo/getdp.html#Integration),
[Gmsh MSH format and structured meshing](https://gmsh.info/doc/texinfo/gmsh.html),
[PETSc MUMPS controls](https://petsc.org/release/manualpages/Mat/MATSOLVERMUMPS/).
