# Native voltage cleanup: execution evidence

Execution date: 2026-09-28. This records the completion work authorized after
the [native voltage description and plan](fem-native-voltage-extraction-plan.md).
The existing production implementation already used native measurements; this
work removes its stale consumers and closes the specified validation gaps.

## Result and scope

Julia builds measurement geometry with the rest of the model, meshes it with
Gmsh and invokes GetDP directly. The shared GetDP formulations evaluate the
scalar reference and, for quasi-full-wave, integrate the transverse vector
potential. No Julia/Python triangle clipping, serialized quadrature, generated
`paths*.pro`, or measurement-preparation subprocess remains in that route.
The detached export uses the same formulations and runs with Gmsh/GetDP only.

The ordinary `_prepare_run_inputs!` remains: it writes model inputs, snapshots
and solver sources. Julia's matrix reduction/inversion remains separate from
voltage extraction. No production equations, runtime controls, tolerances,
mesh/run identity versions or performance policy changed in this cleanup.

Retired executable paths:

- `run_fem_voltage_reference.jl`: removed its old solve/instrumentation modes;
  retained only `--plot STUDY_DIRECTORY`, with the original plotting function.
- `run_fem_terminal_extraction.py`: removed the solve mode, deleted-module
  imports and measurement preparation. `--collect` reads existing raw tables
  and command records and reconstructs the historical summary.
- `fem_static_boundary/mesh.jl`: removed its five-argument preparation mode;
  ordinary native meshing remains.
- `fem_static_boundary/run.py`: removed the cable comparison that required
  preparation. The independent annulus controls remain runnable.
- `fem_cartesian_compactification/run_cable.py`: removed the obsolete launcher.
  Independent wave experiments and their files remain untouched.

Current documentation points to native extraction tests. Historical studies
retain their numbers and old averaging/reference labels, with explicit notices
that those conventions describe the retired implementation. Existing exported
studies were not overwritten. Unrelated dirty-worktree changes were preserved.

## Validation

The active package and extension paths were verified in a fresh Julia process.
Environment: Julia 1.12.7, Gmsh 4.15.2 and the backend's native GetDP 3.5.0
artifact, with Gmsh support:

```text
/home/amartins/.julia/artifacts/51e049a32eeb2ffebcc8a5945bb00f70470a3d29/getdp-3.5.0-Linux64/bin/getdp
```

Focused suite command:

```sh
JULIA_DEPOT_PATH=/tmp/lcm-onelab-depot:/home/amartins/.julia \
julia --project=test --startup-file=no test/runtests.jl \
  fem_native_measurements fem_quasi_full fem_mesh_grading \
  fem_workers fem_resume fem_export
```

Result: **428/428 assertions passed**, 14 maintained test items in six files,
541.05 s elapsed. Counts by file: export 165, mesh grading 45, native
measurements 46, quasi-full/extraction 59, resume 80 and workers 33.
The first combined execution tool reported exit 143 after the complete passing
summary and the runner's `run completed` message. To check shutdown explicitly,
the same tests were repeated in two processes, capturing each shell exit code:
`fem_quasi_full` passed 59/59 in 400.70 s; the other five files passed 369/369
in 221.26 s. Both exited with status 0. The worker cancellation test also passed
33/33 with status 0 on its own. The initial tool termination remains unexplained;
the repeated checks cover all 428 assertions and confirm clean process exits.

| Plan gates | Evidence |
|---|---|
| V01, V10: removal | Source audit has no live preparation helper or deleted-module import. The only executable `PathDataPath` mention is a negative assertion. Fresh direct runs for both physics choices have ordinary input snapshots and native raw columns, with no generated path files or weight arrays. |
| V02: geometry | Existing overhead/buried tests plus the extended rotated-polygon, insulated, disconnected-terminal case check the chosen endpoint, reference and absence of volume embedding. |
| V03, V04: integration/support | Existing exact complex open/closed-path integrals retain `atol=1e-13`. The new piecewise edge-field fixture checks reversal, zero metal contribution and zero lookup outside the field support; native volume-element lookup identifies the outside point as invalid geometry. |
| V05: production extraction | Both actual formulation extraction blocks pass their manufactured, basis-dependent complex reference/vector checks with non-unit source amplitudes and field maps disabled, at the existing `rtol=1e-11`. |
| V06: route preservation | Direct Julia and relocated native export compare primitive P/Z and final Y/Z for both formulations, both bases and 50 Hz/10 kHz on matching meshes. The existing componentwise bound is unchanged: `2e-9*abs(reference) + 100eps(Float64)*maximum(abs, reference_component)`. |
| V07: independent refinement | The new manufactured fixture and the separate physical toy below refine only line elements while requiring identical volume coordinates/connectivity. |
| V08, V09: execution | Existing mesh, restart, worker, factorization and field-map checks remain in the focused suite; an additional geometry assertion rejects a mesh missing its native measurement group. |
| V11: scientific scope | One native path per terminal remains distinct from the historical overhead contour average. Refinement sensitivity below is not a physical/PML accuracy bound. |

The new manufactured fixture is in `test/fixtures/data/fem/native_getdp/support.*`.
In a unit square, four horizontal intervals have vector fields
`(1+2j)*(0,2,0)`, zero metal support, `(1+2j)*(0,4,0)` and
`(1+2j)*(0,6,0)`. Its exact upward integral is `3+6j`. It tests 1, 2, 4 and 8
line subdivisions: the single-element error exceeds `1e-2`; the final two
aligned meshes meet absolute error `1e-13`, and reversal changes the sign.

This analytic fixture sets its known fractional line-node coordinates exactly
after native meshing. Gmsh's transfinite placement was about `1e-12` away from
the prescribed material jumps, which would otherwise test coordinate placement
against a `1e-13` integration bound. Only test-owned standalone line nodes are
set; volume coordinates/connectivity are checked unchanged. No tolerance was
relaxed and no production interpolation/preparation helper was introduced.

## Independent physical line refinement

A fresh native export used 5 mm wires at `(0,0.1)` and `(0.2,-0.1)` m,
earth resistivity 100 ohm m, relative permittivity 10, eight PML layers and
no bundle/Kron/transposition reduction. At 10 kHz, each native line-refinement
level was meshed and both formulations were solved with both source columns
and maps disabled. Volume coordinates and connectivity were exactly identical
at every level; Z was bitwise identical across levels for each formulation.

The table reports `max(abs((M_level-M_5)/M_5))` over individual complex matrix
entries, with elementwise division, separately for P and Y. Level 5 is the
finest tested native line mesh, not an exact reference.

| `VoltageRefinements` | Line elements | Quasi-full P difference | Quasi-full Y difference |
|---:|---:|---:|---:|
| 0 | 337 | 1.141546% | 0.797239% |
| 1 | 674 | 0.090732% | 0.062876% |
| 2 | 1348 | 0.310660% | 0.216759% |
| 3 | 2696 | 0.060693% | 0.042433% |
| 4 | 5392 | 0.089091% | 0.062183% |
| 5 | 10784 | reference | reference |

Quasi-TEM P and Y were bitwise identical at all levels. Its extraction uses a
scalar difference, so measurement-line quadrature does not enter that branch.
The quasi-full changes are nonmonotone. The eight-layer PML and fixed volume
mesh are deliberately small execution fixtures; neither this table nor route
parity establishes engineering accuracy. No performance improvement is claimed.

To reproduce native solves on a fresh export, use these commands for each level
0 through 5, preserving each level's mesh/results before the next run:

```sh
gmsh study.geo -setnumber BuildMesh 1 -setnumber FrequencyIndex 2 \
  -setnumber VoltageRefinements LEVEL -0
getdp study.pro -msh study.msh -solve LineCableModelsFEM \
  -setnumber FrequencyIndex 2 -setnumber Physics 1 -setnumber PlotFieldMaps 0
```

Use `Physics 0` for quasi-TEM. The evidence script additionally verifies exact
volume equality and records the native commands, matrices and line counts.

## Retained research tools and evidence

The saved-output collector reproduced 72 complex values in each historical
formulation summary, 144 total, with maximum difference exactly zero. The saved
voltage-reference plot regenerated PNG/PDF and was visually inspected; its
historical labels and signed conductance scale remain unchanged. These checks
used copies of the saved data.

All four independent static annulus solves passed their original assertions.
Maximum field errors at bulk sizes 0.2 and 0.1 were respectively
`0.01419280718` and `0.00881117316` for both modes. The scalar-field bound
`1e-3` also passed. Modified manual Julia files passed `Meta.parseall`; edited
Python files passed AST parsing. `test/manual/**` remains excluded from normal
test discovery. Relative document file links were checked.

Local evidence is under
`.linecablemodels/fem/native-measurement-cleanup-evidence/`: the pre-change
status/HEAD, per-file snapshots and hashes, task-only diff and manifest,
source audit, fresh Julia input/raw-column copies, focused suite log,
reader/plot outputs, independent boundary controls and `refinement/` meshes,
matrices and sensitivity CSVs. `refinement.jl` records the executable study;
`refinement.log` records the identical volume digest at every level.

The baseline HEAD was `94e4b058706031e937954802cec3e2921c4dd751` with existing
uncommitted work. Evidence uses the active working tree, not HEAD alone. No
files were staged or committed. The full repository suite and documentation
build were not run for this bounded cleanup.
