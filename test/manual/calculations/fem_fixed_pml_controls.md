# Fixed-control directional PML meshing

This feature separates the requested normal interval counts from the grading
shape. It applies to the existing shared quasi-fw/quasi-TEM geometry and to
both Julia-managed computation and detached ONELAB export.

The 96-layer examples below are operating controls, not a qualified replacement
for 192 layers across all cases. The subsequent aerial radius sweep at 0.1 ohm m
exposed low-frequency mutual-conductance sign errors at `(96,96,96)` and
`(96,96,48)`. See [the follow-up diagnosis](fem_pml_low_frequency.md).

```julia
mesh_controls = (pml_layers=(96,96,48), pml_grading=(192/191)*log(1536))
result = compute(problem, fem; options=mesh_controls)
entry = export_data(:onelab, problem, fem;
    file_name="detached-study/study.pro", mesh_options=mesh_controls)
```

The tuple order is `(side, top, bottom)`. Scalars apply to all three directions.
Both controls normalize through `computation_options`; the existing mesh plan,
mesh fingerprint, run inputs, and mesh metadata retain the directional tuples.
Scalar inputs and equivalent tuples have the same identity. Changing either
control invalidates mesh reuse and completed-run resume through the existing
identity checks. An explicit `mesh_path` remains a caller-owned override.

For a thickness `L`, exponent `g`, and `N` intervals, nodes follow
`d_i = L*expm1(g*i/N)/expm1(g)`, with uniform spacing at `g=0`.
Native Gmsh transfinite curves receive `N+1` nodes and ratio `exp(g/N)`;
orientation reverses the ratio. One interval has two endpoints irrespective of
`g`. Exponents whose native ratio rounds to one produce uniform spacing.
Unsupported counts and unrepresentable spacing are construction errors.
Repeated assignments to a shared curve must agree after orientation correction.

Side counts are shared across both sides and both media. The bottom voltage
path uses exactly the same bottom sequence as its adjacent patches. Top
corners inherit side/top counts and bottom corners side/bottom counts; their
triangle count is `4*Ns*(Nt+Nb)`. There is no separate corner control.

The API default remains 128 intervals in each direction. The manual runner used
192 during this feature's qualification; its editable controls are forwarded to
its detached export as well. Its current selection is not an accuracy guarantee.
The fixed default exponent `(192/191)*log(1536)` retains the old normalized
192-interval shape. Other counts intentionally sample that shape instead of
retaining their old count-dependent shapes. Whole unstructured meshes are not
promised to be byte-identical.

Physical-domain dimensions, PML thickness and strength, tangential targets,
conductor prescriptions, materials, quadrature, solver settings, equations,
and complete scalar-plus-circulation voltage extraction are unchanged. All
nine GetDP source files were compared with the pre-edit snapshot. No scientific
acceptance test, adaptive refinement, extra solve, sign correction, or result
comparison was added to production.

## Implementation checks

Focused items cover option types and invalid inputs; native normal nodes at
1, 5, 8, 128 and 192 intervals; uniform/tiny exponents and reversed curves;
unequal directional counts and exponents; exact corner counts; shared field
nodes and the buried voltage path; scalar/tuple identity; metadata, cache and
resume; and exported native geometry. The existing detached numerical parity
item uses `(8,6,4)` counts and `(3,2,1)` grading in both physics, comparing the
individual real and imaginary components of Z, P and Y. That check uses the
same mesh to isolate equation/extraction parity; the separate geometry check
meshes the detached `.geo` and compares its curve nodes with the managed route.

All eight selected items passed: 576 assertions across the final runs. The
full repository suite was not run. Evidence accounting is in
`checks-complete.json`; earlier failed test expressions remain visible in the logs.

Reproduce from the repository root:

```bash
JULIA_DEPOT_PATH=/tmp/lcm-onelab-depot:/home/amartins/.julia \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/home/amartins/.juliaup/bin/julia --compiled-modules=existing --startup-file=no --project=. \
  test/runtests.jl \
  'Engine / FEM computation option ownership' \
  'Gmsh FEM / fixed PML' 'Gmsh FEM / directional PML' \
  'Gmsh FEM / exterior grading' \
  'Gmsh FEM / detached native export and caller ownership' \
  'Gmsh FEM / detached native numerical parity and relocation' \
  'Gmsh FEM / resume requires effective inputs'
```

Run evidence is retained in `.linecablemodels/fem/fixed-pml-controls/`.
`before/` contains the pre-edit files from the already dirty worktree; the
feature diff is relative to that snapshot, not to the branch's unrelated edits.
The first checks exposed two test-expression mistakes: a vector norm where the
criterion was per-node coordinate error, and a new inference assertion stronger
than the existing supplied-options API. The latter reproduced against the
pre-edit function (`inference.log`). The checks retain their original per-node
tolerance, verify concrete normalized tuple types, and preserve the existing
default-options inference assertion.

## Bounded observation

`run_fem_fixed_pml.jl` makes three independent public computations. Each has
exactly two frequencies (0.1 Hz and 1 MHz), two bare copper wires of radius
0.0425 m at `(0,+1)` and `(1,-1)` m, and 100 ohm m earth with relative
permittivity and permeability one. All three use quasi-fw, fixed default grading,
and the active physical/conductor preset. Counts are explicitly:

- `(192,192,192)`;
- `(96,96,96)`;
- `(96,96,48)`.

This is six frequency solves and twelve source columns. It does not expand the
search or select a setting. Resuming uses the existing compatible-run mechanism.
Raw signed R/X/G/B and complete complex P are recorded without clipping. P here
is the FEM inverse-admittance coefficient in ohm m, with `Y=inv(P)`.

```bash
JULIA_DEPOT_PATH=/tmp/lcm-onelab-depot:/home/amartins/.julia \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/home/amartins/.juliaup/bin/julia --compiled-modules=existing --startup-file=no --project=. \
  test/manual/calculations/run_fem_fixed_pml.jl \
  > .linecablemodels/fem/fixed-pml-controls/observation.log 2>&1
```

```bash
tail -n 60 -F .linecablemodels/fem/fixed-pml-controls/observation.log
```

The script writes `observation/observations.md`, measured `meshes.csv`, all-entry
`differences.csv`, and per-configuration raw `matrices.csv` and `cost.toml`.
Compilation is reported separately from the inclusive elapsed time; later calls
use warmed Julia code. Each native run retains its existing timing summary and
solver logs. No discretization/PML error bound on Z, P or Y is available from
these observations; native residuals retain their original meanings.

## Completed observation (2026-09-29)

All six prescribed frequency solves and twelve source columns completed. The
measurement files are under
`.linecablemodels/fem/fixed-pml-controls/observation/`.

| Counts (side, top, bottom) | PML triangles per mesh | Corner triangles | Scan wall, s | Compilation included, s | Slowest native process, s |
|---|---:|---:|---:|---:|---:|
| (192,192,192) | 400128 | 294912 | 107.015 | 32.898 | 70.963 |
| (96,96,96) | 126336 | 73728 | 25.414 | 0.201 | 22.996 |
| (96,96,48) | 103968 | 55296 | 21.344 | 0.000 | 19.199 |

PML and corner counts are identical at both sampled frequencies. Total triangle
counts are 434052/435538, 160262/161746, and 137894/139384 respectively at
0.1 Hz/1 MHz. Native timings exclude Julia compilation; the cold scan total
must not be interpreted as warmed execution. This is one scan per setting.

Compared with the 192-layer baseline, the largest individual G relative change
was 0.04323%; the largest absolute G changes were 2.042e-12 S/m at 0.1 Hz
and 9.338e-9 S/m at 1 MHz. The largest relative change among the separate P
components was 0.29436%, in the small imaginary off-diagonal coefficient.
None of the signed R/X/G/B or P components changed sign at either sampled
frequency. These observations establish neither a tolerance nor unsampled
frequency behavior, and do not select a production preset.

`observations.md` contains the complete measurement summary and raw G matrices.
`differences.csv` retains every R/X/G/B and P-component difference, including
absolute differences near zero; each configuration's `matrices.csv` retains the
complete complex Z/Y/P. `review.json` records native costs and the component
review. Re-running the bounded script reconstructs its primary observations
from ordinary public computations or compatible completed-run resumes.
