# Fixed-node corner diagonals and soil localization

Completed 2026-09-30. Three new serial native solves / six source columns.
The saved aerial reference was reused. Production sources and the user's
existing detached ONELAB bundle were not modified.

## Why the displayed soil mesh was still fine

The existing `.linecablemodels/fem/onelab-two-bare-wires/` bundle was exported
at 14:13, before the both-media localization update. Its first geometry still
contains this obsolete soil condition:

```c
If(MeshWaveEarth < MeshRemoteEarth && MeshRemoteEarth <= MeshScale*133.45547689072069)
```

At its prescribed 0.1 Hz settings, `MeshRemoteEarth` is 603.950545 m and the
extra limit is 133.455477 m. The condition is false, leaving the full-width
soil distance field active. The corresponding air limit is much larger, so
the same bundle localizes air. Remeshing this detached bundle cannot update
its exported field definitions from Julia sources.

Current production `mesh.jl` localizes both media; `export.jl` writes only
`MeshWaveEarth < MeshRemoteEarth` for this condition. Soil and air retain their
own propagation-dependent sizes and decay radii. Identical localization
strategy does not require identical element sizes in both materials.

A fresh public detached export is available at:

`.linecablemodels/fem/pml-corner-mesh/diagonals-and-soil/soil-localization/localized/study.pro`

It contains the manual runner's ten frequencies from 0.1 Hz to 1 MHz, mixed
wires at (0,-1) and (1,+1) m, radius 0.0425 m, rho=0.1 ohm m and Gamma=0.
Only its first frequency was meshed and solved in this focused test. It is
an ordinary detached public export and supports native ONELAB remeshing.
Restart Julia before re-exporting the usual manual runner if its existing
session still holds methods loaded before the source update.

## Isolated soil comparison

Both bundles were freshly exported from current production. The control
disables only the final soil-localization override; air localization remains
enabled in both. Equations, material data, domain, mesh factors, conductor
grading, PML and native voltage-path prescriptions are identical.

| Quantity | Air-only localization | Both-media localization |
| --- | ---: | ---: |
| Soil triangles | 16159 | 6375 |
| Air triangles | 5739 | 5739 |
| PML triangles | 169100 | 169100 |
| Total nodes | 96410 | 91518 |
| DOFs | 381026 | 361458 |
| Native meshing, s | 0.4199 | 0.3463 |
| Native GetDP elapsed, s | 29.14 | 28.21 |
| Assembly, s | 20.2198 | 19.7064 |
| Factorization/solve, s | 6.5277 | 6.2850 |
| Peak RSS, KiB | 1903928 | 1813772 |

Soil triangles decrease by 60.55%; total nodes decrease by 5.07%. PML,
conductor-contour and voltage-path coordinate hashes match exactly. The
maximum raw per-entry R/X/G/B change is **0.832042%**, with **zero sign
changes**, passing the existing 2% preservation limit without clipping or
error floors. The largest change is the aerial conductor's very small G22:

| G22, S/m | Value |
| --- | ---: |
| Air-only localized reference | -4.446181862981705e-25 |
| Both-media localized | -4.483175969756719e-25 |
| Analytical | -9.696067878459929e-25 |

This passes preservation of the retained FEM; it does not establish 2%
agreement with the analytical model. `components.csv` retains analytical
values for every component. The measured native elapsed reduction is 3.19%
in this single matched pair, not a repeated warmed speedup claim. No Julia
compilation is included in native process timings; meshing reported zero
Julia compilation for both calls. PML still dominates the total mesh.

The image `soil-interface-comparison.png` (and SVG) beside the fresh bundle's
parent directory shows the same interface window before and after. Remaining
refinement is concentrated around the cables and their interface footprint;
the narrow boundary-node distribution on the interface is retained.

This focused comparison supplements the completed fifteen-case both-media
qualification in `../fem_local_interface_mesh/both-media-results.md`; it does
not repeat that spectrum or weaken its acceptance policy.

## Fixed-node diagonal experiment

For the aerial 0.1 Hz, Gamma=0 case, all 51912 rectangular corner cells were
split along the opposite diagonal. Node coordinates, node count (91044),
total triangle count (180892) and corner triangle count (103824) are retained.
Physical elements, PML side-strip elements and voltage paths have identical
hashes. All new triangles have positive areas; maximum cell area change is
4.44e-16 relative. Native transfinite coordinates have small interpolation
roundoff, so diagonals are identified relative to each cell's coordinate
spans rather than requiring exactly horizontal/vertical edges.

The worst raw component change is **0.0349019%**, with **zero sign changes**.
G12 changes from -2.3422506890550553e-25 to -2.3430627216740622e-25 S/m.
The native solve took 28.05 s versus the saved 27.80 s reference; this
diagnostic changes no element counts and is not a speed optimization.

Thus simply reversing the diagonal at fixed nodes does not reproduce the
large conductance error in the two anisotropic coarsening candidates. Their
failure is not explained by diagonal direction alone. This does not prove
that all unstructured triangulations are harmless, or that the current
corner count is minimal. The next prescription must control the interpolation
error introduced by moving/removing corner nodes; those earlier candidates
remain rejected and the structured production corners remain unchanged.

Evidence root: `.linecablemodels/fem/pml-corner-mesh/`.

- `air-f0.1-gamma0.0/opposite-diagonals/`: frozen candidate mesh and raw solve.
- `air-f0.1-gamma0.0/components-opposite-diagonals.csv`: all component comparisons.
- `diagonals-and-soil/soil-localization/{baseline,localized}/`: detached bundles,
  meshes, hashes, native solver logs, timing and raw matrices.
- `diagonals-and-soil/soil-localization/components.csv`: soil comparison.
- `diagonals-and-soil/complete.toml`: completion and both acceptance decisions.
- `live.log`: append-only shared execution log.
