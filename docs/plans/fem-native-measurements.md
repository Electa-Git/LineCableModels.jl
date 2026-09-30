# Native voltage measurements and matrix-operation evidence

Validated on 2026-09-28. This records native measurement semantics and the
independent extraction study. Production exporter validation is recorded in
[fem-onelab-export-validation.md](fem-onelab-export-validation.md).

## Implemented measurement replacement

The shared GetDP formulations own field evaluation and integration. Both
execution routes use the same native measurements without external path
preparation or generated quadrature files.

The geometry defines one physical `LCM/voltage_path/NNNN` curve group and one
`LCM/voltage_reference/NNNN` point group for each terminal. A path is vertical,
from the air/earth interface for an overhead terminal, or the outer bottom PML
boundary for a buried terminal, to the lowest existing terminal CAD vertex
(lowest x breaks ties). The buried path has a break at the physical/PML
interface. Native transfinite grading uses the existing volume-mesher grading
rule and PML progression. It avoids imposing fine terminal spacing over the
entire far-field path. These geometric entities remain readable in exports.

Quasi-TEM uses the grouped terminal scalar potential and the native point
reference. Quasi-full-wave additionally stores the pulled-back vector potential
with `StoreInField`, then applies `ComplexVectorField`, `Integral`, a surface
Jacobian for lines in 2D, and four-point Gauss quadrature. The stored field
excludes metal interiors, whose transverse-field contribution is zero. This
requires a GetDP build with Gmsh support, including the distributed binaries
used by the backend and the installed ONELAB bundle. No field-map output file
is required by the measurement operation.

One explicit path replaces the former bare-overhead contour average. Field
equations, current normalization, scalar boundary conditions, gauge constraints,
terminal order and the pulled-back PML circulation convention are unchanged.
This is an intentional measurement-convention change, not a claim of bitwise
equivalence to contour averaging.

Existing meshes without the new physical groups require regeneration. The mesh
fingerprint version was advanced; standalone external-mesh validation checks
the measurement groups. Existing exported bundles, including
`.linecablemodels/fem/onelab-two-bare-wires`, are not rewritten by this change;
re-export them to obtain the replacement. Historical research scripts using
the deleted private clipping helpers are not maintained compatibility clients.

## Matrix inversion: verified native operations

The official [GetDP manual](https://getdp.info/doc/texinfo/getdp.html#Extended-math-functions)
documents `Inv[expression]`; its `Tensor` representation is 3 by 3. A complex
2 by 2 matrix embedded as `diag(P,1)` was inverted inside a `PostOperation` on
the installed GetDP 3.6.0-git build. The maximum entry of `P*Y-I` was
3.1031676915590914e-17. The experiment is retained in
`/tmp/lcm-native-matrix-audit/model.pro`.

For arbitrary terminal count, a native algebraic `Formulation` uses one global
unknown per row, `GlobalTerm` coefficients from P, and `Generate`/`Solve` to
solve `P*y_j=e_j`. `PostOperation` exports each solution column. This delegates
factorization to GetDP's solver; it does not implement elimination in the input
language. The maintained fixture `test/fixtures/data/fem/native_getdp/inverse.pro`
tests a complex nonsymmetric 4 by 4 matrix. On the installed build its maximum
`P*Y-I` residual was 4.440892098500626e-16. The fixture also passes with the
backend's GetDP 3.5 binary.

The production `getdp/line-parameters.pro` uses this native algebraic route
for arbitrary matrix sizes, bundle/Kron reduction and admittance inversion.
The maintained 1-, 2- and 4-terminal fixtures exercise complex coefficients,
all reduction flags, repeated coefficients in one process, and singular failure.

## Validation and numerical scope

The current focused native test command and results are recorded in
[fem-onelab-export-validation.md](fem-onelab-export-validation.md). Native
manufactured extraction also tests excitation-dependent stored fields and
maps-disabled execution, at the existing numerical tolerances.

The analytic edge-field fixture uses `(1+2j)*(1-y,x,0)`. A nonconforming open
line gives `0.24+0.48j`; a closed boundary gives `2+4j`. Both match to rounding
error. Field interpolation at `(0.31,0.44)` is independently checked. This
validates complex arithmetic and orientation without an external path builder.

Nonconforming line quadrature introduces an integration discretization error
when a path crosses jumps in the piecewise finite-element field. To measure it
separately from contour averaging, an offline comparison used the former
triangle-partition quadrature with **one endpoint per terminal**, on the same
solved mixed-wire toy mesh at 10 kHz and eight PML layers. The maximum relative
complex P-entry difference was:

| Line subdivisions relative to initial grading | Maximum relative difference |
|---:|---:|
| 1 | 1.122896% |
| 2 | 0.072016% |
| 4 | 0.291954% |
| 8 | 0.079424% |
| 16 | 0.070371% |
| 32 | 0.018727% |

The refinement is not monotone, as expected for quadrature across element
boundaries. These are extraction errors on a deliberately coarse FEM solution,
not estimates of physical-model accuracy. The native method is not exactly the
old triangle-partition rule; measurement-line refinement remains part of a
convergence study. No production retries or acceptance thresholds were added,
and no existing numerical test tolerance was relaxed.

Local audit artifacts are `/tmp/lcm-native-replacement-tests.log`,
`/tmp/lcm-native-integration-compare.log`,
`/tmp/lcm-native-refine-lines.log`, and `/tmp/lcm-native-refine-lines-fine.log`.
Task-only pre-edit copies and a comparison patch are under
`/tmp/lcm-native-replacement-backup` and `/tmp/lcm-native-replacement.patch`.
