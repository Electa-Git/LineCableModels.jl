# Low-frequency aerial conductance after PML coarsening

The subsequent [cost and conductance assessment](fem_pml_conductance_cost.md)
completed the prescribed-controls investigation. It qualifies `(144,144,96)`
for the recorded sign fixtures and documents the remaining magnitude errors.
The observations below describe the earlier failure and its reproduction.

## Reproduction

The user's two-wire runner uses quasi-fw, copper wires at `(0,1)` and `(1,1)` m,
0.1 ohm m earth, relative earth permittivity/permeability one, and radii
0.001, 0.01 and 0.085 m. Its ten frequencies are
`10.0 .^ range(-1,6; length=10)`. The relevant prescribed controls are
24 soil skin depths, `pml_layers=(96,96,96)`,
`pml_grading=(192/191)*log(1536)`, `mesh_size_factor=3`, and
`exterior_mesh_size_factor=8`; conductor controls are unchanged.

The saved runs confirm those settings, so this is not a scalar/tuple option
normalization failure. The 0.001/0.01 m cases have positive G12 at their first
three samples. The 0.085 m case has positive G12 at its first four samples,
through 21.544346900318832 Hz. Its negative analytical reference is obtained
from the same serialized physical problem with `formula(:unified; options=(Γ=0.,))`.

| Radius 0.085 m, frequency | G12 with 96/96/96, S/m | Earlier 192/192/192, S/m | Analytical G12, S/m |
|---:|---:|---:|---:|
| 0.1 Hz | +6.31264e-25 | -3.68364e-25 | -6.79782e-25 |
| 3.59381 Hz | +6.17340e-22 | -8.45771e-22 | -1.29496e-21 |
| 21.54435 Hz | +5.70728e-21 | -5.27577e-20 | -7.03197e-20 |

The earlier 192-layer radius sweeps have negative G12 at all ten samples.
They establish the observed sign behavior for these fixtures, not a universal
accuracy bound. Even their conductance magnitudes remain imperfect.

## What the retained data establish

The sign problem is present in raw native P and in a direct `Y=inv(P)` evaluation,
before observation clipping or Makie. Inverting those same double-precision
coefficients with 256-bit arithmetic reproduces the wrong signs. The 2-by-2 P
condition number for the 0.085 m case is about 1.69. Increasing the precision of
this final inversion cannot restore information missing from its input.

Increasing only the bottom count from 48 to 96 leaves the aerial defect almost
unchanged. Both saved sweeps retain the same wrong-sign samples. At 0.1 Hz the
0.085 m G12 values differ by approximately 1.75e-30 S/m.

The issue concerns a tiny real component within a predominantly imaginary
admittance. At 21.54435 Hz, B12 is about -6.68e-10 S/m while the conductance of
interest is around 1e-20 S/m. Agreement in the complex norm or the susceptance
does not establish the conductance sign. The evidence points to mesh-dependent
error in the field/voltage coefficients overwhelming this small component;
it does not justify clipping, forcing signs, or substituting the analytical G.

The older 192-layer runs predate the native PML expression optimization. Two
fresh serial controls reuse their existing meshes at 21.54435 Hz and change
only the copied `pml.pro`, isolating that optimization from mesh resolution.
Original runs and production sources are preserved. Their results and logs are
recorded under the evidence root below.

Both controls completed (four source columns, no remeshing):

| Mesh and expression control | G12 at 21.54435 Hz, S/m | Native wall seconds |
|---|---:|---:|
| 96-layer mesh, old PML expressions | +5.707279332948888e-21 | 35.33 |
| 192-layer mesh, current PML expressions | -5.275796312240106e-20 | 67.01 |

The old expression gives the same positive G12 on the 96-layer mesh. Updating
the expression on the 192-layer mesh retains its negative value, differing from
the old result by about 2.52e-25 S/m. Thus the expression optimization does not
explain the sign failure. These are native execution times for diagnostic
controls, not compilation/warmed whole-engine performance measurements.

## Interpretation and next control

The earlier fixed-PML observation covered 0.1 Hz and 1 MHz in one mixed-layout,
100 ohm m case with radius 0.0425 m. Its absence of sign changes cannot establish
accuracy for an aerial, 0.1 ohm m radius sweep. The mesh-control feature is
usable; the 96-layer setting has failed this additional accuracy case.

For reproducing the earlier sign behavior, retain the 192-layer prescription.
Reducing cost requires testing side/top resolution or grading against the raw
real components at the failing low-frequency samples. Increasing the bottom
count alone has already failed to address them. This investigation has not
qualified a cheaper replacement, nor changed the user's runner or API defaults.

The continuous PML matching property does not supply a discrete error bound.
[COMSOL's PML implementation reference](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_definitions.21.137.html)
likewise distinguishes absorption from mesh resolution and discusses grading
for mixed propagating and evanescent components. That general explanation is
consistent with this sensitivity; the local controls are the case-specific
evidence.

## Evidence

`.linecablemodels/fem/pml-low-frequency-diagnosis-20260929/` contains the Julia
readers, the 85-run catalog, raw coefficient/admittance comparisons, matching
analytical references, and the two source-expression controls. No Python or
new production validation machinery is involved.

Relevant saved runs:

- `(96,96,96)`: `run-KgMjbL`, `run-YWmojw`, `run-WFpqS5`, in ascending radius.
- `(96,96,48)`: `run-WkyEy1`, `run-ZeU5x1`, `run-Rk7YKH`.
- Earlier `(192,192,192)`: `run-zXvH76`, `run-zRwwXx`, `run-vn7GJ1`.
