# Static air-boundary experiment

Historical experiment: all executable prototypes, including the independent
annular controls, were archived during the [development cleanup](../fem_development_cleanup.md).
The equations and results below are preserved as research notes; their commands
are not current entry points. The cable comparison used the retired
path-preparation implementation. See the [native voltage contract](../../../../docs/plans/fem-native-voltage-extraction-plan.md).

This isolated experiment tests a first-order static asymptotic boundary on the
upper semicircle of a circular air/earth domain. It does not modify the backend,
the manual runner, or a retained FEM run. All generated meshes, copied solver
inputs and raw results live under the selected output directory.

The lower semicircle retains zero Dirichlet constraints. The air/earth
interface retains the production material transmission equations. No PML is
present. The interior quasi-fw equations, terminal excitation, copper material,
and voltage definition are unchanged. This is an air-boundary diagnostic,
not an exact transparent boundary for the entire layered exterior.

## Boundary equations

For the production convention `exp(j omega t - Gamma z)`, define

```
a = A_z, b = A_t/Gamma, v = phi/Gamma,
h = H_z/Gamma = nu curl(b), kappa = sigma+j omega epsilon,
t = (-n_y,n_x).
```

On a circular air boundary of radius R, the candidate imposes

```
partial_n a + a/R = 0,
partial_n h + h/R = 0.
```

These are exact for the first decaying spatial harmonic of each static scalar
field in a homogeneous exterior. They approximate other harmonics and finite
frequency fields. They are not imposed independently on the gauge-dependent
potentials v and b.

Transverse Ampere gives `kappa*(grad(v)+j omega b) = nu grad(a)-C*h`, where
`C*h=(partial_y h,-partial_x h)`. With `w_t=partial_t v+j omega b_t`, the
second Robin condition is equivalently

```
h = R*(nu*partial_t a-kappa*w_t).
```

Consequently the boundary contributions to the production weak equations are

```
a row:  integral (nu/R)*a*test_a,
b row: -integral h*test_b_t,
v row: -integral h*partial_t(test_v) + integral (nu/R)*a*test_v.
```

The last row follows from the normal current flux and integration by parts
along the boundary. Its endpoint term vanishes because the scalar test function
is zero at the junctions with the earth-side Dirichlet boundary. `air_terms.pro`
substitutes the expression for h into these terms. Boundary entities are
explicitly included in the finite-element supports. The gauge tree starts on
the earth boundary and conductor contours; it no longer pins the air boundary.
The tangential derivative of a uses the trace of the volume curl, with
`curl(a*e_z).n = partial_t a`.

The spatial-harmonic premise is described in David Meeker,
[Improvised Open Boundary Conditions for Magnetic Finite Elements](https://www.femm.info/dmeeker/pdf/TMAG-13-02-0097_R1.pdf),
section II. The coupled quasi-fw weak-boundary reduction above is specific to
this experiment, not a formulation taken from that paper.

## Independent checks and cable comparison

`control.pro` supplies two independent static annulus problems with known
solutions: `(a,v,h)=(0,cos(theta)/r,sin(theta)/r)` and
`(cos(theta)/r,0,-sin(theta)/r)`. Two mesh sizes check both boundary cross terms
as well as the direct Robin term. `run.py` expands the boundary fragment into
the copied control file, since GetDP does not accept Include inside Equation.

The cable comparison uses two bare wires of radius 0.0425 m at (0,1) and (1,1),
rho=0.1 ohm m, and the retained run `run-wjUAdy` as its material/source fixture.
It compares all-zero exterior Dirichlet with the static air condition on the
same mesh, initially at 0.1 Hz and 4641.588833612777 Hz. The radius is the larger
of 5 m and two or four soil skin depths. The prescribed bulk and conductor mesh
sizes are fixed when increasing the radius. New circular meshes are generated;
they are not the original rectangular PML meshes. Receiver-contour quadrature
and voltage-path quadrature are regenerated from each new mesh using the
backend's own `_write_voltage_paths`. The 32-edge conductor polygons are fixed;
the path samples follow the new triangle crossings.

The initial experiment incorrectly reused path samples from the rectangular
source mesh. Its cable outputs under `/tmp/fem-static-boundary` are superseded
by the rerun under `/tmp/fem-static-boundary-corrected`. The independent annulus
controls have no such path extraction and remain valid.

The scalar 2x2 P matrix is read from the native output; `Y=inv(P)` and raw
`G=real(Y)` are reported without clipping or forcing signs. This inexpensive
diagnostic is not a conductor skin-mesh convergence study or a full spectral
acceptance campaign.

## Current independent controls

From the repository root, with Julia/Gmsh and Python/numpy available:

```sh
JULIA_DEPOT_PATH=/tmp/lcm-fem-depot:/home/amartins/.julia \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 10 \
  python3 test/manual/calculations/fem_static_boundary/run.py \
  --output /tmp/fem-static-boundary-controls
```

The script now runs only the independent annular controls. The historical cable
solver and its mesh-dependent preparation mode have been removed. Their
scientific results below remain historical evidence. `--getdp` overrides
the local GetDP artifact path. Each native command and log is retained next to
its output. The scripts use one native solver thread and run cases sequentially.

## Results, 2026-09-26

The independent annulus controls pass. Reducing the prescribed bulk element
size from 0.2 to 0.1 reduces the largest sampled absolute h error from 0.0142
to 0.00882. On the finer mesh, the largest sampled potential error is below
8e-5; the potential that should vanish in the magnetic control remains below
1.7e-6. These controls validate the static boundary assembly and cross-term
signs, not the accuracy of its finite-frequency approximation for layered media.

All eight two-source cable solves completed with regenerated voltage paths:
two frequencies, two radii, and two boundary choices. Total native solve/output
wall time was 120.1 s; individual cases took 6.3–30.9 s with native peak memory
169–650 MB. These costs
exclude mesh preparation and are not a performance comparison at equal error
against the enlarged PML preset.

At four soil skin depths, the raw entries were:

| Frequency [Hz] | Boundary/reference | G11 [S/m] | G12 [S/m] |
| ---: | --- | ---: | ---: |
| 0.1 | Analytical | -1.578139e-24 | -5.953726e-25 |
| 0.1 | Dirichlet | -1.827586e-24 | -6.091587e-25 |
| 0.1 | Static air ABC | -1.420212e-21 | -1.368342e-21 |
| 4641.588834 | Analytical | -2.013703e-14 | -1.792821e-14 |
| 4641.588834 | Dirichlet | -5.244859e-15 | -3.065839e-15 |
| 4641.588834 | Static air ABC | -8.642539e-12 | +9.886701e-12 |

The static closure gives wrong signs for G12 and G21 at the higher frequency
at both radii. It passes only 12/16 entry sign checks. Doubling the radius
reduces its large error but does not resolve those signs. It is rejected as a
replacement for the working preset. Regenerating the path quadrature changes
some reported values but does not change this conclusion or either sign count.
This finding does not reject higher-order
static expansions or a proper layered exterior operator.

The circular Dirichlet controls pass all 16 sampled entry sign checks. Their
amplitudes are not converged: at the higher frequency G11 is about 74% below
the reference magnitude and G12 about 83% below it. This is a useful cheap
control, not evidence of broadband correctness or a validated new preset.

Results: `/tmp/fem-static-boundary-corrected/cable/conductance.csv`,
independent controls: `/tmp/fem-static-boundary/controls/controls.csv`,
plot: `/tmp/fem-static-boundary-corrected/conductance.png`.
`plot.py` draws the unclipped raw conductance samples and existing analytical
reference, with no connecting curves between FEM frequencies. It requires
matplotlib and the reference CSV from the previous manual-case investigation.

No production files, active manual settings, or user-run outputs were changed.
