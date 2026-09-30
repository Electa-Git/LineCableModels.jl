> Historical investigation record. Equations, observations and limitations are
> retained at their recorded source versions. Referenced prototypes and campaign
> launchers may have been retired; their commands and pending-work statements
> are not current execution instructions. See the [cleanup record](fem_development_cleanup.md)
> and [current PML controls](fem_fixed_pml_controls.md).

# Exterior mesh grading investigation

References below to prepared voltage paths and contour-averaged measurements describe the historical extraction implementation. The recorded mesh and scientific results retain that convention. See the [native voltage contract](../../../docs/plans/fem-native-voltage-extraction-plan.md).

The candidate retains the physical half-width of 24 soil skin depths, the
192-layer Cartesian PML, its coordinate stretch, the terminal excitations,
and the existing voltage references. `mesh_size_factor=3` retains the selected
manual preset's conductor and central-medium size targets. The new computation
control `exterior_mesh_size_factor=8` permits coarser remote elements; its
package default remains 1.

The central bulk cap is preserved out to twice the resolution radius
`max(layout_radius, earth_skin_depth)`. Outside it, the prescribed size grows
with the existing slope of 0.2 to a bounded remote size. The air cap retains
the propagation-length bound. The soil keeps its refinement near the
air/soil interface and conductors. Tangential PML divisions use the remote
targets, with geometric grading towards the interface on the side edges.
Normal PML grading is unchanged. Every remesh regenerates the receiver and
voltage-path quadrature.

This is a mesh change, not a new electromagnetic boundary condition. It does
not replace computed signs, clip values, or change the analytical reference.

## Evidence and limits

Evidence is retained under `/tmp/fem-exterior-grading`. The initial wide-wire
case has two aerial wires, radius 0.085 m, separation and height 1 m, and
soil resistivity 0.1 ohm m. At 0.1 Hz the 192-layer graded mesh reduces the
unknown count from 2,050,583 to 856,819. Fresh one-thread native solve/output
time is 203.185 s and reported peak memory 4048.11 MB; the retained ungraded
comparison took 450.805 s and 10540.3 MB on the same host. These measurements
are native costs, excluding Julia compilation and mesh preparation, and are
subject to concurrent host load.

The first two-frequency public compute call took 258.96 s with 47.40 s of
reported compilation. The subsequent warmed 192-layer doubled-thickness call
took 192.73 s with zero reported compilation. The latter uses a different
thickness and is not a same-input warm/cold speed ratio. Per-call compilation,
allocation and wall-time measurements are in each case's `execution.csv`.

Both 192-layer thicknesses preserve all eight conductance signs at 0.1 Hz and
4641.588833612777 Hz. At 0.1 Hz, G12 is -3.775194e-25 S/m with the original
thickness and -2.248051e-25 S/m with twice that thickness, against
-6.797818e-25 analytically. This is sign stability at these samples, not
relative convergence of the smallest conductances.

The 96-layer variant is rejected: both thicknesses give positive mutual G
at 0.1 Hz. Its first-thickness G12 is +8.165909e-25 S/m. No other parameter
is adjusted to hide this failure. The 100 and 1000 ohm m checks at 1 MHz
also retain all four reference signs with the 192-layer graded candidate.

The five standard-radius rho=0.1 frequency checks also pass all 20 signs.
Together with the two wide-wire and two high-resistivity checks, the
original-thickness candidate passes 36/36 entry sign checks. The manual runner
now selects `exterior_mesh_size_factor=8.0` alongside its existing 24/192/3
controls. The package default remains 1.

The actual manual runner's smaller radii, 0.001 and 0.01 m, also pass all
24 entry sign checks at 0.1 Hz, 4641.588833612777 Hz and 1 MHz with rho=0.1.
This brings the targeted original-thickness checks to 60/60. Their evidence
is under `g8-r0.001` and `g8-r0.01`.

The first complete 99-frequency sweep is quasi-fw, two aerial wires,
rho=1 ohm m, standard radius. All 396 entry signs agree with the analytical
reference from 0.1 Hz through 1 MHz. The maximum relative G difference is
1.694%, and the maximum relative B difference is 1.757%. These accuracy
figures apply to this case, not to the small rho=0.1 conductances above.

The full 18-case campaign continues under `/tmp/fem-exterior-grading/dense`,
with six single-thread native workers. It covers both formulations, three
placements, three resistivities and 99 frequencies (1782 case/frequency
points). Completing its first two cases does not establish all-layout acceptance.
Completed CSVs and plots are copied to
`.linecablemodels/fem/exterior-grading-evidence`; `status.txt` records completed
case counts, and `sign_mismatches.csv` records reference-sign differences.
Quasi-tem is compared separately against the independent scalar reference.
The buried rho=1 mutual conductances have a genuine reference zero between
794328.2347242822 and 857695.8985908945 Hz in both references. Differences in
that zero's location must be assessed as a convergence issue; they must not
be removed by imposing a common sign on the spectrum.

## Reproduction

From this checkout, the existing manual validation driver accepts the new
control. Restart Julia before running the updated interactive manual script:
the backend mesh-plan type has gained fields. A targeted offline wide-wire
check is:

```sh
JULIA_DEPOT_PATH=/tmp/lcm-fem-depot:/home/amartins/.julia \
JULIA_LOAD_PATH=@:@stdlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
PML_DOMAIN_SKIN_DEPTHS=24 PML_LAYERS=192 PML_MESH_SIZE_FACTOR=3 \
PML_EXTERIOR_MESH_SIZE_FACTOR=8 PML_THICKNESS_FACTOR=1 \
PML_RADIUS=0.085 PML_RHOS=0.1 PML_LAYOUTS=air_air PML_PHYSICS=quasi_fw \
PML_FREQUENCIES=0.1,4641.588833612777 PML_WORKERS=2 \
PML_OUTPUT=/tmp/fem-exterior-grading/repeated-wide \
julia --startup-file=no --compiled-modules=existing --project=test \
  test/manual/calculations/run_fem_pml_validation.jl
```

Use a separate output directory for each control change. `PML_THICKNESS_FACTOR=2`
and `PML_LAYERS=96` select the two independent comparisons. Omitting the
frequency/layout/physics/resistivity selectors uses the existing 99-frequency,
three-layout, two-formulation, 1/100/1000 ohm m campaign. Its runtime is much
larger than this targeted check.

```sh
python3 test/manual/calculations/summarize_fem_exterior_grading.py \
  /tmp/fem-exterior-grading
MPLCONFIGDIR=/tmp/fem-pml-matplotlib /usr/bin/python3 \
  test/manual/calculations/plot_fem_pml_validation.py \
  /tmp/fem-exterior-grading/g8-wide \
  /tmp/fem-exterior-grading/g8-n96-t1-wide \
  /tmp/fem-exterior-grading/g8-n96-t2-wide \
  /tmp/fem-exterior-grading/g8-n192-t2-wide
```

The summary writes ordinary `conductance.csv` and `costs.csv`. Plots read the
matrix CSVs without JSON or any zero cutoff. Quasi-tem differences from the
unified analytical formula require separate comparison with its scalar-PDE
reference; a sign count alone does not establish model equivalence.

## Regression checks

Mesh, option and resume tests passed 209 assertions before adding export
round-trip coverage. The subsequent grading/export checks passed 58 assertions,
including exact node-count reproduction after reopening the graded `.geo`.
Detached Python/native numerical parity passed 139 assertions using
`LINECABLEMODELS_ONELAB_PYTHON=/tmp/lcm-onelab-export-env/bin/python`.
The default shell Python lacks Gmsh; its initial parity attempt failed before
running the detached solver.

The export check exposed two required corrections: preserve string-valued
mesh-field expressions and export the actual configured global mesh-size cap.
Neither changes the GetDP weak forms.

A large-core preservation check also exposed an interior coarsening error in
the first grading implementation. When a local cable size exceeded the old
bulk cap, lifting that global cap reduced a 3 m-radius core from 50 nodes to
29. The old bulk cap is now retained on all cable materials. The dedicated
regression passes, and the final mesh/resume selection passes 113 assertions.
Mesh-cache generation is incremented to prevent reuse of affected meshes.

The material-cap correction produces byte-identical `.msh` files for the
standard and wide low-frequency wires, the standard intermediate-frequency
case, and the 1000 ohm m/1 MHz case. Logs and SHA-256 comparisons are in
`/tmp/fem-exterior-grading/cap-equivalence.log`. Thus the running thin-wire
campaign and the targeted numerical results remain applicable; their native
input snapshots retain the precise sources actually used.

The manual runner's four-thread setting also passes the wide-wire 0.1 Hz
check (`threads4-wide`). The dense campaign uses one thread per solve to
measure costs consistently. A separate map campaign selects both formulations,
all three placements, rho=100 ohm m, and 0.1 Hz/1 MHz under the same grading
controls. Its outputs are under `/tmp/fem-exterior-grading/maps`.

The final mesh selection passes 37 assertions, including the near-equal
endpoint-size regression. The geometric-spacing logarithm uses `log1p` near
equal endpoints. Across 4752 planned edge distributions (all three layouts,
four resistivities, two radii and 99 frequencies), its final values match
the initially tested distributions exactly. This avoids an unnecessary
last-bit change to the generated meshes; see `grading-roundoff-final.log`.
That check also verifies that the material cap is inactive for every wire
fixture in this campaign: the cable-distance field remains below the bulk
cap everywhere inside each circular core.

The first completed map pair is quasi-fw, aerial, rho=100 ohm m. Both endpoints
match all four analytical G signs. Native material tags show physical air and
earth and the corresponding PML partitions in the expected locations. Magnitude
and phase maps are under `maps-preview/quasi-fw-air-air-rho100-0.1` and
`maps-preview/quasi-fw-air-air-rho100-1MHz`. This checks material ownership and
field behaviour at these samples; it does not by itself measure reflection.

For the standard-radius quasi-tem campaign, the independent scalar reference
can be included in the CSV summary:

```sh
python3 test/manual/calculations/summarize_fem_exterior_grading.py \
  /tmp/fem-exterior-grading \
  /tmp/fem-pml-execution/scalar-reference-full/matrices.csv
```

The scalar-reference fixture has 0.0425 m wire radii. It must not be used
for a different radius without recomputing that independent reference.

## Completed checks and remaining model discrepancy

The aerial rho=1 quasi-tem sweep also completes all 99 frequencies. Its
396 conductance signs agree with the independent scalar reference; maximum
relative differences are 1.410% for G and 1.751% for B. This validates the
discretized scalar equation for this case, not its equivalence to quasi-fw.

All six rho=100 map cases complete: both formulations, all three placements,
and 0.1 Hz/1 MHz. Quasi-fw matches all 24 unified-reference signs. Quasi-tem
matches all 24 independent scalar-reference signs, but one of those signs
differs from the unified reference. The mismatch is G[1,2] for the mixed pair
at 1 MHz: conductor 1 is in air and conductor 2 is in earth.

| Calculation | G[1,2], microSiemens/m |
|---|---:|
| Unified analytical reference | -11.3545525 |
| Quasi-fw FEM | -11.5872940 |
| Existing quasi-tem FEM | +1.2675728 |
| Independent finite-electrode scalar PDE | +1.2684146 |

The scalar FEM differs from the independent scalar calculation by 0.0664%
here. Its opposite sign relative to the unified formula is therefore not
resolved by this mesh/PML correction. Reconciling it requires examining the
scalar formulation and its voltage definition. The coupled quasi-fw result
has the expected negative sign. No weak form or voltage definition is changed
by this exterior-grading implementation.

Independent refinement confirms that this positive scalar result is stable:

| Electrode panels | Spectral quadrature | Scalar G[1,2], microSiemens/m |
|---:|---:|---:|
| 64 | 256 | +1.268414611447263 |
| 64 | 512 | +1.268414611447265 |
| 128 | 256 | +1.268414856385922 |
| 128 | 512 | +1.268414856385935 |

Doubling both resolutions changes G[1,2] by 0.00001931%. This independent
calculation has neither a finite FEM boundary nor a PML. Reproduce it with:

```sh
OPENBLAS_NUM_THREADS=1 PML_LAYOUTS=air_earth PML_RHOS=100 \
PML_FREQUENCIES=1000000 PML_BEM_PANELS=64,128 PML_BEM_QUADRATURE=256,512 \
PML_BEM_OUTPUT=/tmp/fem-exterior-scalar-check \
python3 test/manual/calculations/run_fem_scalar_reference.py
```

The independent calculation also isolates the sign change algebraically.
Let `C` map absolute electrode potentials to currents and let `M` map those
potentials to the requested voltages. Only the aerial row subtracts the
potential on the earth surface, so the mixed case has `M[2,:] = [0,1]` and
`Y = C * inv(M)`. Consequently,

```text
Y[1,2] = C[1,2] - C[1,1] * M[1,2] / M[1,1].
```

At the refined resolution, the real parts of the two terms are -12.1609402
and -13.4293550 microSiemens/m. Their difference is +1.2684149. Thus this
particular positive entry arises when the existing scalar solution is
expressed in its surface-referenced voltages. Dropping the reference term
would change the requested observable; it is not an established correction.
Agreement between two scalar solvers does not establish equivalence to the
unified electromagnetic model.

The mixed quasi-tem conductance plot under
`.linecablemodels/fem/exterior-grading-evidence/maps/plots/`
includes the unified reference and the independent scalar curve. The
`scalar-reference` legend entry identifies the latter. The map campaign has
only two frequency endpoints, so its connecting lines do not locate a zero.
`sign_mismatches.csv` retains disagreement against either reference rather
than suppressing the unified-reference mismatch.

The larger 18-case campaign remains running after these checks. Its remaining
cases are not claimed as passed. Completed plots and raw comparisons update
automatically in the evidence directory named above.
