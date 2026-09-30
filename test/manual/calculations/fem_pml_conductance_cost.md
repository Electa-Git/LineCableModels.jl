# Quasi-fw PML conductance and cost assessment — 2026-09-29

The fixed prescription `(144,144,96)` removes the observed conductance sign
errors in the complete ordinary two-wire study. It also preserves signs in
the selected buried/mixed three-wire and cable-shape fixtures. The manual
two-wire runner now uses this tuple and forwards it to its detached export.
Production equations, voltage extraction, PML stretching and API defaults
were not changed. There is no production acceptance, correction or refinement.

## Prescribed setup

```julia
mesh_controls = (
    domain_skin_depths = 24.0,
    pml_layers = (144, 144, 96),  # side, top, bottom
    pml_grading = (192/191)*log(1536),
    mesh_size_factor = 3.0,
    exterior_mesh_size_factor = 8.0,
    volume_quadrature = 12,
    conductor_geometry_tolerance = 1e-3,
    conductor_skin_depth_elements = 6.0,
    conductor_mesh_growth = sqrt(1.25),
    conductor_skin_depths = 5.0,
    conductor_thickness_elements = 4,
)
result = compute(problem, fem; options=(; mesh_controls...,
    frequency_workers=4, solver_threads=1))
export_data(:onelab, problem, fem;
    file_name="detached/study.pro", mesh_options=mesh_controls)
```

Here `fem` is the existing quasi-fw formulation. These results concern the
near-zero longitudinal propagation constant and the exact physical inputs
below. The tuple is not an accuracy guarantee for arbitrary materials,
frequencies or prescribed propagation constants. In particular, substituting
the API's default physical-domain size does not reproduce this qualification.

## Diagnosis and controls

The error was already present in complete native P before `Y=inv(P)` and
before plotting. Higher-precision inversion and the old PML expression did
not remove it. The change from 96 to 192 intervals was almost entirely an
equal-entry/common-mode real-P offset; see
[the original diagnosis](fem_pml_low_frequency.md).

Eight new serial diagnostic solves were performed, with two source columns
each. All other physical and conductor targets stayed fixed.

| Control | Frequency Hz | G12 S/m | Native seconds | DOFs |
|---|---:|---:|---:|---:|
| Existing 96 mesh, 13-point quadrature | 21.54435 | +5.70733e-21 | 23.58 | 311891 |
| (192,96,96), original grading | 21.54435 | -2.33947e-20 | 39.97 | 532691 |
| (96,192,96), original grading | 21.54435 | -2.31094e-20 | 28.72 | 400787 |
| (96,192,96), original grading | 0.1 | +1.36163e-25 | 29.28 | 400451 |
| (96,96,96), side/top grading 6 | 0.1 | +1.31038e-24 | 22.63 | 311555 |
| (96,96,96), side/top grading 9 | 0.1 | +6.40762e-25 | 22.37 | 311555 |
| (144,144,96), original grading | 0.1 | -1.11209e-25 | 36.53 | 484835 |
| (144,144,96), original grading | 21.54435 | -3.77208e-20 | 36.01 | 485171 |

The two analytical G12 values are -6.79782e-25 and -7.03197e-20 S/m.
Neither higher quadrature nor either tested grading exponent solved the
problem. Side/top count changes followed approximately inverse-square error
reduction. This motivated the intermediate count, rather than an unrestricted
parameter search. The eighth solve was an explicitly recorded extension of
the initial seven-solve diagnostic limit, confirming the successful seventh
solve at the second frequency before running the full selection.

The 0.1 Hz controls preserve all conductor node coordinates. Both probes
preserve physical/PML interface nodes. At 21.54 Hz independent export meshing
moved 17 interior nodes in conductor 1, with maximum nearest-node distance
1.35 mm; conductor 2 matched. Therefore the latter controls alone do not
perfectly isolate exterior discretization. The low-frequency controls and
complete managed sweeps provide the stronger evidence for the prescribed fix.

## Completed scientific selection

One candidate was used without per-case tuning:

| Selection | Frequencies / source columns | G comparisons | Sign disagreements |
|---|---:|---:|---:|
| Aerial two-wire study: 3 radii and 4 soil resistivities | 70 / 140 | 280 | 0 |
| Three buried wires and one aerial plus two buried wires | 6 / 18 | 54 | 0 |
| Screen, tube and sector preservation at 1 MHz | 3 / 9 | 29 | 0 |
| Total managed | 79 / 167 | 363 | 0 |

The two-wire centres are `(0,1)` and `(1,1)` m. The radius sweep is
0.001/0.01/0.085 m in 0.1 ohm m earth. The resistivity sweep is
0.1/1/100/1000 ohm m at radius 0.0425 m. All use ten logarithmically spaced
frequencies from 0.1 Hz to 1 MHz, the original copper material, and relative
soil permittivity/permeability one. Three-wire fixtures retain the established
left-to-right terminal ordering at 0.1 Hz, 21.544346900318832 Hz and 1 MHz.
Every energization/observation pair is retained; no symmetrization is applied.

Wire references are the matching unified analytical formulation with Gamma=0.
Cable references are the existing qualified physical screen, tubular-sheath
and sector designs and saved FEM results, not equivalent round conductors.
Their comparisons are preservation observations, not independent error bounds.

The sign objective is satisfied for this set; magnitude accuracy remains
separate. The worst two-wire G relative error is **84.69%** at the smallest
low-frequency conductance. Maximum G errors in the 1/100/1000 ohm m sweeps are
2.159%, 2.268% and 2.327%. The worst analytical self-R error is **0.4652%**
across all wire selections (0.4548% buried, 0.4311% mixed). Large relative
errors of small mutual components remain visible in the saved tables.

| Selection | Maximum X relative error | Maximum G relative error | Maximum B relative error |
|---|---:|---:|---:|
| All ordinary two-wire cases | 1.620% | 84.69% | 2.936% |
| Three buried | 9.269% | 7.366% | 39.07% |
| Three mixed | 2.433% | 55.99% | 8.080% |

These maxima include mutual entries and are not self-only errors. For cable
preservation, maximum G changes are 0.08158% (screen), 1.606% (tube), and
0.1063% (sector); all retain the saved signs. No new percentage acceptance
threshold has been imposed. Full signed, absolute and relative R/X/G/B errors
are retained without clipping in `qualification/*/components.csv`.

## Cost and detached execution

At the 21.54 Hz diagnostic point, the current-source 192/192/192 baseline
took 67.01 s with 857939 DOFs. The 144/144/96 control took 36.01 s with
485171 DOFs: **46.3% less native wall time and 43.4% fewer DOFs** in this
single controlled comparison. Corner triangles drop from 294912 to 138240.
This is not a repeated whole-study timing claim.

The seven ordinary managed scans total 889.17 s (14.82 minutes), with four
frequency workers and one native thread per worker. The first FEM compute
call includes 26.43 s of compilation; the remaining six recorded zero
compilation time and took 109–141 s each. Package loading, analytical reference
construction and plotting are outside these compute timers. `costs.csv`
separates compilation, GC and accumulated native stages; accumulated worker
wall time must not be confused with elapsed scan time.

One additional exported `.pro` solve used its independently generated mesh
at the failing 0.085 m / 21.54 Hz point. It retained all analytical G signs.
Compared with the separately meshed ordinary managed result, maximum G
difference was 7.08e-26 S/m (about 1.9 ppm of mutual G). That exceeds the
existing stringent same-mesh parity bound and is recorded as such.

To separate execution from remeshing, compare this detached solve with the
already completed managed-worker command on the **same mesh**, verifying its
saved hash. Z and P were identical; Y differed by at most 7.23e-35 S/m in its
real part and 8.28e-25 S/m in its imaginary part. The existing component-wise
parity bound passed unchanged. This comparison reused saved data and added
no solve. Native export and managed execution therefore remain supported;
independent-mesh differences are not presented as bitwise parity.

## Evidence and reproduction

Evidence root: `.linecablemodels/fem/pml-conductance-cost/`.

- `live.log`: single append-only record of decisions, native output and summaries.
- `candidate.toml`, `qualification/*/options.txt`: prescribed controls.
- `qualification/*/{matrices.csv,reference.csv,components.csv,complete.toml}`:
  raw values, errors, exact source-run locations and completion markers.
- `component-summary.csv`, `costs.csv`: component and performance summaries.
- `plots/{aerial-radius,aerial-resistivity,three-buried,three-mixed}/{R,X,G,B}.{svg,png}`:
  public plotting API, `clip=false`.
- `phase2-balanced144-mid/detached/study.pro`: runnable detached problem.
- `detached-parity.csv`: original independent-mesh comparison, including failures.
- `detached-comparisons.csv`, `detached-preservation.toml`: resolved distinction
  between independent meshing sensitivity and same-mesh execution parity.
- `resume.sh`: resumes saved qualification, detached assessment and plotting;
  completed native cases are reused. Run only when no copy is already active.

The focused retained manual regression is
[`run_fem_pml_signs.jl`](run_fem_pml_signs.jl). It covers the three complete
radius sweeps and the two three-wire placements: 36 frequency solves. Its
sign requirement belongs to those fixtures, never to the production engine.
The full ten-frequency definition is retained for the radius cases so that
the mesh planning inputs match the qualified study. This convenience caller
was syntax/load checked; the scientific evidence above comes from the
completed public-compute campaign, not an additional duplicated run.

```bash
tail -n 60 -F .linecablemodels/fem/pml-conductance-cost/live.log
```

The documented preset and focused regression deliver the existing-controls
route in [the execution plan](../../../docs/plans/fem-pml-conductance-cost.md).
They do not introduce a new production mesh algorithm or a universal sign
guarantee.
