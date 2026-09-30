> Historical investigation record. Equations, observations and limitations are
> retained at their recorded source versions. Referenced prototypes and campaign
> launchers may have been retired; their commands and pending-work statements
> are not current execution instructions. See the [cleanup record](fem_development_cleanup.md)
> and [current PML controls](fem_fixed_pml_controls.md).

# Conductor mesh qualification: execution evidence

Current assessment, 2026-09-29: see
[delivery-assessment.md](fem_conductor_mesh/delivery-assessment.md) for the
completed matrices, plots, accuracy/cost table, diagnosed screen dielectric
resolution defect and bounded coupled confirmation. It supersedes the pending
status statements in the historical chronology below.

All 36 cable cases / 108 source columns completed. Primitive Z/P agree exactly
in all six same-mesh native comparisons. Cable self-R refinement differences
are at most 0.062015% (screen), 0.00078242% (tube) and 0.086865% (sector).
All 18 cable figures exist as SVG and PNG. The plot caller now retains complete
R/X and G/B pairs through the supported observation API; numerical results and
captured numerical source are unchanged.

Screen and sector shunt entries remain unresolved by the original conductor-only
refinement. A local GetDP electrostatic control reproduces the screen B23 and
shows a large field-discretization error on fixed polygon boundaries. Its
four coupled 1 MHz controls now confirm the local result; the correction is
implemented in the existing mesh/export owners and passes 125 focused checks.
All eight final public 1 MHz cases completed and were assessed. Delivered screen
Y agrees with the qualified prototype within 0.003100%; sector Y changes by
0.067236% under refinement, with loop R changing by 0.069526%. The conductor
mesh feature is delivered through both compute and detached export. The
bare-wire conductance sign discrepancies remain open; this closes the mesh
delivery scope, not all FEM accuracy questions.

Earlier screen-path and sector-fill recoveries remain recorded under
`cable-path-diagnosis/` and `sector-fill-diagnosis/` in the qualification root.
The former corrected internal voltage-path sizing while preserving conformity
checks. The latter removed redundant same-PVC seams, preserving all 33 metal
shapes and physical coefficients. Their completed results are retained.

## Historical qualification record

The voltage-path preparation timings and contour-averaged extraction discussed
below refer to the retired extraction implementation. This does not change the
recorded geometry and conductor-mesh results. See the live
[FEM measurement description](../../../docs/src/fem.md) and the
[realigned integration contract](../../../docs/plans/fem-conductor-mesh.md#first-delivery-and-current-integration-contract).

Scope update, 2026-09-28: the first end-to-end delivery is the existing
three-bare-conductor fixture, followed by the two-wire manual caller. Broader
shape work is deferred from that MVP. The current shared GetDP formulations use
conforming field-edge voltage paths with direct `BF_Edge` trace integration;
neither historical prepared paths nor the intermediate nonconforming
stored-field integration are the execution contract. No new full-line solve or
production conductor-mesh change was made during this realignment. Earlier
results below retain their captured source and measurement conventions.

2026-09-27. The [implementation plan](../../../docs/plans/fem-conductor-mesh.md)
is being executed. Q0 source capture is complete; Q1 is in progress, round
Q2-A reference/refinement controls have passed. The captured-equation full-line
spectrum and follow-ups are complete; their preservation on current extraction
code remains to be checked. Corrected screen refinement/reference comparisons
are complete. Shape/contact, native-parity and whole-engine cost work remains.
No production feature implementation has started.
Only this main agent executes the campaign.

## Source and environment

The pre-feature dirty source, tracked diff, source hashes, git status, revision
and manifests are retained in
`.linecablemodels/fem/conductor-mesh-qualification/q0-20260927/` (302 captured
source/input files). Existing unrelated changes were preserved. The experiments
reuse `/tmp/lcm-onelab-export-env/bin/python`, with Gmsh 4.15.2 and NumPy 2.4.6.

The [manual experiment scripts](fem_conductor_mesh/README.md) write normal logs,
CSV measurements, native geometry and meshes. No engine validation or retry
machinery has been added. The following are mesh observations, not electromagnetic
accuracy or runtime qualifications.

## Q1 observations

| Control | Measured observation |
|---|---|
| Two disconnected circles, different normal sizes | Both active fields work simultaneously. Targets 0.01 and 0.03 produced median height/target 0.999995 in both regions. The surrounding-medium first-layer heights stayed unchanged. |
| Annulus, wall/delta = 0.2, 2, 20 | Both faces receive their prescribed sizes; median height/target 0.999995. All interface edges have one triangle on each material side. |
| 96-segment circles/annuli | Relative conductor area error 0.0713794%. All inspected triangle Jacobians are positive. |
| Two touching rectangular materials | Independent first-layer sizes work on both sides of the shared curve, with median height/target 1 and no interface-edge errors. |
| 24 disconnected screen wires, alternating sizes | All 24 fields act on the intended conductor sides. API mesh: 16492 nodes total; median height/target 0.999995 in every wire. No interface-edge errors. |
| 270-degree sector, fan at the reentrant apex | API mesh: 852 nodes; median height/target 0.999995; no interface-edge errors and positive Jacobians. |
| 120-degree sector with direct boundary-layer fields | Failed with `Edge not recovered`. Forcing fans at all corners, using no fans, refining radial curves, and using an apex fan did not resolve this fixture. These failed recipes remain recorded. |
| 120-degree sector, explicit native transfinite strips plus a small core | Meshed with 2567 nodes and 4438 conductor triangles. Relative area error 0.0167732%; maximum first height/target 0.904. No interface-edge errors; positive Jacobians. All partitions retain one conductor physical group, and internal cuts are excluded from its measured contour. |
| Current mixed-pair exterior/PML | All 399360 PML triangles and 200720 unique PML nodes preserved (coordinates within absolute 1e-10 m); physical memberships identical. Conductor layers added 1141 nodes, or 0.528% of the geometry-only mesh. Probe deltas, not actual 0.1 Hz conductor skin depths. |
| Thin annulus, wall/radius=0.001, wall/delta=0.2, 2, 20 | Positive Jacobians and conforming interfaces; area error 0.0713794%. Initial meshes had 677, 869 and 1637 nodes. |
| Quarter-annular bent strips, wall/radius=0.1 and 0.001 | All six wall/delta controls meshed, without fans. Positive Jacobians, conforming interfaces, area error 0.0713794%. |
| Physical 3 mm copper wire, 50 Hz / 10 kHz / 1 MHz | 432 / 651 / 980 total nodes in the initial local fixtures. No boundary layer at 50 Hz; conductor finite section resolved directly. At 1 MHz, prescribed first size 22.03 micrometres was honored. These counts exclude the large production exterior. |

The 120-degree construction is an explicit alternative qualification experiment,
not a fallback added to production. Its fixed counts/ratios are fixture inputs;
general sizing from user controls and electrical convergence remain unqualified.
These failures show why direct boundary layers cannot yet be selected as the
universal construction for every supported shape.

## Native-file activation

Gmsh's `.geo_unrolled` writer retained field definitions but omitted their
boundary-layer activation. Reopening the raw file reproduced the unrefined
interiors: for the pair, median height/target became 5.87 and 1.96. Explicit
`BoundaryLayer Field = <id>;` statements restored both prescribed normal sizes.
The same omission and correction were observed on the annulus, touching-material
pair, 270-degree sector and screen.

Native activated and API meshes preserve the measured sizes and interfaces, but
some unstructured node counts differ slightly (pair: 1891 versus 1887; screen:
16492 versus 16490). This is not a claim of identical mesh numbering or matrices.
The existing export owner already emits selected constraints at full precision;
the eventual feature must similarly preserve layer activation and its parameters.

## Layer depth and thin-wall coverage

The initial normal-column probe required exact radial alignment and missed
Gmsh's slightly deflected layer nodes on circles. This was an instrumentation
failure, not a failed mesh. `q1-normal-probe` corrects it using a band of 0.001
radius about the 45-degree column, much narrower than its angular neighbours.
Unstructured controls need not have an aligned column; absence is reported.

Measured growth is 1.25, but the initial `h=delta/3, Thickness=5*delta` recipe
only emitted six complete layers, ending at about 3.755 delta. Gmsh did not add
the next layer because it exceeded the requested extent. The thin-wall weak-skin
control similarly needed an explicit complete-layer prescription.

The separate `--fit-extent` experiment computes the layer count before meshing:
`n=ceil(log(1+(g-1)*extent/h)/log(g))`, then decreases the first size to
`extent*(g-1)/(g^n-1)`. This is fixed input-based construction, not assessment and
retry. A small floating-point inclusion margin on `Thickness` admits the last
layer. At strong skin effect it gives seven layers, measured growth 1.25 and
normal-column depth 5.00268 delta (the radial/normal direction difference causes
the small excess). The weak-skin thin annulus has five elements across its wall.
API and explicitly activated native files agree on these measured properties.
The 3 mm, 1 MHz local mesh then has 1065 nodes. General corner sizing and the
transition from graded layers to unstructured bulk remain to be qualified.

## First isolated electrical controls

`qualify_impedance.py` and `internal_impedance.pro` solve an independent,
conductor-only diffusion fixture with GetDP. This is a qualification problem,
never an alternative production solver. It uses the same first-order triangle
family and experimental Gmsh grading. The outer axial electric field is 1 V/m;
a tube's inner boundary has zero normal derivative. The exact cylindrical Bessel
solution supplies the reference; DC uses rho/area. Full line equations, earth
return, voltage paths and PML are absent and need their separate Q2-B coverage.

Eight pilot solves finished successfully in `q2-pilot`. For the 3 mm solid wire:

| Frequency / mesh | Internal R error | Complex internal Z error |
|---|---:|---:|
| DC, 96 segments | 0.0714304% | 0.0714304% |
| DC, 192 segments | 0.0178509% | 0.0178509% |
| 1 MHz, candidate (96 segments, delta/3, growth 1.25) | 0.65444% | **1.04310%** |
| 1 MHz, combined refinement (192 segments, delta/6, growth sqrt(1.25)) | 0.16312% | 0.25947% |

For the 3 mm outside-diameter tube with inner/outer radius 0.8, the corresponding
1 MHz complex errors are 0.71320% and 0.21772%. Integrated Joule loss agrees with
Re(Z)*abs(I)^2 to relative 1.8e-15 or better in the pilot. These are solved values,
not polygon-area corrections. The initial candidate narrowly misses the locked
1% complex-impedance target on the solid wire and is **not qualified**. Refinement
choices and convergence still need the spectrum results.

Before launching the spectrum batch, SHA-256 comparison found zero changes in
all 298 production `src/` and `ext/` files captured at Q0.

## Completed round-conductor spectrum

`q2-round-spectrum` completed all 252 solves: seven sections (five solid radii
and two tubular sections), DC plus eight AC frequencies, and four mesh levels.
The original thresholds remain unchanged.

| Mesh level | Cases meeting reference targets | Worst R error | Worst complex Z error | Largest local mesh |
|---|---:|---:|---:|---:|
| 96 segments, delta/3, growth 1.25 | 43/63 | 0.67744% | 1.10080% | 1779 nodes |
| 192 segments, delta/3, growth 1.25 | 43/63 | 0.66437% | 1.10574% | 3430 nodes |
| 96 segments, delta/6, growth sqrt(1.25) | **63/63** | **0.19712%** | **0.26715%** | **3123 nodes** |
| 192 segments, delta/6, growth sqrt(1.25) | 63/63 | 0.17893% | 0.27076% | 6118 nodes |

Doubling tangential geometry resolution did not cure the failed strong-skin
controls. Refining the normal grading did. This supports delta/6 and growth
sqrt(1.25) as the next experimental candidate; it does not yet freeze production
defaults. The largest resistance change from delta/3 to delta/6 was 0.4872%,
so that pair of normal levels does not satisfy the final 0.25% convergence gate.
A further independent normal refinement is reported below. The delta/6 meshes
changed by at most 0.0546% in resistance when only geometry was doubled.

All 252 meshes had conforming interfaces and positive triangle Jacobians.
Integrated loss/terminal power mismatch stayed below 2.1e-14. DC errors were
0.0714304% at 96 segments and 0.0178509% at 192. The largest GetDP peak RSS was
28892 KiB for the 96-segment delta/6 level; this is the isolated diffusion problem,
not a production cost comparison. Full exterior and coupled-system costs remain
unmeasured for this candidate.

`q2-round-summary/levels.csv`, `refinement.csv`, and `reference-errors.png/.svg`
are generated directly from the retained CSV by `summarize_impedance.py`.
The error plot was visually inspected. No solver was rerun to generate it.

`q2-round-normal-fine` subsequently completed 63 controls with 96 segments,
delta/12 and growth `1.25^0.25`. All 63 meet the same independent reference
targets. Relative to delta/6, the largest final resistance change is
**0.128335%**, below the locked 0.25% refinement target. The largest fine-mesh
reference errors are 0.0764673% in R and 0.103129% in complex Z. This completes
the round/tubular isolated reference and normal-refinement gate; it does not
establish full-system accuracy, proximity loss or production runtime.

## Owned sector geometry and native partitions

`export_shape_controls.jl` exports five actual `Sector` geometries through the
existing Julia geometry owner: sharp 60/120 degree, offset 120 degree, and rounded
60/120 degree sections at 3 mm back radius. It does not reconstruct their contact
geometry in Python or alter the production owner. The successful export is
`q1-owned-shapes-v2`; the initial script spelling error is retained separately.

The native transfinite experiment makes a 0.1-scale internal copy around the
known centroid, connects corresponding boundary curves with triangular
transfinite strips, and keeps a small unstructured core. All strips and the core
remain in one conductor physical group. The original contour is unchanged;
artificial cuts are not terminal boundaries. A common radial count follows the
longest spoke and the prescribed first size/growth. This is a bounded construction
for the convex, single-contour fixtures, not a universal CAD partitioner or runtime
fallback.

The first experiment exposed noisy `gmsh.model.getCurvature` estimates on these
small built-in GEO arcs, causing needless subdivisions. Using their actual circle
geometry removed that experimental sizing error. The production owner already
has exact arc radii and should use them directly. The earlier inflated meshes
remain in `q1-owned-partitions`.

With exact arcs, all 15 combinations of five shapes and delta/back-radius
0.02, 0.2, 5 meshed successfully. Area errors ranged from 0.0632% to 0.0786%,
all interface-edge checks passed, and all triangle Jacobians were positive.
Maximum first normal heights were no larger than their prescribed targets.
At delta/radius=0.02, total local node counts were 956, 1600, 1302, 2937 and 3349.
These results are in `q1-owned-partitions-exact-arcs`.

Direct boundary-layer controls failed with `Edge not recovered` on four of the
five actual shapes. The offset sharp control returned a mesh but included a
negative signed triangle Jacobian. They are not qualified alternatives.
Ending layers at the sharp pie-sector vertices also failed (`q1-sector-ends`).

Eight sector diffusion pilots (`q2-sector-pilot`) checked sharp and rounded
120-degree sections at DC and delta/radius=0.02. Normal/tangential combined
refinement changed resistance by at most 0.0506% and complex impedance by at
most 0.0728%; loss balance was below 8.2e-14. This establishes promising internal
diffusion convergence only. The complete line-current problem, proximity loss,
material contacts and shared-terminal screen tests remain separate requirements.

The completed follow-up covers all five sectors at DC and delta/back-radius
5, 0.2 and 0.02 with independent normal/tangential refinements: 80 graded
solves in total, plus five strong-skin locally isotropic references. The largest
R change is 0.023660% for normal refinement at 96 segments, 0.203202% for
tangential refinement, and 0.028701% for normal refinement at 192 segments.
Each is below the 0.25% final-refinement target. All 85 meshes have conforming
interfaces and positive Jacobians; worst relative loss balance error is 4.23e-13.

At delta/back-radius=0.02 (approximately 1.213 MHz), the 96-segment, delta/6
candidate compares with the independent isotropic meshes as follows:

| Sector | Candidate nodes | Isotropic nodes | R difference | Complex Z difference |
|---|---:|---:|---:|---:|
| sharp 60 degrees | 956 | 218917 | 0.543570% | 0.420241% |
| sharp 120 degrees | 1600 | 437343 | 0.227717% | 0.199279% |
| offset 120 degrees | 1302 | 300693 | 0.317287% | 0.259256% |
| rounded 60 degrees | 2937 | 126155 | 0.007461% | 0.043418% |
| rounded 120 degrees | 3349 | 306377 | 0.004795% | 0.046197% |

The graded GetDP runs used 0.05–0.13 s and 21996–29300 KiB peak RSS; isotropic
controls used 9.48–48.62 s and 541176–1995280 KiB. Meshing took 0.009–0.031 s
and 8.23–37.41 s, respectively. These measurements include only the local
fixture, not the full exterior/PML or Julia compilation. Only one isotropic
level was run; its values are independent discretization references, not exact
sector solutions. The prescribed-boundary limitation above still applies.

## Evidence directories

All paths are beneath `.linecablemodels/fem/conductor-mesh-qualification/`:

- `q1-circles-annuli/`: initial API controls and measurements.
- `q1-shapes-roundtrip/`: native activation experiment and the all-corner fan failure.
- `q1-corners-interfaces/`: no-fan convex failure; successful reentrant, contact and screen controls.
- `q1-sector-controls/`: failed dense-radial/apex-fan variants.
- `q1-sector-partition/`: native transfinite alternative and round trips.
- `q1-pml/`: completed exterior preservation check.
- `q1-thin-and-physical/`: thin shells and bent strips; failed exact-ray probes retained.
- `q1-normal-probe/`: corrected probe for circular and physical-scale wire controls.
- `q1-fitted-depth/`: complete-layer depth prescription and native round trips.
- `q2-pilot/`: eight solved DC/1 MHz solid and tubular conductor controls.
- `q2-round-spectrum/` and `q2-round-summary/`: completed spectrum and derived report.
- `q1-owned-shapes-v2/`: actual Julia-owned sector geometry and exact areas/centroids.
- `q1-owned-partitions-exact-arcs/`: 15 successful convex-sector mesh controls.
- `q2-sector-pilot/`: eight solved sector diffusion refinement controls.
- `q2-round-normal-fine/`: 63 completed final normal-refinement controls.
- `q2-sectors-remaining/`, `q2-sectors-weak-transition/`, `q2-sectors-independent/`:
  72 additional graded sector controls, completing the 80-case matrix.
- `q2-sectors-isotropic/`: five completed independent strong-skin references.
- `q2-line-inputs/`: twelve native two-wire exports and matched analytical matrices.
- `q2-line-pilot/`: full-formulation mixed-pair pilot; inspect `runs.csv` for completion.

Each batch has an adjacent `.log`. Later batches also retain the exact script
used as `experiment.py`. The script prints fixture starts, failures and measured
completion summaries; successful prior cases are not rerun for progress reporting.

## Full-system pilot and remaining work

The serial batch recorded in `q2-followup.sh` completed all 140 new solves.
Its completed log is `q2-followup.log`; previous eight sector pilots were reused.

`export_line_controls.jl` exported twelve cases through the existing public
exporter: three placements and four soil resistivities, each with the eight
specified frequencies. Both analytical and FEM reductions are explicitly off.
The exterior preset remains domain_skin_depths=24, pml_layers=192,
mesh_size_factor=3 and exterior_mesh_size_factor=8.

`qualify_line_pair.py` changes only the round-conductor mesh in these native
exports. It uses fixed quarter-arc counts, conductor-restricted bulk sizing,
and the same fitted normal grading as the isolated candidate. Native expressions
read copper sigma, mu and frequency from the exported arrays. The global lower
size floor is removed so it cannot defeat conductor-local sizing; existing
exterior fields and transfinite PML prescriptions are retained. The harness
compares PML node coordinates and checks conductor conformity/Jacobians.

The first pilot is mixed placement, 100 ohm m, 10 kHz, with baseline, normal
(96 segments, delta/6), and combined_fine (192 segments, delta/12) in both
quasi-fw and quasi-tem: six solves. The first candidate mesh has 218647 nodes.
The existing detached driver owns paths, equations and matrix extraction.
Each successful solve appends its raw-result directory and measured costs to
`runs.csv`; `q2-line-pilot.log` records starts and completions. No maps are
requested in this first cost/control pilot; selected maps remain required.

Partial pilot evidence (first three of six solves complete): quasi-fw candidate
self-R errors against the matched analytical reference are approximately
0.209% (aerial) and 0.245% (buried), versus 1.091% and 1.073% on the baseline.
All four real and four imaginary Y components have the reference signs at this
single frequency on both meshes; this does not qualify crossings over the spectrum.
Candidate nodes increase by 1.44%, GetDP wall time by 8.64%, and GetDP peak RSS
by 2.81%. However, complete detached invocation time increases from 227.25 s to
366.79 s (61.4%); the additional time precedes/follows GetDP, so the total-time
budget is not yet established. The candidate ran first, and a warmed breakdown
of measurement-path preparation is still needed. Finer-control results and the
baseline quasi-TEM solve were still pending at that checkpoint.

At the next checkpoint, four of six solves were complete; the combined-fine
quasi-fw run was active. Baseline quasi-TEM used 94.54 s / 1759232 KiB versus
82.47 s / 1781876 KiB for the candidate. Retained `paths.pro` shows the aerial
receiver has 24 quadrature samples on the baseline and 192 on the candidate.
Existing detached `write_paths` calls `path_points` separately for each sample;
that routine converts the entire triangle array and recomputes its bounding
boxes on every call. Non-GetDP elapsed time rises from 16.08 s to 137.38 s,
consistent with this eightfold increase in path preparation. This is a source-
and scaling-based diagnosis, not yet a profiled attribution of every second.
Do not attribute the extra time simply to first-use compilation. Preparation
cost/parity must be measured before the full runtime gate can pass; do not
reduce contour accuracy or change voltage quadrature to hide this cost.

The pilot subsequently completed all six solves. The candidate-to-combined-fine
self-R changes are 0.002015% (aerial) and 0.000909% (buried) in both formulations.
The largest complex Y change is 0.014450% for quasi-fw and 0.013421% for quasi-TEM.
All four real and four imaginary Y entries retain the matched reference signs
on all three meshes. This is one frequency/placement/resistivity, not spectrum
qualification. The remaining approximately 0.21–0.25% self-R disagreement is
stable under this conductor refinement; further conductor resolution does not
remove it. See `q2-line-pilot-summary/matrices.csv` and `refinement.csv`.

| Mesh | Nodes | Quasi-fw GetDP time / RSS | Quasi-TEM GetDP time / RSS |
|---|---:|---:|---:|
| baseline | 215551 | 211.17 s / 4019656 KiB | 94.54 s / 1759232 KiB |
| normal | 218647 | 229.41 s / 4132612 KiB | 82.47 s / 1781876 KiB |
| combined_fine | 226614 | 371.27 s / 4310300 KiB | 94.71 s / 1798308 KiB |

The fine quasi-fw invocation took 617.95 s including preparation. Its higher
GetDP time includes more voltage-path postprocessing, not only solving the
larger matrix. The normal level remains the candidate; the fine level is a
qualification control, not a proposed default.

## Path-preparation qualification and native integration

`q2-path-bounds` completed twelve measurements: the six retained pilot meshes /
formulations, each called twice. `qualify_path_bounds.py` reuses triangle arrays
and bounds within one path-writing call, filters candidates without changing
their order, and delegates clipping to the existing implementation. Every output
is byte-identical to the original `paths.pro`. Candidate quasi-fw preparation
took 15.20 s on first call and 14.74 s on the repeat; combined-fine took 29.23 s
and 28.76 s. These are Python preparation measurements, not Julia compilation or
a complete remeasured compute time. The experiment remains in the harness and
is selected explicitly with `--reuse-path-bounds`; production is unchanged.

The user asked whether GetDP can do this integration natively. A small independent
complex H(curl) projection in `q1-native-line-integral` confirms that it can
integrate a field along conforming physical line groups using
`Integral { [{b}*Tangent[]]; ... Jacobian Lin; Integration ...; }` and
`Print[Circulation[Path], OnRegion Path, ...]`. The line groups must be included
in the edge field's support. With only volume support the test returned zero;
adding the existing shared mesh edges to the support preserved the DOF count.
`Trace` is rejected in postprocessing by the installed GetDP 3.5.0; that failed
attempt is retained separately.

For the exact field `(1+2j)*(3-0.7*y, -2+0.7*x)`, the native vertical-path
integral is `-1.649999999999999-3.299999999999999j` versus exact `-1.65-3.3j`.
The counterclockwise rectangle circulation is
`2.800000000000001+5.600000000000001j` versus exact `2.8+5.6j`.
This small 232-DOF test establishes the native capability, not equivalence of
the full cable extraction. Current receiver paths cut across volume elements;
they are not physical mesh-line groups. A native replacement must preserve those
paths, averaging, metal exclusions and PML convention, and qualify the resulting
mesh/field support and extracted P/Y. Uniform `OnLine` sampling alone does not
provide element-crossing quadrature. Do not substitute it silently.

Sources: [GetDP postprocessing](https://getdp.info/doc/texinfo/getdp.html#PostProcessing),
[postoperations](https://getdp.info/doc/texinfo/getdp.html#Types-for-PostOperation),
and [tangent function](https://getdp.info/doc/texinfo/getdp.html#Miscellaneous-functions).

The user subsequently clarified that preparation is to be retained, including
the byte-preserving bounds reuse. The native-integration question was not a
request to change extraction. The small native experiment is informational;
replacing preparation is not a gate or an implementation task in this campaign.

## Completed initial full spectrum — 2026-09-28

`run_line_spectrum.sh` starts 190 new candidate solves in `q2-line-spectrum`,
serially, using the byte-preserving path-bounds experiment in the harness. It
reuses the two completed normal-level pilot solves at mixed/100 ohm m/10 kHz.
Together these cover the planned 192 cases. Both formulations write selected
field maps for mixed placement, 100 ohm m, 1 MHz. The native extraction equations
remain unchanged; replacing them with native mesh-line integration requires its
own equivalence check and is not mixed into this conductor-mesh comparison.

All 190 new solves finished successfully at 04:47:57 local time on 2026-09-28.
Together with the two retained pilot results, all 192 cases are complete, unique
and contain finite raw Z/P/Y. No equations or production files changed during
this qualification: all 298 captured `src/` and `ext/` files match Q0.

The outer `q2-line-spectrum.log` contains meshing/stage updates. Actual GetDP
output is appended across every case to `q2-line-spectrum/solver.log`, which can
be followed with `tail -f`. Completion rows are flushed to each case's `runs.csv`.
No agent polls the batch while waiting for the user to resume. All 298 captured
production files still match Q0 at launch.

`summarize_spectrum.py` reads these artifacts without solving and writes
`q2-line-summary/{matrices,signs,self_resistance,costs,algebraic_residuals,crossing_brackets}.csv`.
The first conductor of the mixed fixture is aerial; the second is buried.
Each row below represents 32 frequency/resistivity cases and 256 individual
real/imaginary Y entries. All reductions are off on both references.

| Physics | Placement | Raw G/B signs matching analytical | Worst self-R relative error |
|---|---|---:|---:|
| quasi-fw | both aerial | 256/256 | 0.339459% |
| quasi-fw | mixed | 256/256 | 0.698626% |
| quasi-fw | both buried | 256/256 | 0.675943% |
| quasi-TEM | both aerial | 254/256 | 0.339459% |
| quasi-TEM | mixed | 245/256 | 0.698626% |
| quasi-TEM | both buried | 256/256 | 0.675943% |

These are sampled sign comparisons, not a certification of signs below numerical
uncertainty or between samples. The analytical mixed and buried mutual terms
also cross zero. The 65 component/model/formulation crossing records collapse
to 11 distinct placement/resistivity/frequency intervals. They need refinement
before any crossing conclusion. No sign is forced, clipped or repaired.

The aerial quasi-TEM differences are the two mutual G entries at 1000 Hz,
rho=0.1 ohm m: approximately +9e-19 S/m in FEM versus -4.38e-16 S/m in the
reference. The mixed quasi-TEM differences concern G12 (and at higher
frequencies B12), reaching non-negligible magnitudes. For example at rho=0.1,
100 kHz, FEM G12=-5.4910e-8 S/m versus analytical +1.7019e-5 S/m. This cannot
be classified as tiny roundoff. Compare the original and finer meshes before
attributing it to the new conductor discretization. Quasi-TEM model equivalence
is not assumed by this mesh feature.

The worst self-R discrepancy is the buried wire in the mixed rho=0.1 case at
100 kHz (0.698626%). The earlier 10 kHz/rho=100 pilot is already stable under
conductor refinement; the other cases still need their selected controls.
384 logged basis residuals were retained across both formulations. The largest
residual/RHS norm is 2.23068e-6 (quasi-fw, both buried, rho=100, 1 Hz).
This global algebraic residual is not a componentwise uncertainty bound on Y.
Per-case DOFs, wall times and peak GetDP memory are retained in `costs.csv`.
This initial spectrum alone does not establish the baseline-relative warmed
runtime gate, integrated-loss gate or native-edit parity.

`plot_spectrum.jl` produced 18 R/G/B matrix pages (PNG and SVG) in
`q2-line-summary/plots`, using public `LineParameters`, `ObservedResult(...;
clip=false, length_unit=:base)` and `LineCableModels.plot`. The script defaults
to interactive GLMakie. `FEM_PLOT_BACKEND=cairo` was explicitly selected for the
saved report. The rendered aerial and mixed quasi-fw G pages were inspected;
tiny values remain visible and sampled crossings are retained. The ordinary
manual FEM runner is unchanged.

The retained mixed/rho=100/1 MHz/quasi-fw/basis-1 maps were also rendered using
the existing `plot_fem_pml_fields.py`, with no additional solve. Magnitude,
phase, material ownership and an unaveraged vertical cut are in
`q2-line-summary/fields-quasi-fw-basis1`. The magnitude page was inspected;
physical and PML panels have separate colour scales. This rendering is not an
integrated-loss or field-convergence test. Both bases and both formulations'
native selected maps remain available in the original runs.

## Fixed crossing and preservation follow-up

`export_line_followup.jl` exports six selected cases through the same public
exporter. It contains 11 logarithmic bracket midpoints plus 3.16228 Hz near the
conductor skin transition for each placement at rho=0.1. Original frequencies,
baseline meshes and completed normal results remain reusable.

`run_line_followup.sh` prescribes 108 new solves, serially:

- 38 combined-fine solves at all 19 distinct crossing endpoints, both physics;
- 56 normal/combined-fine solves at the 14 intermediate/transition frequencies;
- 14 original-mesh quasi-TEM controls at the anomaly interval endpoints.

The finer conductor mesh uses 192 contour segments and delta/12; the candidate
uses 96 and delta/6. Each frequency retains the same prescribed exterior across
its mesh comparisons. All use the existing preparation and byte-preserving bounds
reuse. This is a fixed qualification schedule, with no retries, fallback or
engine-level acceptance logic. The outputs are `q2-line-followup-inputs` and
`q2-line-followup`; GetDP appends to `q2-line-followup/solver.log`.

Read componentwise real/imaginary changes before deciding whether each sign is
resolved. An unresolved root bracket may need more explicit test frequencies;
this first follow-up does not promise an exact zero location. Existing quasi-TEM
differences must be preserved/converged or explained separately. The native
integration experiment above is not part of this schedule.

## Completed crossing follow-up — 2026-09-28

The IDE interruption left 79 completed follow-up records. The explicit recovery
script `q2-line-followup/resume-20260928.sh` ran only the remaining 29, in fresh
subdirectories while preserving the original results and interrupted artifacts.
All 108 follow-up controls completed at 10:50:02 local time on 2026-09-28.

`summarize_line_followup.py` combines these with the initial spectrum and all
six pilot solves: 304 unique solves, or 912 finite Z/P/Y matrices. It writes
raw tables, refinement changes, sign comparisons, baseline values, costs and
PML coordinate comparisons to `q2-line-followup-summary`.

- All 68 normal/combined-fine pairs retain the PML coordinates to 1e-10 m.
- Across the initial and added frequencies, quasi-fw matches all 880 sampled
  analytical G/B signs. All 272 entries with a paired finer mesh also match.
- The maximum change in self-R between normal and combined-fine meshes is
  0.140194065%, for the second buried conductor, rho=0.1, 1 MHz. This is a
  conductor-refinement measurement, not the analytical reference error.
- Quasi-TEM matches 863/880 normal-mesh signs and 255/272 paired-fine signs.
  All 13 original anomalous entries have the same anomalous sign on the
  original, normal and finer meshes. Four further differences occur at new
  midpoint frequencies. The original anomalies are therefore not introduced by
  conductor grading. Their values and changes remain in `signs.csv`.

For example, mixed rho=100 at 1 MHz has quasi-TEM G12 values
`+1.267573e-6`, `+1.266605e-6`, and `+1.266536e-6 S/m` on baseline, normal,
and combined-fine meshes, against analytical `-1.135455e-5 S/m`.
This persistent formulation/reference difference is separate from conductor
mesh accuracy. No clipping, sign repair or reciprocity enforcement is applied.

The paired FEM values and reference values all exceed the observed componentwise
conductor-refinement changes in these comparisons. This is not a rigorous total
error bound: exterior discretization, roundoff and algebraic conditioning have
not been bounded by that test. Root positions are bracketed, not exactly located.

The public API plots were regenerated with the actual per-case frequency grids
in `q2-line-followup-summary/plots` (18 PNG and 18 SVG pages, `clip=false`). No
zero-valued placeholders are inserted for frequencies absent from another case.
The mixed quasi-fw conductance page was inspected, including the additional
crossing samples shared with the analytical reference.

## Coupled screen qualification

`qualify_screens.py` and `coupled_conductors.pro` are test-only local magnetic
controls. They use the A_z/terminal-voltage weak form of the maintained axial
formulation, with displacement and the separate transverse problem omitted.
Copper conductors are surrounded by lossless air; A_z=0 on a circular exterior.
A 6 mm radius core carries +1 A; all explicit 3 mm diameter screen strands share
one terminal voltage and a total -1 A return current. The conductor interiors
remain resolved. This exercises proximity/current sharing and native integrated
losses, unlike the earlier prescribed-Ez diffusion test. It does not replace
full layered Z/Y or arbitrary-material qualification.

The fixed-gap fixtures keep 0.5 mm between adjacent strands while changing
24/48/96 strand count. The separate fixed-envelope fixtures keep the strand
centre radius at 60 mm. All use the same 250 mm exterior radius. Normal and
combined-fine levels retain the existing candidate parameters. Isotropic
references use native Distance/MathEval/Restrict/Min fields within the metals,
with boundary spacing <=delta/6 or delta/12 and gradually larger interior
elements; no isotropic refinement is imposed throughout the air. The coarse
12-segment control reports overhead against a coarse local benchmark, not
warmed production-engine performance.

The shared experimental `configure()` previously replaced the background Min
for each bulk-restricted region, leaving only the last active. Before these
multi-wire bulk tests it was corrected to combine all restrictions in one Min.
Earlier single-region electrical controls are unaffected. This campaign's
change is only in `test/manual/`; it makes no production edit.

Two completed setup controls establish the current-driven benchmark:

| Control at 1 MHz | Nodes | GetDP time / peak RSS | Loss/power relative difference | Current constraint error |
|---|---:|---:|---:|---:|
| single 3 mm wire | 2301 | 0.05 s / 26004 KiB | 4.60e-11 | 5.29e-13 A |
| core + 24 strands, fixed gap | 38898 | 2.69 s / 224360 KiB | 2.00e-9 | 3.12e-11 A |

The single-wire resistance is 0.02835979 ohm/m versus the independent cylindrical
value 0.02830136 ohm/m: 0.206455% error. The 24-strand loop resistance is
0.01354179 ohm/m, with integrated loss 0.01354179 W/m at unit RMS drive.
Its normalized algebraic residual is 9.89e-14. All conductor area errors are
below 0.1%, shared edges conform and triangle Jacobians are positive. This pilot
does not yet establish screen convergence or relative cost against a reference.
The pilot's optional `jz.pos` was initially emitted as a GetDP Table; subsequent
map output uses a separate native Gmsh-format postoperation. The numerical
metrics and weak form are unchanged by that format correction.

`run_screens.sh` prescribes 56 new serial controls in `q2-screens`, reusing the
two pilots. It covers the independent single-wire checks, both screen layouts,
50 Hz/10 kHz/1 MHz, conductor refinement, isotropic references and coarse
controls. Actual GetDP output is appended to `q2-screens/solver.log`. Its
qualification results must be assessed before production implementation.

At the subsequent source check on 2026-09-28, independent live-checkout changes
were detected in 13 captured production files, with `voltage_paths.jl` and
`onelab_export/measurements.py` removed. These edits were not made by this
campaign and are left untouched. Their current hashes are recorded in
`q0-20260928-live-drift.csv`. The completed full-line results refer to the
captured/exported solver versions, not automatically to this updated backend.
The local screen benchmark is self-contained and can proceed; production
integration/native parity must be reconciled with the live backend afterward.

## Screen reference ordering diagnosis — 2026-09-28

The first isotropic 24-strand/1 MHz reference reached 224562 DOFs, then emitted
no output after `Solve[S]` at 11:21:37. The process remained on one CPU core,
with approximately 613 MiB resident memory and no process swapping. After over
35 minutes, two debugger samples located it in
`findIndMultisecs -> shrinkDomainDecomposition -> constructSeparator ->
SPACE_ordering -> mumps_pord -> MatLUFactorSymbolic_AIJMUMPS`.
The bottleneck was PORD's ordering/analysis, before numerical factorization.
CPU activity alone was insufficient evidence of acceptable progress.

The cumulative and per-case logs were live: the harness writes and flushes each
received line to both. GetDP/MUMPS was emitting no new lines. Only this campaign's
identified GetDP process was terminated; eight completed records, the unfinished
mesh and all other FEM jobs were preserved.

An explicit test in `q2-screen-ordering/isotropic-amd` copied the mesh and `.pro`
byte for byte, then requested AMD using `-mat_mumps_icntl_7 0`. Native MUMPS
messages were enabled with `-mat_mumps_icntl_3 6 -mat_mumps_icntl_4 2`.
The solve completed in **47.39 seconds**: MUMPS reported 8.0002 s analysis,
0.4092 s numerical factorization and 0.0644 s solve. The reported effective
ordering was AMD. Peak GetDP RSS was 1015304 KiB. Loss/power mismatch was
3.47e-9 relative, maximum current error 5.24e-11 A and normalized residual
1.87e-11.

Same-mesh comparisons on the already solved normal and combined-fine controls
changed R by only 2.33e-11 and 5.76e-11 relative, respectively; complex loop Z
changed by less than 8.5e-13 relative. Their AMD invocations took 1.56 s and
8.10 s. These are local same-mesh solver measurements, not a claim about warmed
whole-engine or Julia performance.

The recovered isotropic result is indexed by
`q2-screens/isotropic-amd-recovery/measurements.csv`; its run directory points
to the retained ordering experiment. Its Python mesh timing is recorded as NaN
because the interrupted process had not flushed that timing. The original Gmsh
log records 8.05631 s meshing. No finished solve was discarded.

`qualify_screens.py` now accepts an explicit `--mumps-ordering` and emits native
MUMPS phase messages. The fixed schedule selects AMD. This is a prescribed
qualification setting, with no retry or automatic solver fallback; no production
solver default changed. `q2-screens/resume-amd.sh` runs the remaining 47 cases
serially and appends to the same `q2-screens/solver.log`.

The screen accuracy gate remains open: at this point loop R is 0.01354179,
0.01246196 and 0.01282818 ohm/m for normal, combined-fine and isotropic meshes.
The finer isotropic and remaining controls must distinguish conductor/proximity
discretization effects. A small algebraic residual or loss balance does not
establish mesh convergence.

Native option references: [PETSc MUMPS controls](https://petsc.org/release/manualpages/Mat/MATSOLVERMUMPS/)
and [GetDP solver options](https://www.getdp.info/doc/texinfo/getdp.html#Frequently-asked-questions).

Q1 still needs qualification of partition integration with contacts/composed
regions, layer-to-bulk transition and native editable expressions. Q2 needs
completion/assessment of shape/proximity/screen controls and cost measurements
and native remeshing parity. No production implementation may precede these gates.

## Fine screen reference assembly diagnosis — 2026-09-28

The 415092-DOF isotropic-fine reference exposed a separate assembly bottleneck,
before MUMPS. An earlier invocation exited with status 143; its termination
source is unknown. On explicit restart, GetDP remained alive at `Generate[S]`
80%, using one CPU core, approximately 1.06 GiB RSS and no process swap. Two
debugger samples showed `memmove -> MatSetValues_SeqAIJ ->
LinAlg_AddComplexInMatrix -> Dof_AssembleInMat`, rather than symbolic ordering.

Inspection of the unchanged mesh found 415282 finite-element nodes and 830370
triangles. Its conductor-local isotropic sizing leaves large triangle fans in
air: 130 nodes have more than 90 incident triangles, with a maximum of 955
(956 adjacent nodes including self). GetDP's ordinary-row preallocation is 100;
the two terminal rows are already marked nonlocal in the saved `.pre` file.
The local air rows exceed that capacity and force expensive PETSc matrix
storage reallocations. This is an execution problem of this reference mesh,
not evidence of a physical sign change or an adequate exterior discretization.

The explicit native option `-petsc_prealloc 1024` is tested on byte-identical
normal-control mesh/equations in `q2-screen-preallocation/normal-1024`.
The maximum absolute difference across its terminal voltage, current, loss and
area records is 3.54e-14 (mixed units; not a relative scientific error bound).
The control completed in 2.05 s with peak RSS 924080 KiB. This setting trades
memory for avoiding repeated matrix copies; it is not an engine performance
claim, a mesh refinement, or a production solver default.

`qualify_screens.py` exposes the prescribed `--petsc-prealloc` option and records
it with each result. Native GetDP progress now uses `-p 1 -cpu`. The explicit
remaining-case schedule `q2-screens/resume-prealloc-20260928.sh` preserves nine
completed controls, selects 1024 for the remaining isotropic references and
keeps AMD ordering. It has no adaptive sizing, retry, or fallback. The ordinary
graded cases retain their previous row preallocation. All interrupted
directories remain intact. Air-transition accuracy and the final screen
convergence/cost gates remain open.

The restarted 415092-DOF reference passed assembly and reached `Solve[S]`
at GetDP wall time **16.6317 s**. Native MUMPS output confirmed AMD analysis
of 3733212 nonzeros. At this startup check the process used 3941984 KiB RSS
and zero swap. The reference solve and remaining serial batch were left running;
this observation is not a claim that the reference accuracy gate has passed.

## Completed screen batch and isolated air-mesh defect — 2026-09-28

The serial screen batch finished at 13:03:16 local time: all 56 scheduled
controls plus the two retained pilots are available, with no duplicate completed
case keys. `q2-screen-summary` contains the raw merged table, 22 refinement pairs,
the single-wire references, and inspected PNG/SVG resistance and cost plots.
All terminal results and diagnostic quantities are finite. Maximum loss/power
mismatch is 4.72e-8, current-constraint error 1.62e-10 A, and normalized algebraic
residual 7.72e-11. These establish solved discrete systems, not mesh accuracy.

Only 10 of the 22 normal/combined-fine pairs meet the <=0.25% R-change criterion:
all four singleton frequencies and the six 50 Hz screens. At 1 MHz the fixed-gap
24/48/96-strand R changes are 8.665%, 4.399% and 2.895%; fixed-envelope changes
are 1.915%, 2.451% and 2.394%. The 24-strand isotropic/isotropic-fine values differ
by 15.47% relative to the finer value, so they cannot yet serve as a qualified
accuracy-matched reference. No target was relaxed.

The 415092-DOF finer reference ultimately completed in 75.70 s with peak RSS
4354224 KiB. The largest retained case, the 96-strand isotropic reference, used
791725 generated nodes, 274.65 s and 9482336 KiB peak GetDP RSS. These are
observed standalone costs with recorded ordering/preallocation settings, not
qualified whole-engine performance ratios or equivalent-accuracy speedups.

The qualification fixture used conductor-local size fields, disabled boundary
size extension, and prescribed a 20.833 mm air bulk cap around 3 mm wires with
0.5 mm gaps. This generated abrupt triangle fans. The preallocation correction
removed their execution bottleneck but did not improve their approximation.

A bounded experiment in `q2-screen-air-pilot` adds Gmsh's native `Extend` field
with `Power=1`, `SizeMax=outer_radius/12`, and
`DistMax=outer_radius/(12*air_growth)`, restricted to the air surface and excluding
its boundary. The pilot uses `air_growth=0.25`. This grades outward from actual
contour-edge lengths while retaining the remote size cap and domain/BCs.
It does not modify skin-depth sizes or conductor geometry. The construction
uses the documented [native Extend field](https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-mesh-size-fields).

An exact comparison of the saved old/new meshes, after canonical node
renumbering, found identical conductor coordinates and triangle connectivity:
37864 conductor nodes and 73278 conductor triangles. The identical canonical
hash and per-mesh incidence counts are in `conductor-preservation.csv`.
The initial comparison by raw node IDs failed because 196 IDs moved; no numerical
tolerance was introduced to obtain the coordinate/connectivity equality.

| 24 strands, fixed gap, 1 MHz, normal conductor mesh | Previous air | Extended air |
|---|---:|---:|
| Loop R [ohm/m] | 0.01354179044 | 0.01081521277 |
| Loop X [ohm/m] | 0.8632214021 | 0.8930010654 |
| Generated nodes | 38898 | 49241 |
| Maximum incident triangles at a node | 40 | 9 |

Changing only the air mesh changes R by **20.1345%** relative to the previous
value. This confirms an important air-discretization error contribution in this
benchmark; it does not establish the new value as converged. The pilot used
0.602 s meshing, 2.04 s GetDP and 214224 KiB peak GetDP RSS; those are individual
observations, not a warmed performance comparison. The native current-density
map and an inspected old/new air-mesh figure are retained with the pilot.

The next fixed schedule, `run_screen_air.sh`, has 17 new controls on this same
24-strand fixture, reusing the completed pilot. At 1 MHz it varies normal and
tangential conductor refinements independently at air growth .25 and .125,
compares two isotropic levels at both air settings, and adds a third air level
.0625 for normal/combined-fine conductors. Four additional controls check 50 Hz
and 10 kHz. Conductor mesh preservation must be audited across air-only pairs.
The script runs serially with ordinary logs; it does not self-refine or retry.

This is a correction to the qualification fixture. The production exterior
owner and full-line Y equations remain untouched. The screen accuracy/cost gates
and the remaining shape/native parity gates are still open.

## Air-refinement results and reference mesh correction — 2026-09-28

The 17 follow-ups finished at 13:35:36 local time. Together with the retained
pilot, they give 18 controls in `q2-screen-air-summary`. For the 24-strand
fixed-gap fixture at air growth .125, normal/combined-fine R changes meet .25%
at all three tested frequencies. At 1 MHz the change is **0.161178%**; independent
normal refinement at 96 contour segments changes R by approximately 0.113%,
and doubling contour resolution at the same normal grading changes R by
approximately 0.0482%. These are observed mesh changes, not absolute error bounds.

On the normal conductor mesh, reducing air growth .25 -> .125 changes R by
0.072705%, and .125 -> .0625 by **0.027983%**. Canonical comparison of all saved
meshes confirms exact conductor coordinates/connectivity across air growth for
normal, isotropic and isotropic-fine levels. It does not for normal-fine,
combined and combined-fine layered reconstructions: those have additional
interior tessellation changes. The failed blanket-preservation assertion is
retained in the audit rather than weakened. For example, combined-fine at .125
has eight more conductor nodes and 16 more triangles than at .25. Such pairs
cannot be described as strictly air-only. The exact normal-level comparison
still isolates the observed air effect without that ambiguity.

The independent reference gate remains open. At air growth .125, isotropic and
isotropic-fine R are 0.01079758548 and 0.00961124245 ohm/m. The discrepancy lies
in the core loss (0.00696831 versus 0.00578692 W/m at unit RMS drive), while
screen loss agrees much more closely. Inspection found **67302 core triangles
in isotropic but only 9426 in isotropic-fine**. Region areas and terminal current
constraints remain correct; those checks did not expose missing skin resolution.

A mesh-only experiment reopens the saved native geometry and meshes only the
visible core with the same prescribed fine size field. Frontal-Delaunay
(algorithm 6) produces 9548 core triangles and 95th-percentile radial span
0.580756 delta among triangles within two delta of the surface. Native Delaunay
(algorithm 5) produces **148514 triangles**, with that span reduced to
**0.257622 delta** and maximum node incidence 9. Meshing took 0.2485 and 1.9750 s.
This isolates an algorithm/size-field failure in the reference construction;
the mesh-only counts do not claim byte-identical reproduction of the complete
fixture. [Gmsh tutorial t10](https://gmsh.info/doc/texinfo/gmsh.html#t10) explicitly
identifies Delaunay as better suited to strong size-field gradients.

The separate attempt to locally remesh a saved MSH2 reference failed before any
solve: that file omits unphysical conductor-boundary line elements, and the
GEO+MSH import lost the curve classification needed by the remeshing operation.
Gmsh repeatedly reported edge-recovery warnings. The process was stopped; its
script and logs remain under `q2-screen-reference-mesh/core-only-delaunay`.
No result was accepted and the unusable helper was removed from the manual
entrypoints. This is not a production fallback or change to external-mesh use.

`run_screen_air_spectrum.sh` prescribes 42 new serial controls: 30 missing graded
controls over the remaining 24/48/96-strand layouts and frequencies, plus both
native-Delaunay isotropic reference levels for all six layout/count combinations.
It reuses six completed fixed-gap/24-strand graded controls. Air growth stays
.125 and AMD ordering stays explicit. The corrected air field gives low node
incidence, so these references use the ordinary row preallocation. The field,
algorithm and equations are prescribed once; no acceptance-driven retry or
self-refinement runs in this batch or in production.

## Power-outage recovery of the final screen reference — 2026-09-28

After the host reboot, 41 of the 42 prescribed `q2-screen-air-spectrum`
controls were complete. All 41 CSV records were checked against finite saved
terminal results, matching loop R, nonempty mesh/solution/cost artifacts and
completed GetDP logs. The last completed case was the fixed-envelope 96-strand
isotropic reference at 1 MHz, saved at 16:21:29 local time. The fixed-gap
96-strand isotropic-fine reference had also completed, at 15:00:11, taking
1893.48 s in GetDP. None of these completed results needs rerunning.

Only `fixed_envelope-n96-f1e+06-isotropic_fine` was interrupted. Its directory
contains the native geometry but no saved mesh or solver result: the outage
occurred during meshing. No qualification process survived the reboot.
The temporary Python Gmsh installation was also lost; the same pinned
Gmsh 4.15.2 wheel was restored under `/tmp/lcm-screen-recovery-env` and native
initialization verified. NumPy 2.4.6, mpmath 1.3.0 and GetDP 3.5.0 remain in use.

The retained artifact script `q2-screen-air-spectrum/resume-after-power-outage.sh`
starts just that final case, with unchanged parameters and the captured harness
sources. The four captured source files match their live counterparts byte for
byte. Recovery writes into `q2-screen-air-spectrum/power-recovery-20260928`,
preserving the interrupted directory and all original CSVs. The original outer
log and `solver.log` are appended, so the existing live-log command remains valid.
Meshing startup was checked once; the ordinary serial process was left running.
Final assessment must merge the recovery CSV with the original 41 records and
the six reused graded controls. Production feature code remains untouched.

A second host reboot interrupted the recovery during meshing. At 17:04 local
time the host uptime was under three minutes and no qualification process was
running. All original 41 completed CSV records still matched their saved finite
terminal results and completed solver logs. The first recovery has only native
geometry, no mesh or result; there is no saved numerical stage to resume inside
that final case. The previous-boot kernel journal was not readable by this user,
so the reboot cause has not been established.

`q2-screen-air-spectrum/resume-after-second-reboot.sh` restarts only that same
case into `power-recovery-20260928-b`, retaining both interrupted directories.
Pinned Gmsh 4.15.2 is now installed under the qualification artifact root's
`python-env`, outside `/tmp`, so a reboot does not erase the dependency. Native
initialization and the unchanged captured harness were checked before launch;
the environment details are saved in `recovery-environment.txt`. The recovery
still uses one ordinary foreground process and appends the existing two logs.
Meshing startup was confirmed once and the agent stopped polling. No production
feature changes or scientific parameter changes were made.

## Completed corrected screen qualification — 2026-09-28

The final recovery finished at 17:46:09 local time. The non-solving
`summarize_screen_air.py` combines all 42 completed new controls and six reused
graded controls, verifies unique prescribed case keys and saved terminal/log
consistency, and writes `q2-screen-air-spectrum-summary`. This retains the raw
48-row table, 18 refinement comparisons, 18 fine-reference comparisons, observed
costs, and PNG/SVG resistance/reference/cost figures. Interrupted directories
remain intact and contribute no results.

All 18 normal/combined-fine R changes meet 0.25%; the maximum is **0.161178%**.
The six coarse/fine isotropic reference pairs change R by at most **0.142800%**.
Candidate-to-fine-reference R differences are **0.15994–0.20971%**; combined-fine
differences are **0.02983–0.04845%**. These are numerical reference differences,
not exact-error bounds. Maximum region area error is 0.0713794%, relative
loss/power mismatch 1.42e-7, current-constraint error 3.91e-10 A, and relative
algebraic residual 1.21e-10. Thus the completed screen evidence now supports the
graded construction; the previous failure is not silently waived.

| Fixed-envelope 96-strand screen, 1 MHz | Normal graded | Combined fine | Fine isotropic reference |
|---|---:|---:|---:|
| DOFs | 234973 | 715510 | 3518369 |
| Gmsh meshing [s] | 4.73 | 11.98 | 256.71 |
| GetDP execution [s] | 10.20 | 85.25 | 2071.92 |
| GetDP peak RSS [GiB] | 0.938 | 2.848 | 13.546 |

These observations exclude some preparation and do not establish warmed
whole-engine ratios. The expensive isotropic references remain harness-only.
The fixture uses 12 mm diameter copper core, 3 mm diameter explicit screen
wires, +1/-1 A terminal currents, air and an A_z=0 boundary at 0.25 m radius.
It does not qualify earth/PML admittance or quasi-TEM sign behavior. Remaining
work includes coupled noncircular shapes and material contacts, full-formulation
shape placement controls, native editable remeshing parity, and fair warmed
cost comparisons. Existing full-line evidence must also be reconciled with
the concurrently updated native voltage-extraction owner before release.

## Coupled-current sector proximity pilot — 2026-09-28

The next harness `qualify_sector_proximity.py` reuses the captured Julia-owned
sector contours, experimental transfinite partition function and maintained
local screen-control A_z/u equations. Two reflected copies face across a
0.5 mm gap; their currents are prescribed +1/-1 A and their terminal voltages
are solved. Air uses the native Extend field, growth .125 and outer radius
25 mm with A_z=0. The test checks that internal partition edges disappear from
the combined terminal contour and that every material interface remains conforming.
It does not prescribe Ez around the conductor surface, unlike the earlier
isolated sector diffusion controls. No full-line or contact claim follows.

The sharp 120-degree pilot at delta/back-radius .02 (1.21311 MHz) completed
normal and combined-fine controls, with native Jz maps, in
`q2-sector-proximity-pilot`. R is 0.08983001708 and 0.09009111530 ohm/m,
a **0.289816%** change relative to the finer result, above the .25% criterion.
Both meshes preserve positive Jacobians and conforming interfaces. Their node
counts are 4853 and 13601; GetDP takes .53 and 1.70 s with peak RSS 36048 and
73200 KiB. These maps include a field-output cost. Loss balance is below 1.68e-11
and current-constraint error below 1.74e-12 A; neither qualifies discretization
accuracy. The pilot is recorded as an unresolved refinement result, not accepted.

Four prescribed controls in `q2-sector-proximity-controls` separate normal-only
and tangential-only refinement and compare two independent locally isotropic
Delaunay meshes. They retain the same geometry, material, frequency and BCs.
This bounded diagnosis precedes a broader shape campaign; it does not change
production defaults or relax the criterion. The completed pilot is reused.

All four diagnostic controls completed. Normal refinement at the nominal
96-segment full-circle resolution changes
R by only **0.016808%**, while tangential refinement at the same normal grading
changes R by **0.305053%**. At nominal 192 resolution, halving the normal step changes R
by **0.015284%**. The two independent Delaunay references differ by **0.139284%**;
the nominal-192 delta/6 and delta/12 values differ from the fine reference by
**0.027076%** and **0.042354%**, respectively. Actual contour counts, including
straight sides, are retained in `regions.csv`; the nominal resolution is not a
claim that every sector has 96 or 192 boundary edges. This pilot identifies tangential
resolution as the main contribution to the initial .290% difference. It supports
testing the nominal-192 construction for coupled sectors; it does not change the
96-segment screen/wire prescription or prove a universal sector error bound.

`run_sector_proximity.sh` schedules 84 new serial controls, reusing these six
completed strong-skin sharp120 results. Four independent graded levels cover
all five captured sectors at DC and delta/back-radius 5, .2 and .02. Two
independent isotropic levels cover each shape at .02, giving 90 total controls
including reuse. The preset levels are requested once and no scientific
acceptance decision drives retries. Remaining material-contact, tube/strip,
native-edit and full-line/whole-engine cost requirements are unchanged.

The serial sector batch was launched into `q2-sector-proximity-spectrum`, with
outer output in `q2-sector-proximity-spectrum.log` and native GetDP output in
`q2-sector-proximity-spectrum/solver.log`. Startup was confirmed through completed
DC/AC records and the subsequent sector meshes. The agent then paused under the
agreed user-monitored execution protocol. No completed screen or sector-pilot
case is repeated in this batch.

## Normalized native three-wire comparison and production integration — 2026-09-28

The 18-case native selection completed: three prescribed frequencies (0.1,
10 k and 1 M Hz), baseline/normal/combined-fine meshes, both formulations and
all three source columns. The existing mixed `air_1` fixture uses radius
42.5 mm, soil resistivity 100 ohm m, domain size 24 skin depths, 192 PML layers,
bulk size factor 3 and exterior factor 8. The current GetDP edge-trace extraction
is retained. There is no voltage preparation stage.

`summarize_native_wires.py` verifies finite retained matrices and writes
`summary.csv` / `matrices.csv` beneath
`.linecablemodels/fem/conductor-mesh-qualification/native-three-wire/`.
Normal versus combined-fine results are:

| Frequency | Physics | Maximum self-R change | Maximum complex Y-entry change | Normal/baseline solve time | Normal/baseline peak RSS |
|---:|---|---:|---:|---:|---:|
| 0.1 Hz | quasi-fw | 0.051896% | 0.145825% | 0.902 | 0.998 |
| 0.1 Hz | quasi-tem | 0.051896% | 0.145736% | 0.900 | 1.006 |
| 10 kHz | quasi-fw | 0.003943% | 0.107035% | 0.889 | 1.005 |
| 10 kHz | quasi-tem | 0.003943% | 0.107303% | 1.384 | 1.028 |
| 1 MHz | quasi-fw | 0.003137% | 0.094159% | 1.070 | 1.006 |
| 1 MHz | quasi-tem | 0.003137% | 0.092315% | 1.191 | 1.017 |

At 0.1 Hz the quasi-fw self-R error against the matched Gamma=0 analytical
calculation falls from 4.5636% for the baseline to at most 0.06210% for normal
grading; combined-fine gives at most 0.01020%. This supports the diagnosed
conductor-geometry contribution. It does not prove every admittance component
accurate, eliminate every sign crossing or establish convergence between these
three frequency samples. Native process costs exclude Julia compilation and
mesh-validation work. They are not warmed Julia performance claims.

The public plotting API generated R/G/B figures for both physics choices from
the normal retained matrices with `clip=false`; see
`native-three-wire/plot-inputs/mixed-air_1-normal/<physics>/plots/`.
The quasi-fw G figure was inspected: it retains signed values down to about
1e-21 S/m, and the buried self and mutual curves retain their visible analytical
differences. The source tables contain every complex entry without clipping.

Production work now resolves optional conductor controls in `options.jl` and
`model.jl`, retains region-surface ownership in `geometry.jl`, applies native
transfinite arc counts and restricted BoundaryLayer/bulk fields in `mesh.jl`,
and emits editable equivalent requests through `export.jl`. The cache key
includes the controls and resolved conductor material coefficients. Disks and
annuli are implemented; sector integration and the three cable-fixture runs
remain pending. The two-wire manual caller exposes the controls and passes
them to its single detached export, preserving its physical inputs and GLMakie.

An integration check exposed a native lifecycle detail: removing a Gmsh field
does not remove its ID from the active boundary-layer list. Reusing that ID
at low frequency initially left unnecessary layers active. This is visible in
[Gmsh's field removal](https://raw.githubusercontent.com/live-clones/gmsh/master/src/mesh/Field.cpp)
and [activation list](https://raw.githubusercontent.com/live-clones/gmsh/master/src/mesh/Field.h).
Inactive prescribed layers now exclude every surface, in Julia and native
exports. A high/low/high-frequency check matches fresh low-frequency conductor
meshes (314/310 triangles in the two regions, previously 896/896 after reuse).
No failed scientific comparison triggers this operation; it implements the
requested frequency-dependent construction directly.

Focused software tests cover the 3 mm strand, dissimilar materials, optional
wall divisions and uniform growth, frequency reuse, native edits, conforming
paths, caller ownership and option grammar. Tests of exported structured PML
counts and native size caps replace incidental equality of unstructured node
counts: repeated equal targets can differ by one interior node. No scientific
accuracy threshold was relaxed. Scientific acceptance remains in the harness.

`run_bare_wire_mvp.jl` performs public `compute` integration. Its first selection
repeats the six normal native cases and adds two same-mesh detached comparisons.
Only after these pass does the serial launcher run the existing five three-wire
placements at the fixture's nine frequencies, both physics, then the two-wire
manual caller's current two radii and ten-frequency grid, both physics. It
retains raw matrices, matching analytical tables, costs and native run paths;
`plot_bare_wire_mvp.jl` renders the saved results through the public API.
The total is 138 frequency/formulation solves and 374 source columns. This is
production integration, not a repetition of the isolated conductor grids.

The launcher is
`.linecablemodels/fem/conductor-mesh-qualification/run-bare-wire-mvp.sh`.
Its harness-only executable wrapper streams actual GetDP output and fixes
MUMPS AMD ordering to match the earlier native control. Computation uses one
frequency worker and one solver thread. Completed scans and valid native
column checkpoints resume after interruptions. Monitor with:

```sh
tail -n 60 -F .linecablemodels/fem/conductor-mesh-qualification/bare-wire-mvp.log
```

### Public MVP completion and plot recovery — 2026-09-29

All 138 solves / 374 source columns completed. Every public scan has a completed
run state and finite stored matrices. Production Z/Y entries differ from the
retained native normal controls by at most 1.075e-8 relatively. The two detached
same-mesh checks give identical Z and a maximum relative Y difference of
6.27e-16. Across this selection, maximum self-R disagreement with the matched
Gamma=0 analytical calculation is 0.9292%; summed public-compute wall time is
2.562 hours, including each call's orchestration but excluding detached parity
and plotting. This is not a first-use versus warmed speedup claim.

Conductance sign disagreements remain in the raw spectra. The files
`bare-wire-mvp/comparison-summary.csv` and `G-sign-disagreements.csv` retain
per-case costs, self-R errors and individual sign-disagreement magnitudes.
For example, the all-earth quasi-fw 1-to-3 entry at 10 MHz is about +3.028e-6 S/m
against -1.800e-16 S/m analytically. Successful execution and native parity do
not establish convergence of that small mutual term. Do not clip or dismiss
these differences while reporting the conductor self-R improvement.

The final plotting pass initially revisited SVG files already saved after the
initial checks; the public export API correctly refused overwriting. The manual
plotter now skips complete six-file plot sets and generates only missing files
in partial sets. Plot-only recovery completed without running FEM: all 48
R/G/B figures exist as SVG and PNG, and a second invocation skipped all 16
completed sets. The 92 previously retained result/state/cost/plot files were
verified unchanged by SHA-256 and modification time. The recovered mixed-case
conductance figure was inspected. Sector integration and the three remaining
cable fixtures are still outstanding in the overall plan.


## Sector integration and the three cable fixtures, 2026-09-29

The production geometry owner now retains a 0.1-scale core and native triangular
transfinite strips for convex Sector conductors. Exact owned circle radii and
straight-side lengths prescribe tangential counts; the normal count follows the
longest spoke, material skin depth and user growth. Default sector resolution is
192 angular divisions; disks and annuli retain 96. Internal cuts stay inside the
same terminal/material and are excluded from its contour. Detached GEO exports
recompute the same counts after native material/frequency edits. The normalized
GetDP formulation, BCs and voltage extraction are unchanged.

A reproduced same-endpoint arc failure is fixed by skipping adjacent merged
point tags. The original sector fixture also had overlapping PVC sleeves inside
identical PVC bedding; it now represents that union as one fill. The screen's
identical PE layers similarly use one fill, removing artificial tangent seams
that broke native boundary layers. Metal shapes, materials, locations, terminal
ordering and outer cable dimensions are preserved. These are fixture
representation corrections, not an engine boolean-union fallback or a claim of
arbitrary contact support. See `fixtures.md` for the physical definitions.

The conductor mesh suite passes 100 checks, including 56 sector checks and two
checks for the merged arc. The wider mesh/export selection passes 285 checks.
Its ONELAB frequency scan initially could not bind a Unix socket inside the
sandbox. Outside the sandbox it passes 42 checks and fails two exact TSV-text
comparisons; the unchanged pre-sector snapshot reproduces the same two failures.
Observed changes are in floating-point solve digits. That pre-existing assertion
is retained, not silently relaxed. Same-mesh native numerical parity tests pass.

Geometry/mesh controls at 1 MHz and eight PML layers passed for all three cable
fixtures, with 94,708 / 22,384 / 89,530 total nodes for screen / tube / sector.
Their four-case native exports and material/frequency edits also passed; these
are integration checks, not full-exterior solve costs or scientific accuracy
claims. Evidence and logs are under
`.linecablemodels/fem/conductor-mesh-qualification/sector-integration-checks/`.

`run_cable_fixtures.jl` implements exactly the prescribed 36-case, 108-column
serial selection. Normal sweeps cover 0.1/50/10k/1M Hz in both formulations;
refined controls and same-mesh detached comparisons are at 1 MHz. Detached runs
also retain field maps. Public Z/Y, primitive Z/P, exact area/conductivity DC
references, run paths and costs are saved. `summarize_cable_fixtures.py` writes
componentwise refinement and DC-loop comparisons without clipping or convergence
claims. `plot_cable_fixtures.jl` reads those matrices through the public plotting
API; it runs no solver and writes only missing SVG/PNG files.

The batch launcher is
`.linecablemodels/fem/conductor-mesh-qualification/run-cable-fixtures.sh`.
The live log is `cable-fixtures.log` in the same directory. Completed public scans
are reused; interrupted scans use the engine's existing checkpoints. Source
identity is fixed by the harness, and one process lock prevents duplicate
launches. No production validation, fallback or self-refinement was added.
