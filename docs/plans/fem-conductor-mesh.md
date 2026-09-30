# Conductor geometry and skin-depth mesh implementation plan

Current status, 2026-09-29: the conductor sizing feature is implemented end to
end in the existing owners, with optional controls and matching detached export.
Its bounded delivery assessment is complete. The
bare-wire selection and all 36 cable-fixture cases / 108 source columns are
complete. All 18 cable R/G/B figures now exist as SVG and PNG after correcting
the manual caller's incomplete observation pair; no numerical rerun was needed.

Quasi-fw is the primary application target. Keep its coupled results explicit;
quasi-TEM and electrostatic controls are supporting comparisons, not substitutes
for quasi-fw qualification, particularly for small real admittance entries.

The authoritative current findings and limitations are in the
[delivery assessment](../../test/manual/calculations/fem_conductor_mesh/delivery-assessment.md).
Maximum cable self-R refinement changes are below 0.087%, but screen and sector
Y changes are 11.0562% and 3.32555%. Exact native Z/P parity does not close those
questions. A small electrostatic control on the saved screen meshes reproduces
its full quasi-TEM B23 and isolates deficient internal-dielectric discretization.
The four coupled screen controls completed before implementation. The correction
is now in the existing geometry/mesh/export owners: local foil edge sizes plus
a native restricted Extend field in the insulation. Full shunt refinement change
is at most 0.296%, and net-zero loop R change at most 0.113%; all 125 focused
integration checks pass. All eight final 1 MHz public cases completed:
screen/tube normal, sector normal/refined, both formulations. Delivered screen
Y reproduces the qualified prototype within 0.003100%; delivered sector Y
changes by at most 0.067236% under refinement, with loop R changing by at most
0.069526%. The tube comparison and recorded costs are in the assessment.
Do not repeat the completed historical grids.

The conductor mesh delivery scope is complete. Keep all qualification and
scientific acceptance outside production. Bare-wire G sign disagreements remain
separately open; no mesh refinement result establishes their resolution. The following
sequence and binding ownership requirements remain the plan; older qualification
matrices below record the original protocol, not additional pending campaigns.

## Delivery sequence: implementation is the outcome

The outcome is a usable conductor mesh feature in the existing FEM backend,
available to the manual runners and detached exports. Fixture creation is
complete for the selected shapes; it is not another open-ended research stage.
The sequence below takes precedence over the original expanded qualification
matrix later in this document.

1. **Close the remaining bare-wire viability question.** Apply the already
   qualified round-conductor construction to a fresh native export of the
   existing three-wire fixture. Check the current conforming voltage-path
   integration and compare baseline/graded/refined results and cost with fixed
   exterior settings. Keep this experiment in the manual harness. Reuse the
   completed isolated conductor evidence; do not restart geometry, screen or
   sector grids. If this fails, address the observed defect or cost bottleneck
   in that case before implementing the recipe.
2. **Implement the bare-wire feature.** Resolve conductor skin depth from each
   material and frequency in `model.jl`; apply local boundary resolution and
   normal grading through `geometry.jl`/`mesh.jl`. The selected round recipe is
   96 full-circle segments at the initial geometry target, first spacing at
   most delta/6 and growth sqrt(1.25), with bounded graded depth and bulk fill.
   Expose the consumed controls, preserve shared field edges, and update mesh
   identity/reuse. No equation or voltage-extraction redesign is included.
3. **Deliver both execution routes.** Emit the same constraints and editable
   material/frequency-dependent sizes through `export.jl`. Verify the actual
   feature with the existing three-wire placements and two-wire runner, both
   physics choices, raw matrices and R/G/B plots. Report accuracy and cost;
   do not promise removal of every sign crossing or analytical discrepancy.
   This is the first runnable MVP, before the remaining cable-shape extension.
4. **Extend to the three defined cable fixtures.** Reuse the same material sizing
   for screen strands and both faces of tubular walls; integrate the qualified
   sector construction at its geometry owner. Diagnose and correct the sector
   CAD arc failure without changing its physical dimensions. Use the bounded
   checks in `fixtures.md` before shipping these additional constructions.
   Do not turn the three fixtures into another Cartesian product campaign.

Completion means public optional controls, runnable native and detached routes,
documented defaults/limits, and retained end-to-end evidence for the delivered
scope. Qualification remains outside production; production neither rejects
scientific error estimates nor retries, falls back or self-refines. Long solves
remain serial with the agreed log-and-return protocol. No new batch was launched
while defining this sequence.

## First delivery and current integration contract

The existing `CurrentScenarios.three_bare_wires_problem` is the first fixture:
three distinct terminals, all-air/all-earth placements and each of the three
mixed placements. Preserve every source/receiver coefficient, including each
self term. Three source excitations produce the full 3 by 3 matrix; the nine
entries are not nine independent excitations. Reuse its physical inputs and
frequency grid rather than inventing another fixture. Cover both `quasi_fw`
and `quasi_tem` explicitly; the existing baseline caller currently selects only
`quasi_fw`, so its output alone cannot establish both-physics coverage.
The existing two-wire manual caller is the next end-to-end target. Preserve its
user-selected study settings unless a qualification run records an explicit
override. Run the qualification serially.

This scope takes precedence over the expanded shape matrix later in this plan.
Screens, tubular shells and sectors followed the bare-wire delivery. Their
bounded fixture selection is now complete; composed/material-contact cases
remain outside the qualified scope. Reuse
completed isolated round-wire geometry, DC and skin-depth evidence. Complete the
remaining applicable native-integration and cost checks on this fixture before
implementing the round-wire feature. Do not claim general-shape qualification
from that implementation.

The user's subsequent request for a minimal remaining fixture set is captured in
[three passive cable constructors](../../test/manual/calculations/fem_conductor_mesh/fixtures.md):
the Tutorial 2 wire screen plus thin Al foil, the Tutorial 3 thick Pb sheath,
and the linked three-sector cable plus concentric neutral. These replace a new
combinatorial shape campaign as the next integration selection. Their defined
frequency/refinement/parity checks are complete; see the current assessment. Construction and
geometry/mesh status is recorded with the fixtures. The zero-length arc defect
is corrected; sector grading and native control edits pass their focused checks.
The bare-wire MVP was delivered first.

The live owners, verified on 2026-09-28, are:

- `geometry.jl` creates `LCM/voltage_path/NNNN` curves and reference points with
  the model. Paths are vertical to the lowest terminal CAD vertex. Air receivers
  reference the local earth surface; buried receivers reference the outer bottom
  PML boundary. Metal intervals are excluded. Physical-domain paths are embedded
  in their host surfaces; buried PML paths share transfinite-block boundaries.
- `mesh.jl` requires each path element to be an edge of the electric field mesh.
  Preserve that connectivity, physical membership and PML construction when
  adding conductor constraints. A detached line is an invalid field trace,
  not another integration approximation to qualify or recover automatically.
- `getdp/quasi-full.pro` evaluates the solved `BF_Edge` trace directly, using
  native `Integral` postquantities and the four-point line rule `I2`. For source
  current `I_j`, it writes
  `P_ij = (v_i - v_ref,i + j*omega*integral(bt_y dy))/I_j`.
  The upward direction determines the circulation sign, independently of CAD
  edge numbering. No additional PML stretch multiplies this pulled-back trace.
- `getdp/quasi-tem.pro` writes the referenced grouped scalar potential divided
  by `UnitTransverseSource`. It has no vector-potential line-integral term.
- Julia `compute.jl` writes ordinary inputs, runs GetDP, and parses the native
  normalized Z/P columns. `results.jl` applies the requested reductions and
  computes `Y = P^-1`; there is no extra `j*omega` factor. Detached exports use
  the same field formulations and native `line-parameters.pro` matrix algebra.

There is no separate measurement-preparation stage, Python detached driver,
triangle clipping, generated path weights, stored-field lookup, contour average
or independent measurement mesh in this contract. Do not restore any of them.
Use public `compute` and fresh `export_data(:onelab, ...)` bundles with the
maintained Gmsh/GetDP entrypoints. `qualify_line_pair.py` and its preparation
profilers describe captured historical bundles, not the current MVP route.

The live [FEM documentation](../src/fem.md) describes the conforming edge-trace
contract. The earlier native-cleanup validation notes describe an intermediate
nonconforming stored-field implementation. Preserve those recorded results, but
do not present its 428 passing checks, fixed-volume line-refinement results or
the older prepared-path full-line comparisons as validation of today's trace
implementation. Conforming path refinement changes the field mesh; independent
nonconforming-line refinement is no longer an MVP convergence variable. Use the
current maintained trace/extraction tests for integration preservation, and
record fresh source identity with the fixture's end-to-end results.

## Binding boundary: qualification belongs to the test harness

**The engine runs the requested case. The user owns its scientific soundness.**
This rule governs every sizing rule, tolerance, check and gate in this plan.

- Complete the viability and scientific qualification experiments applicable to
  the delivery scope above before implementing that production feature.
  Experimental Gmsh scripts, exported-case
  variants, reference solutions, convergence studies and cost comparisons belong
  under `test/manual/` or in the test harness. Production code must not import,
  include, package or invoke that experimental/qualification machinery.
- Production implements the qualified mesh construction with optional,
  user-configurable refinement parameters. It resolves the prescribed sizes,
  generates or reuses the requested mesh, solves the case and returns the computed
  results. Parameters are meshing requests, not guarantees of scientific accuracy.
- Cheap sanity checks may emit warnings, for example when already available
  geometry/mesh statistics show a requested tolerance was missed. Such warnings
  are informational: they must not reject a result, mark a completed computation
  as scientifically failed, suppress results or change the requested calculation.
- No production reference comparisons, convergence campaigns, scientific
  acceptance/rejection, fallback mesh recipes, automatic tolerance tightening,
  remesh-and-resolve loops or self-refinement in pursuit of exactness. A missed
  tolerance must not trigger another mesh, solve or alternate formulation.
- Input-contract checks and genuine Gmsh/GetDP execution failures remain ordinary
  errors. Missing required entities, invalid connectivity that prevents assembly,
  or a failed solve are execution problems; disagreement with a reference, a
  coarse approximation or a missed scientific accuracy target is not such a problem.
- Qualification thresholds and cost budgets below are test-harness criteria for
  development decisions on controlled fixtures. They are never engine acceptance
  gates. After implementation, regression tests verify the faithful integration
  of the already qualified method; implementation is not the experiment used to
  discover whether the method is viable.

Deterministic geometry/material/frequency sizing from the user's controls is part
of mesh construction. Evaluating the result and changing those controls to seek a
better answer is a user or test-harness action. Ordinary Gmsh meshing algorithms
are not a license to add an engine-level scientific adaptation loop.

## Outcome and evidence

Resolve conductor boundaries accurately and grade the conductor interior using
its own material and frequency. Support solid wires, explicitly modelled wire
screens, tubular shells, sectors and the existing composed shapes. Preserve
headless Julia computation and the editable, detached ONELAB workflow.

The immediate verified defect is geometric. The retained mixed-pair mesh uses
12 straight segments around a circular conductor of radius 0.0425 m. Its area
deficit is 4.507034%; the resulting conductor DC resistance bias is 4.719755%.
The predicted absolute resistance error agrees with the observed self-resistance
error to 0.17%. See [the mixed-case diagnosis](../../test/manual/calculations/fem_mixed_comparison.md).

A mesh-only experiment with `Mesh.MinCircleNodes=96` reduced the predicted DC
bias to 0.0714304%. Total nodes changed from 215978 to 216247 (+0.12455%), with
identical PML element count. Evidence is retained locally in
`.linecablemodels/fem/conductor-geometry-evidence/`. This establishes a useful
geometry candidate, not broadband solver accuracy or a measured speed benefit.
Enabling global curvature sizing was ineffective with the existing size floor;
removing that floor made the alternate experiment expensive enough to stop.

Current `_fem_region_mesh_size` is geometry based and frequency independent.
Existing exterior/soil grading does not resolve conductor skin depth. Julia mesh generation and detached native exports use the same owned geometry
construction and Gmsh mesh settings.

## Fixed scope

- Retain the current physical equations, materials, source and terminal
  conventions, voltage extraction, boundary conditions, exterior and PML rules.
- Retain the current first-order triangles and contour segments. Preserve native
  measurement paths as shared field edges. No quadrilateral recombination or
  higher-order element migration in this feature.
- Preserve shared interfaces, terminal ordering, material coverage and the
  existing rules for merging fully occupied formations. Do not homogenize
  explicit screen wires or merge distinct electrical terminals to save elements.
- Separate geometry error, conductor discretization error and exterior error in
  validation. Hold the exterior fixed during conductor convergence experiments.
- Compare identical primitive matrices and explicitly identical reduction
  settings. The earlier analytical ideal-transposition mismatch is documented;
  no test may reproduce that comparison mistake. Do not require symmetric Y or
  remove legitimate zero crossings. No clipping, sign repair or area correction.

The separate quasi-TEM analytical discrepancy is not attributed to conductor
meshing without evidence. Both formulations must preserve their established
behavior and converge under mesh refinement; this feature is not a promise to
eliminate every remaining model discrepancy.

## Meshing decisions

### 1. Geometry fidelity is independent of coarse bulk sizing

Offer an optional conductor geometry target, initially `1e-3` for the qualification
fixtures. In the test harness, require relative area error at most this value for
**each conducting region**, not just the sum over a terminal. Also check boundary
approximation and preservation of positive gaps and wall thickness: an area match
alone cannot qualify a shape. In production, this parameter prescribes geometry
construction; missing the target may produce a cheap warning, never a rejection
or an automatic refinement.

For circles, use the exact inscribed-polygon area relation
`A_mesh/A_exact = N*sin(2*pi/N)/(2*pi)` to select a boundary count. Start with
96 segments for the default target. Use native curve discretization constraints
in the production geometry owner, distributed over the actual arcs and their
shared junctions. The successful global `MinCircleNodes` probe is evidence, not
a reason to impose 96 segments on every unrelated curve in the domain.

For other smooth boundaries, prescribe local curvature/chord-error control and
qualify it with area checks in the harness. A starting chord-error target is
`1e-3` of local curvature radius, combined before meshing with known thin-feature
and gap dimensions. For concentric shells, match angular subdivisions on the
inner and outer faces; verify wall thickness and annular area in the harness.
Applying unrelated polygon approximations to two nearly equal areas is
unacceptable as a design choice. These checks must not become a production
measure-adjust-remesh loop.

Native arcs and already polygonized curves need different treatment. Subdividing
a straight chord cannot recover an ellipse offset or curved strip that was
approximated too coarsely during geometry construction. Refine that approximation
at its existing geometry owner. Keep true polygon corners and straight sides;
do not silently round sectors or manufacture fillets.

The 3 mm diameter case is mandatory. At 96 segments its boundary chord is about
98 micrometres, and its relative polygon-area error is the same as the larger
wire. The count itself does not grow as radius decreases. Explicitly modelling
many strands does increase total cost, which must be measured.

### 2. Resolve decay normal to the conductor surface

For each material and frequency, use the constitutive values already evaluated
by the model. With the existing harmonic convention, let
`q = sqrt(im*omega*mu*kappa)` and select the decaying branch, `real(q) >= 0`.
Use `delta = 1/real(q)` when there is finite attenuation; for good conductors this
reduces to `sqrt(2/(omega*mu*sigma))`. Do not reuse soil skin depth or add another
material law. Retain a local phase-resolution bound if a supported conducting
material is outside the good-conductor limit. Verify the convention against the
existing constitutive implementation before emitting any mesh formula.

The original delta/3 candidate failed the completed isolated spectrum. The
selected round-conductor recipe for the remaining native integration check is:

- first normal element size no larger than `delta/6`;
- normal growth ratio no larger than `sqrt(1.25)`;
- graded refinement extending to approximately `5*delta`, where the conductor
  is thick enough, followed by a smooth transition to its existing bulk cap;
- at least four elements through a thin wall when opposing skin regions overlap.

These are candidates to qualify, not universal error bounds. The full conductor
remains in the PDE: `5*delta` is a refinement extent, not a truncation boundary.
Tangential spacing is controlled separately by geometry, proximity effects and
corners. Reducing the normal step must not force an isotropic fine grid throughout
the conductor or surrounding air/earth.

When delta exceeds the conductor thickness/radius, resolve the finite section
directly; do not construct a nominal skin region outside it. Copper at 20 C with
resistivity `1.7241e-8 ohm m` has delta about 9.35 mm at 50 Hz and 66.1 micrometres
at 1 MHz. Thus a 3 mm wire needs no thin skin region at 50 Hz, but does at 1 MHz.

Apply grading to real material interfaces, using each side's own material.
Disconnected wires belonging to one terminal remain separate geometric objects.
Artificial cuts inside the same material must not become artificial skin walls.
At material contacts, preserve the existing conforming interface and electrical
assignment; never duplicate coincident faces to obtain independent layers.

### 3. Use Gmsh's boundary-layer and transfinite capabilities

The primary implementation uses Gmsh `BoundaryLayer` fields, activated through
`setAsBoundaryLayer`, with conductor-side surface exclusions and triangular
layers (`Quads=0`). Existing scalar background fields continue to control bulk
and exterior sizes. Boundary-layer activation is separate from a scalar `Min`
background field. Corner fans and explicit transfinite curves/patches are the
native tools for constrained local topology. These are documented in the
[Gmsh manual](https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-mesh-size-fields) and its
[versioned boundary-layer example](https://raw.githubusercontent.com/live-clones/gmsh/gmsh_4_15_2/examples/api/naca_boundary_layer_2d.py).

**Start with a small capability experiment in the test harness, before any
production feature implementation.** On the current built-in GEO kernel,
demonstrate:

| Geometry | Required observation |
|---|---|
| Two disconnected conductors with different materials | Both receive their own normal step and depth in one mesh; no field overrides the other. |
| Annulus, with wall/delta ratios 0.2, 2 and 20 | Refinement acts from both faces; the narrow-wall mesh meets in the interior without overlap, missing layers or slivers. |
| Sector and bent strip, including a sharp/reentrant corner | Native fans/local patches preserve the specified corner and form valid triangles. |
| Adjacent material regions and explicit screen wires | Shared curves remain conforming; no refinement leaks into unrelated faces or changes material/terminal membership. |
| Any one of these inside the current exterior | PML transfinite constraints and prescribed remote sizes survive unchanged. |

Inspect actual element heights, adjacency, coverage and Jacobians; successful
Gmsh return status or a field declaration is insufficient. Check multiple active
boundary-layer fields explicitly. The upstream field manager stores multiple
IDs, but that alone does not prove the required material-side behavior.
On the first strongly graded solves, also check solver residuals and terminal
current normalization against existing tolerances. High aspect ratio alone is
not a defect in an intentional anisotropic layer; invalid elements or degraded
numerical conditioning cannot be dismissed as the price of fewer nodes.

During qualification, compare boundary-layer and conforming transfinite strip
constructions for thin regular shells/strips. Establish the prescribed
construction for each supported thickness/topology before implementation. The
production selection follows that explicit input-based rule once; it must not
try one construction, assess its accuracy and fall back to the other. Test that
internal partitions disappear from terminal contours.
For a corner, qualify the local fan/patch using integrated losses and impedance;
convergence of a singular pointwise peak current is not an acceptance condition.

If Gmsh cannot provide the required general shape coverage with these bounded
constructions, record the exact failing geometry and revise the experiment/design
before production implementation. Do not hide failure with global isotropic
refinement, a CAD kernel rewrite, a custom mesher or a runtime fallback.

## Integration and ownership

Expose optional refinement controls through the existing FEM computation options.
Users can override the qualified defaults without editing code. The proposed
values below are experimental candidates, not scientific acceptance thresholds:

| Control | Initial candidate | Meaning |
|---|---:|---|
| `conductor_geometry_tolerance` | `1e-3` | Requested relative geometry construction target; it does not guarantee or enforce achieved area/impedance accuracy. |
| `conductor_skin_depth_elements` | `6` | Reciprocal first normal step: `h_normal <= delta/value`; it does not claim this many uniform cells through every delta. |
| `conductor_mesh_growth` | `sqrt(1.25)` | Requested ratio between successive normal element sizes. |
| `conductor_skin_depths` | `5` | Requested graded extent in skin depths, limited by the known section geometry. |
| `conductor_thickness_elements` | `4` | Requested element count through a thin wall when opposing graded regions meet. |

Normalize and check parameter types/ranges once in `src/engine/options.jl`.
Values outside the qualification fixtures' refinement levels are not grounds for
rejection: valid user-selected coarser/finer settings are honored. Choose the
defaults from the pre-implementation evidence, document their tested range and
limitations, and keep topology rules in the geometry owner. No selectable
mesh-backend framework or new formulation option is needed.

Existing `mesh_size_factor` and exported `MeshScale` scale their bulk/proximity
targets before applying the user-requested conductor size constraints. Coarsening
those targets must not undo the explicit geometry or skin-depth controls. Document this
tightened precedence and test it with the manual runner's factor 3 and a larger
factor. Avoid leaving a global minimum size or post-applied global scale that
silently defeats the local constraints.

| Existing owner | Required change |
|---|---|
| `model.jl` | Resolve material/frequency conductor sizes and local thickness information from existing region/material plans. Keep ordinary geometry values out of types. |
| `geometry.jl` | Own accurate curve sampling, material-side boundary membership and any necessary conforming local partitions. Retain exact mesh-constraint metadata alongside existing PML metadata. |
| `mesh.jl` | Apply those constraints and native fields in the existing sequence; preserve exterior grading and global Gmsh state. Reapply/reset conductor constraints when changing frequency. |
| Mesh fingerprint/cache | Include the new controls, resolved constraints and meshing revision so old coarse meshes cannot be reused as corrected meshes. Preserve per-frequency compatible reuse. |
| `export.jl` and native export assets | Serialize the same constraints, field activation, side exclusions and precision; extend the existing allowlist and global-state restoration where necessary. |
| Existing tests/docs | Own all qualification criteria, references and comparison logic. Keep experimental recipes and expensive numerical studies under `test/manual/`, outside shipped feature code. |

Prefer fields added to the existing resolved plans and geometry metadata over a
second planning pipeline. Reuse region-to-entity information from construction;
do not infer conductor identity from coordinate tolerances or physical tag order.
Preserve checks for structurally missing/duplicated material entities; do not turn
the new geometry-accuracy targets into runtime material-coverage rejection rules.
Any new cheap tolerance diagnostic reports the observed value and requested target
through an ordinary warning and leaves the calculation and finite results intact.

The detached exporter must preserve editable material/frequency semantics.
Changing a conductor's native conductivity/permeability and remeshing must update
its local skin resolution. A frozen export-time delta is insufficient. Emit the
corresponding visible native derived expressions from the existing export owner,
using the same constitutive inputs and equation as Julia. The detached runtime
uses native Gmsh/GetDP exclusively. Check the native parser expressions against
Julia on the qualification cases.

For externally supplied meshes, preserve the existing no-remeshing contract.
Report that the new generation controls do not certify an externally supplied
mesh; do not silently replace it. Keep numerical convergence checks in tests and
manual qualification, rather than introducing production refinement/retry loops.

## Pre-implementation qualification: test harness only

Complete the delivery's applicable A, B and C checks using experimental
meshes/native case files and current solver entrypoints before changing
production feature code. The broader shape matrix below is retained for later
generalization; it does not override the bare-wire MVP priority. Keep all reference
formulae, error thresholds, refinement sweeps and cost decisions in the harness.
Experiments may compare recipes and rerun finer meshes explicitly; none of that
driver logic is to be transferred into the engine. Only the selected construction
and user-configurable sizing rules are implemented after qualification.

Use three mesh levels: candidate, finer geometry/tangential resolution, and finer
normal resolution. Vary the latter two independently before the final combined
refinement. Record actual element sizes and geometry areas at each level.

### A. Geometry and isolated conductor controls

1. Circular conductors with radii 1, 1.5, 10, 42.5 and 85 mm. Check region area and
   low-frequency conductor-only resistance against `rho/A_exact`.
2. Solid and tubular round conductors with independent analytical internal
   impedance references. Include weak, transitional and strong skin effect, and
   both thin and thick walls. Do not use total earth-return R to hide an internal
   conductor error.
3. Sector, bent-strip and explicit screen fixtures: compare integrated losses
   and terminal Z against independently refined meshes. Include close wires,
   dissimilar materials, material contact, and supported shared-terminal cases.
4. Screen scaling at 24, 48 and 96 explicitly modelled 3 mm diameter strands,
   holding local strand geometry, material and minimum gap fixed for the cost
   scaling study. Also retain a fixed-envelope physical screen comparison.

Initial release acceptance targets are per-region area error <=0.1%, isolated
DC resistance error <=0.2%, and isolated AC internal resistance/complex impedance
error <=1% against the independent round-conductor references. Require <=0.25%
change in integrated resistance under the final refinement on those controls.
These are fixture-specific harness assertions; they do not relax existing stricter
tests or claim the same percentage bound for complete line matrices. If the initial
grading fails, improve the experimental recipe before implementation instead of
loosening the targets. Do not copy these assertions into the feature.

### B. Full formulations and spectrum

Run both `quasi_fw` and `quasi_tem`, each with two conductors in air, two in earth,
and one in each layer. Cover soil resistivities 0.1, 1, 100 and 1000 ohm m.
Start with the representative 42.5 mm radius pair and the existing geometry and
exterior preset, with explicit equal reduction settings on both references.

Use 0.1, 1, 10, 100, 1e3, 1e4, 1e5 and 1e6 Hz for the initial coverage: 192
frequency cases across the 24 combinations, retained/reusable rather than rerun
as one opaque job. Add intermediate frequencies near skin transitions and every
suspected sign crossing. Qualify other radii/shapes first with A, then use low,
transition and high-frequency placement controls instead of multiplying the
entire sweep by every geometry.

Retain raw complex Z, P and Y, loss/normalization quantities, mesh statistics,
solver logs and selected current/field maps. Produce R, G and B versus frequency
plots through the plotting API with its explicit unclipped setting. Use raw
numeric tables to evaluate signs and errors; a signed-axis display is not a test.

For every crossing, compare converged values on both sides and the crossing
bracket. Where a component is smaller than the observed mesh-convergence error,
report its sign as unresolved, rather than force agreement or divide by nearly
zero. Require no new resolved-sign disagreement with the correctly matched
quasi-fw reference. Report existing quasi-TEM reference differences separately;
require convergence and preservation, not an unsupported model-equivalence claim.
An unexplained material regression in any placement or formulation leaves the
experimental recipe unqualified for implementation. This is a development decision
made from harness evidence, never an engine decision about a user's result.

### C. Computational feasibility and native parity

Measure one-worker runs on the same machine, with fixed exterior and solver
settings. Report first-use Julia compilation separately from warmed mesh/solve
time, peak resident memory, DOFs and element counts by conductor/bulk/PML.

- Geometry-only correction on the representative pair: target <=5% additional
  total nodes; the existing mesh-only result is comfortably below this.
- Complete grading on the pair: target <=1.5 times warmed wall time and peak
  memory of the current preset, while meeting the accuracy gates. This is a
  qualification budget, not a performance claim already established.
- Screens and thin shells: compare against an accuracy-matched locally isotropic
  reference, not only an inaccurate coarse mesh. Require no increase in peak
  memory or warmed runtime against that reference; seek substantial savings in
  the strong-skin cases. Also report absolute cost and overhead against today's
  mesh so the tradeoff is visible.
- Record the 24/48/96-strand growth in conductor DOFs and absolute solve cost.
  Approximately linear local mesh growth is the target; direct-solver time need
  not be linear. Never infer constant total cost from the two-wire experiment.
- Run the qualification serially first. Use existing caches and completed-case
  artifacts; do not launch simultaneous dense sweeps or introduce a second cache.

A failed cost target leaves that experimental preset unqualified for implementation.
Investigate where nodes and factorization cost accumulate; change the experimental
mesh design without weakening the accuracy requirements. Do not quietly coarsen a
strand, drop a shape from coverage or impose an accuracy-destroying node cap. If
accuracy and cost cannot both be met, report that limit with the failing case and
measurements. Production does not enforce these fixture budgets on user cases.

Before implementation, use experimental native files to check Gmsh API/native
parser equivalence for the proposed constraints and sizing expressions. Check
material memberships, local sizes and resulting Z/Y on representative circle,
thin-shell, sector and screen cases. Exercise changing frequency and native
material values before remeshing, using the existing detached workflow. Retain
same-mesh solver parity checks separately from remeshing comparisons, where node
numbering need not be identical. After feature integration, repeat the affected
parity, material-edit and relocated/no-Julia checks against the actual exporter;
these verify integration of the qualified method, not its initial viability.

## Long-running qualification: visible terminal, manual resumption

**Run the qualification commands in a visible terminal with solver output, then
end the agent turn. The user watches the terminal and reports completion.** Keep
all necessary cases and frequencies. Waiting for Gmsh/GetDP requires no model
activity, automatic notification or separate agent.

Only the main implementation agent works on this campaign. Do not create
subagents, parallel agent runs, forked conversations, replacement `codex exec`
sessions, timers, daemons or a job-supervision framework.

1. Prepare the next qualification batch and its exact invocation. Record the
   inputs, controls and output directory using the existing study/run artifacts.
   Start with one numerical worker so terminal output remains readable.
2. Run the ordinary foreground Gmsh/GetDP command or manual study in the visible
   IDE terminal, with useful verbosity and stdout/stderr displayed. Leave the
   command attached to that terminal; do not launch an orphaned background child.
   If the agent execution pane cannot provide a persistent user-visible terminal,
   give the user the exact command to launch in the IDE's integrated terminal.
3. Keep normal result files and solver logs. A shell `tee`, with the solver's exit
   status preserved, is sufficient when displaying and saving the same output.
   Gmsh's existing verbosity controls enable terminal output. The current Julia
   GetDP launcher redirects stdout/stderr to a log, so raising its verbosity alone
   does not display that log: for qualification, run the generated solver command
   directly or show the existing log in the terminal. Do not redesign production
   logging just to conduct these experiments.
4. Confirm startup once, identify the terminal command and output directory, and
   end the turn. The user leaves the terminal open and sends `continue` or reports
   that the command has finished. No repeated agent polling, log-tail calls,
   heartbeat messages or automatic goal continuations while waiting. If a goal is
   active, honor the user's pause-while-waiting instruction through its pause control.
5. When the user resumes this same conversation, inspect the exit status, logs and
   saved outputs, then interpret the results. Reuse compatible completed records
   through the existing native resume contract. Do not rerun finished cases merely
   to restore context. If the job is still running, report that once and stop again.

The output directory and a short note of the command/current batch are the
checkpoint. No new status-file schema, process supervisor, terminal-marker
protocol or notification subsystem is required. A solver finishing means its
execution ended; scientific qualification remains a separate harness assessment.

The earlier failed detached-child experiment does not qualify or disqualify a
normal foreground terminal run. Retain its evidence, but remove detached-job
survival and automatic wake-up tests as prerequisites. For the next small Q1
case, confirm that the selected terminal displays solver output; use the same
ordinary execution workflow for the longer batches.

## Execution order and completion gates

| Stage | Deliverable and gate |
|---|---|
| Q0: capture and terminal setup | Record current versions, source diff and fixed baseline inputs/results. Prepare the ordinary visible-terminal invocation and output directory. Preserve all unrelated work; no background-supervision experiment is required. |
| Q1: capability experiments | Under `test/manual/`, establish the Gmsh boundary-layer/transfinite constructions and native expression parity on the small topology cases. No production feature changes. |
| Q2: scientific and cost qualification | Complete A, B and C using the experimental recipes: geometry/DC/AC/shape controls, both formulations and all placements across the spectrum, sign investigation and measured cost comparisons. Retain tables, plots and maps. No production feature changes. |
| Q3: freeze the qualified method | Record the selected input-based construction rules, optional controls, defaults, fixture coverage and limitations. Resolve failed qualification experiments before starting implementation. |
| I1: implement | Implement only the qualified construction and user-selected controls in the existing geometry/model/mesh owners, with reset/reuse semantics and cache invalidation. No qualification machinery or runtime scientific acceptance logic. |
| I2: integrate export | Serialize the same constraints and editable derived sizes in the detached workflow. Verify faithful native/parser/remesh/material-edit integration and relocation. |
| I3: regression and documentation | Run affected option, mesh, partition, voltage-path, cache and export tests against the qualified fixtures. Document controls and precedence, tested defaults and limitations; update the manual runner only as needed. |

Q0-Q3 precede I1-I3 for the delivery scope, reusing applicable completed evidence.
The first delivery does not wait for all deferred shapes. Integration tests
remain necessary after implementation, but
they must not substitute for pre-implementation viability and scientific
qualification. If integration exposes a new scientific uncertainty, investigate
it in the harness; do not add an engine validator or recovery policy.
Each long batch follows the stop-and-resume protocol above. Waiting for a batch
does not justify another reasoning turn or a reduced qualification matrix.

Use existing tests such as `test/extensions/fem_mesh_grading.jl`,
`test/extensions/fem_partition.jl`, `test/extensions/fem_export.jl` and
`test/unit/engine/options.jl`; introduce a focused conductor-mesh test file only
where that gives a coherent owner. Put the numerical/cost study in
`test/manual/calculations/run_fem_conductor_mesh.jl` with a companion results
document. Keep it excluded from ordinary automated test discovery.

Include focused regression coverage showing that a cheap diagnostic for a missed
target warns while the requested case still runs and returns its results, without
retry, fallback or changed refinement parameters. Review the production import/
include graph to ensure it has no dependency on qualification scripts, reference
fixtures or harness acceptance thresholds.

Completion requires evidence for every gate, including limits and residual model
differences. A 96-segment mesh screenshot, a successful solve, or agreement in
only the aerial pair is not completion of this feature.
