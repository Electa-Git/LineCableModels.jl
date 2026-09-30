# PML mesh from physical scales and a smaller physical domain

Locked and execution authorized on 2026-09-29. This implements the plan accepted
in the conversation, not an extension of the completed fixed-count campaign.
One implementation agent; Julia, Gmsh and GetDP only. Evidence and the single
append-only live log: `.linecablemodels/fem/pml-physical-mesh/`.

## Objective and preserved contracts

Derive the exterior node distribution from propagation/attenuation, the
coordinate stretch, and variation of the transformed coefficients. Assess
moving the physical/PML interface to roughly two soil skin depths, subject to
layout clearance. Reduce actual unknowns and native execution cost, not only
the physical dimensions. Target at least 30% less native time than the retained
144/144/96 prescription; this is a target, not a promised result.

The reference prescription is the completed 24-skin-depth, 144/144/96 campaign
in `fem-pml-conductance-cost.md`. Reuse its raw results. Preserve the quasi-fw
equations, conductor geometry/grading, physical material ordering, excitation,
and complete native GetDP voltage extraction. Preserve every observation pair.
The initial design keeps the existing continuous stretch; changing its profile
is a separate decision supported by modal evidence, never an unrecorded tweak.

Scientific qualification belongs solely in the manual harness. Production
constructs the requested mesh and solves once: no reference comparisons,
scientific rejection, repeated trial solves, fallback, clipping, sign correction
or self-refinement. User controls specify resolution, extent and thickness;
cheap diagnostics may report unmet prescribed resolution. Julia-managed and
detached `.geo`/`.pro` execution must share one resolved mesh prescription.

## 1. Cheap design and numerical feasibility, before feature code

Use material propagation constants, complex prescribed Gamma, interface
distance and the PML map. Resolve phase advance, decay and transformed-tensor
variation over each proposed cell. Include propagating and evanescent modes,
small normal propagation constants and the zero-Gamma limit. Mode ranges and
any tail cutoff are explicit; a bulk normal plane wave alone is insufficient.

A Julia one-dimensional finite-element harness compares the discrete PML
boundary response with the known continuous finite-layer and outgoing
responses. Separate finite-wall truncation from discretization error. Keep
this harness outside production. This is a cheap design filter, not a
certificate for the coupled layered cable problem or tiny real admittance.
Near-cutoff modes and shared air/earth side stretching are reported explicitly.

Derive one candidate distribution. Express it using a small number of native
conforming transfinite strips, with physical breakpoints, counts and geometric
progressions. Avoid manual post-mesh node repair or one CAD patch per cell.
Side/top/bottom intersections and the native reference path share nodes.
Count resulting corner and strip elements before any large solve. The first
deliverable is the modal response, resolved cell measures and predicted cost.

The old exponent (192/191)*log(1536) is a comparison distribution, not a physical
design target. All numerical accuracy controls retain explicit meaning; no
grading is fitted directly to analytical cable conductance.

## 2. At most eight exploratory coupled frequency solves

All exploratory coupled solves are serial. Keep the complete original problem
frequency definitions when preparing geometry, solving only selected points.
Reuse saved reference runs and inspect each step before advancing.

| Step | Change | Frequency solves |
|---|---|---:|
| A | New distribution; retain physical domain, PML thickness and stretch | 2 |
| B | Move PML entrance to about 2 delta; retain PML thickness | 2 |
| C | Reduce PML thickness using the modal design's attenuation requirement | 2 |
| D | Check 100 ohm m aerial and mixed three-wire cases | 2 |

A-C use the 0.085 m, 0.1 ohm m aerial case at 0.1 and
21.544346900318832 Hz. D uses the retained 100 ohm m aerial case and mixed
three-wire fixture at 1 MHz. The geometric clearance floor stays explicit.
No per-case tuning. If a step fails, remaining exploratory solves diagnose that
failure instead of advancing to a full batch. The eight-solve limit is firm;
cheap modal/mesh-only calculations do not consume coupled solves.

Mandatory scientific gate: retain the reference conductance signs. Report
signed and absolute R/X/G/B errors independently, self-R and mutual components,
and changes from the current reference prescription. Preserve existing fixture
requirements; invent no new percentage accuracy policy. Do not adopt a result
with unexplained accuracy degradation merely because it is faster.

Record actual mesh counts, DOFs, assembly/factorization time and peak memory.
Separate compilation from execution. Stop before production changes/full
qualification if modal response is poor, signs fail or savings are negligible.
Document a negative feasibility result rather than silently broadening the
search or relaxing checks.

## 3. Production implementation only after successful qualification

Implement the successful construction in the existing mesh-plan and geometry
owners. One resolved prescription drives both native managed meshing and
detached export. Keep definitions passive and introduce only controls with
demonstrated numerical meaning. Preserve caller-owned explicit mesh overrides,
cache/resume identity, native terminal paths and shared interface constraints.
Qualification helpers, mode sweeps and result gates do not enter feature code.

Check native geometry generation/parsing before coupled numerical runs. Verify
both source paths produce the prescribed nodes and conforming corners. Document
the new prescription, its controls and the measured validity limits.

## 4. One final selection for one surviving prescription

Reuse the established physical fixtures and comparison readers:

- Ordinary two-wire study: seven parameter cases, ten frequencies each (70).
- Three buried and mixed three-wire cases at three frequencies (6).
- Screen, tubular shell and sector preservation at 1 MHz (3).
- One independent detached solve and a same-mesh execution comparison, reusing
  saved worker results where possible; distinguish remeshing from solver parity.

Only this final selection may use the established four frequency workers with
one native thread. No full batch for rejected candidates. Compare raw G signs
and all component errors; no observation clipping. Preserve the complete
frequency spectrum and terminal ordering. Produce tables and conductance plots.

## Execution and references

One live log, native stdout, per-case completion markers and saved controls.
Restart from completed work. No parallel agents, Python, per-trial approval
requests or new production monitoring machinery. After launching a long batch,
report its log and checkpoint rather than burning context on frequent polling.

- Johnson, [Notes on Perfectly Matched Layers, sections 6-7](https://ocw.mit.edu/courses/18-303-linear-partial-differential-equations-analysis-and-numerics-fall-2014/0ad128a4b3d9dbb860e83a59a47b1b01_MIT18_303F14_pml.pdf): evanescent stretching, discrete reflection and grazing incidence.
- Druskin, Guddati and Hagstrom, [On generalized discrete PML optimized for propagative and evanescent waves](https://arxiv.org/abs/1210.7862): discrete boundary-response design. Its optimal-grid construction is not assumed to be directly applicable to this coupled layered formulation.

## Execution record

2026-09-29: plan locked; modal design and cost assessment started. Production
sources and the working manual preset remain unchanged at this stage.

Source-scope clarification at qualification start: the quasi-fw backend was the
first-order Gamma-to-zero reduction (`quasi-full.pro:6`), and did not accept a
finite Gamma input. Finite-Gamma modal probes are retained as theoretical
stress tests, not represented as production coverage or a requirement to add
finite-Gamma physics during this mesh task. Adoption uses the actual Gamma=0
equations. The complete finite-Gamma formulation remains outside this task.
Concurrent work subsequently added it; this campaign's coupled results still
cover Gamma=0.

2026-09-29, feasibility decision: eight exploratory native solves completed.
The scalar-interpolation candidate failed G signs. Three controls isolated
side/top distribution error and found negligible improvement from native
quadrature 12-to-13. Adding coefficient-variation density reduced that error;
the tighter 0.12 prescription passed both original sign probes. Native time
improved by 19.2% at 0.1 Hz and 32.4% at 21.54 Hz. This is sufficient to
assess the single surviving prescription in the final retained selection.

The exploratory budget was consumed by diagnosis. Therefore the final
selection retains 24 skin depths and the original thickness/stretch; it does
not silently substitute an untested 2-skin-depth domain. Smaller-domain
production support is not qualified by this campaign. Run the final selection
on the isolated in-memory prototype before implementing feature code, honoring
the user's qualification-before-implementation requirement. No full selection
is run for either rejected prescription. Actual resolved counts and native
strips are saved per case; the prototype's runtime/cache is isolated from
production so its temporary geometry override cannot contaminate normal runs.

Implementation ownership, conditional on final qualification:

- Add an optional `pml_resolution` prescription containing interpolation-cell
  density and coefficient-variation controls. Retain explicit fixed
  `pml_layers`/`pml_grading` as a distinct caller choice; never silently ignore
  explicitly supplied fixed controls when physical resolution is requested.
- Resolve normalized strip endpoints, interval counts and native progression
  ratios once in mesh planning. Geometry and detached export consume those
  same resolved strips. Cache/resume identity includes that complete layout.
- Keep only the prescribed density/strip construction in the engine. Scalar
  FE response assessment, mode-convergence checks, reference comparisons and
  all qualification decisions remain in the manual harness.
- Preserve physical-domain and conductor targets, shared interfaces and
  native voltage-path segments. A fixed one-strip prescription must retain
  its existing geometry and native count semantics.

The original smaller-domain objective remains outstanding. An asynchronous
question requests authorization for four additional serial probes (the
original B/C controls at 0.1 and 21.54435 Hz), because the eight-solve limit
was exhausted. Do not execute those extra probes without the answer; do not
claim the original smaller-domain objective is complete from the retained
24-skin-depth selection. Final qualification and implementation assessment
can continue independently.

2026-09-29, final selection and implementation: all 79 frequencies / 167
source columns completed before feature code. All 363 conductance signs matched
their references. Seven ordinary scans took 561.6 s versus 889.2 s for the
retained fixed mesh; first-call compilation was 28.7 s versus 26.4 s. This is
a 36.8% measured elapsed reduction in one batch each, not a universal speedup.

The optional `pml_resolution` construction is implemented in the existing mesh
plan and native geometry owners. It contains no trial solve or scientific
acceptance machinery. All 79 resolved layouts match the frozen qualification
prescriptions (6874 assertions). Public managed compute, same-mesh detached
`.pro` execution and explicit resume passed. The ordinary manual runner uses
the new control for both managed studies and its one detached export.

Scientific limits remain explicit. The largest two-wire G relative error is
76.9%; self-R stays at about 0.466%. One tiny screen mutual G degrades by
24.8% relative to its retained reference (9.83e-15 S/m absolute). This does not
establish an all-fixture magnitude improvement. Fixed controls and global API
defaults are retained. The new prescription is an optional cost/sign tradeoff,
not a replacement accuracy standard. Detailed components, plots and timings
are recorded in `test/manual/calculations/fem_pml_physical_mesh.md` and its
linked evidence directory. Smaller-domain qualification is still pending.

2026-09-29, user authorization: "yes proceed" approves the four prepared
smaller-domain field solves. The exploratory allowance is now twelve total:
the original eight plus B/C at 0.1 and 21.54435 Hz. Run them serially using
their frozen native exports and meshes. Inspect B before C; keep all raw
component errors and reject a reduced-domain preset if signs regress. This
extension does not authorize another tuning campaign or alter the production
scientific-validation boundary.

2026-09-29, smaller-domain assessment complete: all four authorized probes
ran serially using the prepared native bundles. Both 2δ/24δ and 2δ/2δ
controls fail all four G signs at both frequencies (16 failures / 16 entries).
Self-R errors remain below 0.043%. The thick-PML controls take 20.57/16.44 s;
the thin-PML controls take 27.01/22.12 s versus retained 29.51/24.34 s.
Saved DOFs, peak RSS, assembly/solve times and full component comparisons
are under `domain-probes/`. The thinner PML improves the wrong-sign error
but does not remove it. Neither reduced-domain control is adopted.

The bounded plan ends with the optional physical mesh prescription at the
qualified 24δ extent and a negative two-skin-depth feasibility result. The
original eight exploratory solves plus the four explicitly authorized probes
are complete; the single final retained-fixture selection and both execution
paths are complete. No further parameter campaign is implied. Accuracy limits,
including the screen mutual-G degradation and finite-Gamma coverage boundary,
remain part of the delivered result rather than being silently relaxed.

## Completion audit

Evidence root: `.linecablemodels/fem/pml-physical-mesh/`. Numerical ledgers
were reread from native matrices/comparison rows and completion markers;
`completion-audit.toml` records the totals. Completion means delivery of the
bounded mesh feature and assessment, including a rejected domain reduction;
it does not mean every component's analytical percentage error is resolved.

| Requirement | Authoritative evidence and outcome |
|---|---|
| Physical scales, propagating/evanescent/cutoff modes, static limit | `modal.jl`, `strips.jl`, `modal.csv`, `strips.csv`; finite-layer and outgoing responses separated. Finite-Gamma modal limits reported, not promoted to coupled coverage. |
| Native strips and conforming interfaces, without node repair | `geometry.jl` consumes `FEMMeshPlan.pml_strips`; 2880 current focused assertions include native exported nodes, shared interfaces, buried paths, cache identity and fixed controls. |
| Qualification before implementation | Completed prototype selection and frozen `resolved/` records; production layouts match all 79 frequencies in 6874 checks. |
| Bounded exploratory work | Initial eight controls plus four explicitly authorized domain probes; all serial. Full original frequency definitions preserved. |
| Roughly two-skin-depth domain assessment | Four native completion markers and component tables under `domain-probes/`; 16/16 G signs fail, so neither control is adopted. |
| One complete retained selection | 79 frequencies / 167 source columns; all 363 G signs retained. Includes seven ordinary scans, buried/mixed three-wire cases and screen/tube/sector fixtures. |
| Magnitudes, self-R, cost, memory and limitations visible | `components.csv`, `component-summary.csv`, `costs.csv` and per-probe `solve.toml`; screen mutual-G degradation explicitly retained. No global default or cable preset is replaced based on that degradation. |
| Meaningful computational savings | Seven-scan elapsed 889.2→561.6 s (36.8%); summed native worker time 2822.2→1659.8 s (41.2%). Compilation reported separately; one batch each. |
| Shared managed and detached execution | `detached-preservation.toml` distinguishes independent remeshing from same-mesh parity; `feature-execution/complete.toml` verifies public compute, detached native algebra and resume. |
| Public controls, ownership and numerical policy | Optional `pml_resolution` in computation/export; strip metadata and source digest participate in cache/resume identity. `pml_mesh.jl` constructs a fixed prescription, with no field solve, result gate, retry, fallback or clipping. |
| Manual caller and deliverables | Two-wire runner selects the new control in both managed and detached calls. Sixteen figures exist as SVG and PNG (32 artifacts); raw component tables and qualification notes are retained. |
| Resumability and execution protocol | One `live.log`, per-native completion markers and `resume.sh`; completed domain probes are reread without new solves. One agent, no Python. |

The remaining magnitude-accuracy and arbitrary finite-Gamma questions are
documented limits of this prescription, not silent claims of successful
qualification. The delivered manual preset remains the tested 24δ choice.
