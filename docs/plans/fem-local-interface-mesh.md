# Cable-local physical-domain mesh grading

Execution authorized 2026-09-30. This work changes surrounding-medium size
fields only. Preserve the conductor prescription, Cartesian physical domain,
PML tensors and strips, all PML nodes, and native voltage-path quadrature nodes.
Both Julia-managed runs and detached native ONELAB exports must be delivered.

## Qualification before feature implementation

The existing wave-size fields measure distance to every cable exterior **and
the whole air/soil interface**. Their minimum keeps a fine band across the
complete domain. The first candidate replaces that unbounded interface source
with finite projected cable footprints. For cable centre `(xc,yc)` and outer
radius `r`, the footprint is the interface segment with halfwidth
`abs(yc)+r`. Its distance is
`sqrt(y^2 + max(abs(x-xc)-(abs(yc)+r),0)^2)`.
Take the minimum of these distances and the existing cable-contour distance.
Keep the existing medium-specific wave sizes, decay radii and remote ceilings.
The footprint captures the geometric lateral scale of the interface forcing;
it is a mesh prescription, not an exact truncation of the transmitted field.
The existing decay radii retain material, frequency and prescribed Gamma
dependence. Native Gmsh MathEval/Distance/Min/Threshold/Restrict fields suffice:
https://gmsh.info/doc/texinfo/#Gmsh-mesh-size-fields

Start with the existing mixed two-wire 0.1 Hz, 0.1 ohm m, 4.25 cm case.
Snapshot owned export files into the qualification directory, never modify the
interactive bundle. Compare an independently meshed baseline and candidate on
identical native solver sources. Record full unclipped R, X, G, B, node counts,
DOFs, assembly/solve seconds and peak RSS. Check exact PML/path/contour coordinate
hashes before solving. Initial qualification target: no new nonzero component
sign reversals and at most 1% entrywise change of each nonzero R/X/G/B component
relative to the baseline. A zero baseline is reported separately, never clipped
or silently assigned a relative tolerance. This is a test-harness gate, not a
production acceptance rule or a claim of 1% physical accuracy.

If successful, cover air/air, soil/soil and mixed at 0.1 Hz, 1 kHz, 1 MHz with
Gamma=0 (9 cases), plus each placement at 0.1 Hz and 1 MHz with
Gamma=0.99*gamma_earth (6 cases), including the initial case. Retain independent
analytical comparisons where applicable; do not mistake agreement between two
meshes for independent physical validation. Use mesh preservation checks for
screen, tubular and sector fixtures already qualified by earlier work.

If isotropic localization fails, diagnose the failing entries before choosing
a separate candidate. Tangential coarsening with preserved normal resolution
is a possible follow-up, not an automatic engine fallback. Do not execute the
remaining campaign on a demonstrably failing candidate.

## Delivery

After qualification, implement the prescribed fields at their existing owners,
with consistent optional refinement controls in managed computation and native
ONELAB export. Update mesh fingerprints, focused tests and documentation. No
production self-validation, retries, refinement loops or sign correction.
Measure compilation separately from warmed execution before claiming speedup.
Use Julia and native tools only, one serial resumable job and one live log:
`.linecablemodels/fem/local-interface-mesh/live.log`.

Completion requires qualified implementation in both execution paths, recorded
limitations and measured costs. A mesh-only pilot or completed first solve is
not completion of this task.

## Execution record

The first mixed 0.1 Hz pair passed the prescribed comparison: maximum component
change 0.7397%, no sign changes; self-R changes approximately 0.0011%. Nodes
101290 -> 91518, DOFs 400546 -> 361458. Native process seconds 30.21 -> 28.29
and peak RSS 2003160 -> 1815564 KiB are single-run observations. PML, native
paths and conductor contours have exactly matching coordinate hashes. The
baseline itself has substantial relative disagreement with the analytical
reference for the tiny aerial self-G; this optimization is not a correction
of that pre-existing discrepancy. Spectrum qualification is proceeding.

The second case (air/air, 0.1 Hz) retained all signs but exceeded the initial
1% component-change gate: G12 changed 1.6239%. The remaining spectrum was not
run. Air-only localization passed (0.5463% maximum change); soil-only failed
(2.2273%). Their node counts were 95952 and 95936 versus baseline 100844.
These results motivate preserving more of the soil interface resolution.

A separate native anisotropic metric retained the old normal target and used
the localized tangential target with per-surface BAMG. It was rejected on cost
before any GetDP solve: BAMG generated 900000 vertices on one physical face,
far exceeding the entire reference mesh. The process was interrupted; this
implementation is not a viable feature candidate.

The next isotropic candidate expands each interface footprint by the existing
medium wave-decay radius, already computed from material/frequency/Gamma and
capped by the resolution radius. Thus halfwidth is
`abs(yc)+r+wave_decay_radius[medium]`. This retains a lateral diffusion-scale
neighbourhood of the transmitted field before coarsening. No new fitted length,
changed PML, or altered conductor target is introduced. Test it first on the
failed aerial probe, then the mixed probe and remaining spectrum if successful.

The expanded footprint also exceeded the aerial comparison target (1.4511%).
An identical-mesh control changing only MUMPS factorization ordering to AMF
changed components by at most 6.2232e-7 relative, ruling out this particular
factorization-roundoff explanation for the 1-2% sensitivity.

The next candidate uses the existing wave bound to limit localization:
replace the full interface with geometric cable footprints only when that
medium's remote size is no larger than its existing `wave_size_limit`.
Otherwise retain the full interface source unchanged. The remote air cap is
already bounded this way; the earth cap generally is not at Gamma=0. This
distinction explains why air-only coarsening passed and earth-only coarsening
failed. The prescribed rule also accounts for finite Gamma through the
existing q-dependent limits. It does not inspect computed solutions or retry.
Qualification still precedes feature implementation and retains the same gate.

Conditional implementation for this surviving candidate: add one prescribed
`interface_refinement_factor` (default 1, at least 1) multiplying the geometric
footprint halfwidth. Increasing it retains more interface refinement. A medium
whose remote cap exceeds its wave-size limit retains the full interface source.
Expose the same value in detached ONELAB and evaluate that condition again when
the native physical/exterior mesh factors change. Keep the existing remote,
conductor, path and PML controls authoritative. Include the new factor in mesh
identity. This is a direct change in the existing size-field owner, not a new
meshing backend, physical formulation or solution-dependent selection.

All nine Gamma=0 comparisons passed the unchanged 1% component-change gate,
without new sign reversals. The mixed 0.1 Hz finite-Gamma case also passed
(0.3565%). The mixed 1 MHz finite-Gamma case exposed a distinct limitation:
its air wave field has equal minimum/maximum sizes, and its earth cap does not
permit localization. Thus this optimization has no active size-field change.
Nevertheless independently remeshing the untouched prescription changed one
weak mutual component by 56.1% (G12 by 23.1%), with two extra mesh nodes.
The physical medium triangle counts, PML/path/contour coordinates and every
prescribed size were retained. Rewriting the unnecessary constant field gave
73.9% maximum change. These observations are preserved, not counted as passes.

Do not mutate constant fields. For qualification points where no size field is
changed, verify that the source files match and reuse the exact baseline mesh
for the numerical comparison. This tests preservation without conflating it
with pre-existing remeshing sensitivity. Active localization still requires
the original componentwise 1% gate and unchanged signs. The final report must
distinguish these cases and explicitly retain the finite-Gamma remeshing
limitation; this feature does not establish finite-Gamma mesh convergence.

The final 15-case selection completed with the documented no-change mesh
controls. Screen/tube/sector resolution, PML and measurement-boundary checks
completed before feature implementation. The prescribed rule and optional
factor are now implemented in the existing managed and native-export owners.
Maintained checks passed 2408 assertions; public zero/finite-Gamma managed,
detached and resume checks passed 20 more. Both compilation and warmed public
execution were measured separately. The matched native warmed timing repeat
completed: mixed 0.1 Hz, Gamma=0 took 30.12 -> 29.18 s, with peak RSS
2005408 -> 1904052 KiB. The original pair took 30.21 -> 29.05 s. This is a
modest measured saving, not a general performance guarantee. Implementation,
qualification, execution-path verification and documentation are complete;
the finite-Gamma remeshing limitation above remains explicitly unresolved.
No scientific feature gates were added to production.

## Both-media correction authorized 2026-09-30

The user correctly rejected the preceding air-only outcome as incomplete: the
wave-cap condition leaves a full-width soil band at Gamma=0. Keep the frozen FEM
as reference. Qualify actual localization in both media with the now-authorized
2% per-component R/X/G/B comparison and no new sign reversals. Start with compact
footprints, reuse the two completed candidates, and complete the established
15-case selection. If compact fails, check the wider decay-length footprint on
that case first before continuing its selection. All acceptance remains in the
manual harness. Baseline meshes, sources and solver results are read only.

Only after qualification replace the production wave-cap condition in both
managed and native-export paths. Preserve conductor targets, Cartesian geometry,
PML nodes and measurement nodes. Deliver before/after mesh images, cable-shape
constraint checks, component comparisons and repeated native timing. Completion
requires actual localization in both media; another air-only prescription does
not fulfil the request. Use the same serial append-only live log.

Compact localization passed all nine Gamma=0 points and mixed finite-Gamma
0.1 Hz. It failed mixed 1 MHz, Gamma=0.99 gamma_earth (57.5% maximum component
change); the wider footprint failed there too (66.7%). Neither changed signs.
The four remaining finite-Gamma cases were not run. The full 2% qualification
is unmet; production is unchanged. Both-media screen/tube/sector mesh contracts
passed. A matched repeat established lower memory, but no dependable speedup.
The failed case has mathematically unchanged wave-size targets throughout both
physical media, exposing the pre-existing remeshing sensitivity.
See `test/manual/calculations/fem_local_interface_mesh/both-media-results.md`.
This does not complete the feature.

## Finite-Gamma investigation and correction

The user selected investigation before release. Fixed-mesh controls ruled out
factorization reuse and exposed sensitivity to LU ordering and refinement.
An exact metal-drive substitution, `w = u - Gamma^2 V`, cancels identical metal
basis contributions before assembly and restores the physical drive for output.
On the failing mixed 1 MHz case this reduced compact remeshing differences from
57.5% to 0.0003461%; wider remeshing gives 0.0014414% and changing LU ordering
gives 0.0001717%. No mesh, PML, gauge or voltage path was changed in that control.

Before production edits, all six finite-Gamma independent-mesh comparisons
passed (maximum 0.234699%, no new signs). The Gamma=0 control was bitwise
identical with and without substitution. The earlier nine Gamma=0 mesh controls
also pass the 2% limit. Original references remain immutable; the finite-Gamma
comparison with corrected equations on both meshes is explicitly separate from
preservation of the formerly unstable coefficients. This is not a tolerance
relaxation or a claim of analytical accuracy.

The qualified substitution and both-media footprints are now applied in the
shared native formulation and existing managed/export mesh owners. Mesh cache
identity changes accordingly. Public execution/resume (42 assertions), six
maintained regression items, native finite-Gamma current-map integration and
production screen/tube/sector mesh-preservation checks passed. The first public
high-frequency check exposed a pre-existing 2.0655% difference between managed
and detached construction. A captured pre-edit managed builder established
0.0001547% before/after localization change; the public preservation check uses
that matched reference. The original cross-route failure remains documented,
and the 2% limit is unchanged. Implementation and verification are complete;
analytical errors and independent-mesh discretization differences remain explicit.
See `test/manual/calculations/fem_local_interface_mesh/finite-gamma-diagnosis.md`.
