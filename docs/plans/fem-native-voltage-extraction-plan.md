# Native voltage extraction: description and Julia completion plan

Status: cleanup completed and validated on 2026-09-28 after execution was
authorized. See the [execution evidence](fem-native-voltage-cleanup-validation.md)
for the removals, 428/428 focused checks and independent P/Y refinement results.
The description, pre-execution inventory and acceptance criteria below preserve
the original plan; they should not be read as outstanding implementation work.

The objective is to eliminate all Julia-owned, mesh-dependent voltage-path
preparation. Gmsh owns measurement geometry and its one-dimensional mesh;
GetDP owns field evaluation, quadrature, reference subtraction and normalized
voltage extraction. Julia builds the model, launches GetDP and consumes its
completed numerical columns. Both the Julia and detached ONELAB routes use the
same maintained GetDP formulations.

The current working tree already implements this production boundary. This is
not merely an exporter demonstration waiting to be ported: the Julia worker
also uses the native measurements. Completion therefore means confirming that
boundary, removing stale executable consumers and documenting and validating
the numerical contract. Do not implement a second native path.

## 1. Verified current state

| Owner | Source evidence | Current behavior |
|---|---|---|
| Geometry | [`geometry.jl`](../../ext/LineCableModelsGmshExt/geometry.jl), `_build_geometry!`, measurement block around line 1502 | Creates one oriented vertical path and one reference point per terminal, before meshing |
| Mesh contract | [`mesh.jl`](../../ext/LineCableModelsGmshExt/mesh.jl), `_expected_physical_groups`, `_inspect_loaded_mesh`, `_mesh_fingerprint` | Requires path/reference physical groups; mesh partition version is 7 |
| Julia orchestration | [`compute.jl`](../../ext/LineCableModelsGmshExt/compute.jl), `_headless_solve!` | Geometry → mesh selection → ordinary input writing → GetDP workers → result validation |
| Worker launch | [`workers.jl`](../../ext/LineCableModelsGmshExt/workers.jl), `_getdp_command`, `_frequency_job` | Passes mesh, model data and requested bases; no `PathDataPath` or measurement-preparation call |
| Quasi-TEM extraction | [`quasi-tem.pro`](../../ext/LineCableModelsGmshExt/getdp/quasi-tem.pro), `FEMAppendRaw` | Terminal scalar potential minus the receiving terminal's scalar reference |
| Quasi-full-wave extraction | [`quasi-full.pro`](../../ext/LineCableModelsGmshExt/getdp/quasi-full.pro), `bt_mesh`, `ReVoltageLine`, `ImVoltageLine`, `FEMAppendRaw` | Scalar difference plus native complex vector-potential circulation |
| Removed Julia implementation | `ext/LineCableModelsGmshExt/voltage_paths.jl` | Absent on disk and marked deleted relative to HEAD; no extension include or adapter-source registration |
| Run identity | `compute.jl`, `_fem_input_record`, `_same_fem_inputs` | Input schema 8 and solver protocol 3; old input identities are not reusable as current results |
| Remaining Julia algebra | [`results.jl`](../../ext/LineCableModelsGmshExt/results.jl) | Parses native P, applies established matrix processing and uses Julia linear algebra for Y; this is separate from voltage-path preparation |

The function `_prepare_run_inputs!` still exists. It copies/verifies GetDP source
assets, writes `problem.json` and `model_data.pro`, creates output directories
and initializes raw-table headers. It neither reads mesh triangles nor creates
path samples, interpolation weights or reference quadrature. Keep these
necessary input and reproducibility operations. Deleting or renaming everything
containing “prepare” would not advance the objective.

This inventory describes the live checkout. It does not establish that another
branch, an already loaded Julia session, an installed package, or a retained
run contains the same implementation.

## 2. What the removed preparation used to do

The previous implementation reopened each completed mesh in a temporary Gmsh
model, read mesh nodes and triangles, classified metal versus field media and
constructed voltage measurements from that mesh. Its responsibilities included:

1. Selecting conductor-contour sample positions. The last pre-removal version
   averaged bare overhead disk receivers using normalized contour quadrature;
   other receivers used a lowest contour node.
2. Locating overhead references on the air/earth interface.
3. Intersecting vertical measurement segments with volume triangles, excluding
   metal interiors, checking coverage and assigning integration midpoints and
   oriented displacement weights.
4. Writing `RefStart`, `RefCount`, `RefX`, `RefWeight`, `PathStart`, `PathCount`,
   `PathX`, `PathY`, `PathDX` and `PathDY` into `paths-fNNNN.pro`.
5. Passing that generated file to GetDP as `PathDataPath` for an explicit sum
   over externally prepared evaluations.

The deleted implementation and its callers are recoverable from version
history; the immediately pre-removal copy was also inspected under
`/tmp/lcm-native-replacement-backup/ext/LineCableModelsGmshExt/voltage_paths.jl`.
That temporary copy is evidence, not a future runtime or validation dependency.

The new implementation eliminates that owned mesh traversal, clipping,
quadrature serialization and parameter transport. It does not eliminate the
mathematical need to define a path or the numerical work of evaluating a field.

## 3. The implemented native measurement

### Geometry and orientation

Coordinates are x horizontal, y vertical and z along the line; the interface is
y = 0. The receiving terminal owns the reference for every excitation column.

For each terminal, geometry construction inspects its existing CAD contour
vertices and chooses the one with the smallest y, breaking ties by smallest x.
Call this endpoint C_i. This is a CAD vertex, not a sample selected from the
finished triangle mesh, and not necessarily an arbitrary shape's exact minimum
over every point of its contour. Disconnected pieces assigned to one terminal
still produce one chosen endpoint and share the terminal's scalar unknown.

An overhead terminal has R_i = (x(C_i), 0). A buried terminal has
R_i = (x(C_i), -H - d_bottom), at the outer bottom PML boundary. The overhead
path has one geometric segment. The buried path has two: bottom boundary to
the physical/PML interface at y = -H, then that interface to C_i. Every segment
is oriented from reference toward conductor.

The `.msh` contains physical groups:

```text
LCM/voltage_path/0001         dimension 1
LCM/voltage_reference/0001    dimension 0
```

The corresponding GetDP regions are `VoltagePath~{t}` and
`VoltageReference~{t}`. Physical tags come from the existing model tag owner.
The curves are not embedded into the two-dimensional surfaces. Their line
elements can cross triangle interiors; they do not alter the solved volume
mesh merely to align an integration path.

Gmsh meshes these curves using native transfinite progression. Physical-domain
segments reuse the existing exterior grading rule; buried PML segments reuse
the PML layer progression. This avoids applying the smallest conductor spacing
uniformly over a potentially very long path.

### Quasi-TEM: use the scalar solution directly

Here E_t = -grad(v). Therefore a voltage difference already represents the
electric-field line integral:

```math
v(C_i)-v(R_i)=-\int_{R_i}^{C_i} E_t\cdot d\ell
             =\int_{C_i}^{R_i} E_t\cdot d\ell.
```

For excitation j with transverse current q_j [A/m], the implemented coefficient
is:

```math
P_{ij}=\frac{V_i^{(j)}-v^{(j)}(R_i)}{q_j}.
```

`V_i` is the grouped equipotential terminal unknown. An overhead reference is
evaluated natively using `OnRegion VoltageReference~{i}, Dimension 2`. The
dimension selects the solved two-dimensional field for the point evaluation.
A buried reference is zero by the outer Dirichlet condition. The code does not
numerically integrate the gradient to reproduce an already available scalar
difference.

The air/earth interface is not grounded by this operation. Its overhead
reference value is part of the solved field. Both formulations retain zero
scalar potential on the prescribed outer boundary, including the outer air
boundary; extraction does not change those boundary conditions.

### Quasi-full-wave: add native circulation

This formulation uses the normalized variables A_t = Gamma b_t and
phi = Gamma v. It solves the established first-order longitudinal reduction,
not a finite-Gamma eigenproblem. The normalized transverse electric field is

```math
E_t/\Gamma=-\nabla_t v-\mathrm{j}\omega b_t.
```

With imposed axial current I_j [A], native extraction computes

```math
P_{ij}=\frac{v_i^{(j)}-v^{(j)}(R_i)
       +\mathrm{j}\omega\int_{R_i}^{C_i}b_t^{(j)}\cdot d\ell}{I_j}.
```

The scalar potential alone therefore does not determine the complete voltage
in this branch. A buried terminal also retains the vector integral, even
though its scalar reference value is zero. The sign agrees with integrating
the electric field from conductor toward reference.

After each solved excitation, `FEMAppendRaw` executes the following native
operations, regardless of whether user-visible field maps are enabled:

1. Store the current solution's `bt_mesh` on `DomainMedia_Ele` in the private
   in-memory field `FEMVoltageField = 990001`.
2. Evaluate that complex field at GetDP's line quadrature coordinates.
3. Dot it with the oriented line tangent and integrate over the receiver's
   physical curve group, separately retaining real and imaginary parts.
4. Evaluate the scalar reference and terminal unknown, normalize by the actual
   source amplitude, and write the completed P column.

The operative syntax already present in the formulation is:

```c
Print[bt_mesh, OnElementsOf DomainMedia_Ele, Format Gmsh,
  File "", StoreInField FEMVoltageField];

{ Name ReVoltageLine; Value { Integral {
  [Re[ComplexVectorField[XYZ[]]{FEMVoltageField} * Tangent[]]];
  In VoltagePaths; Jacobian VoltageLine; Integration VoltageGauss;
} } }

Print[ReVoltageLine[VoltagePath~{response_terminal}], OnGlobal,
  File "", StoreInVariable $FEMLineRe];
```

`ImVoltageLine` is the corresponding imaginary expression. `VoltageLine` uses
`Jacobian Sur` for the line measure in this two-dimensional problem;
`VoltageGauss` uses four Gauss points per line element. The returned complex
integral L is combined as `j*omega*L`: its real contribution is
`-omega*Im(L)` and its imaginary contribution is `+omega*Re(L)`.

GetDP documents the in-memory field store, coordinate-based complex field
evaluation, line tangent and native integrated postquantities in its
[PostOperation](https://getdp.info/doc/texinfo/getdp.html#Types-for-PostOperation),
[function](https://getdp.info/doc/texinfo/getdp.html#Miscellaneous-functions) and
[PostProcessing](https://getdp.info/doc/texinfo/getdp.html#PostProcessing)
references. The call sequence above is established by this repository's owned
source, not inferred from an `OnLine` example.

There is no `.pos` round trip. `File ""` plus `StoreInField` supplies an internal
GetDP/Gmsh field; `PlotFieldMaps` controls separate inspectable output files.
The store occurs inside each basis's extraction, so later excitations do not
reuse the first excitation's field. This requires GetDP built with Gmsh support,
as supplied by the native binaries already used by this stack. It requires no
ONELAB server for a Julia-driven solve.

### Metal, PML and normalization

`DomainMedia_Ele` contains air, earth, passive material and their applicable PML
regions, excluding equipotential metal interiors. The stored b_t field has no
support in metal. Native field lookup contributes zero where the selected
field has no value; it supplies the intended zero transverse contribution
through metal. This is why no external triangle classification is needed for
the integral. The same lookup behavior also means a malformed path outside
all field support must not be mistaken for a valid metal crossing; geometry
and coverage validation must exercise that distinction.

In the Cartesian PML, b_t is already the pulled-back one-form. Integrate it
against the mesh-coordinate line displacement once; do not add another
stretch tensor or volume determinant to the circulation. The finite outer
boundary approximates the remote reference, with PML truncation error assessed
separately.

P is the complex inverse-admittance coefficient [ohm m], not the scalar
potential field. Y = P^-1 [S/m]; no additional j*omega factor belongs in that
inversion. The separate quasi-full `Pscalar.tsv` contains the terminal scalar
trace normalized by current. It is diagnostic and is not the complete P.
Julia's established reductions/inversion can remain in Julia; relocating them
is not required to remove voltage-path preparation. Detached export already
has its own native algebraic GetDP execution for those operations.

## 4. Numerical changes and cost

The replacement made two distinct choices. First, it replaced external
triangle-partition integration by native quadrature on independent line
elements. Second, it replaced bare-overhead contour averaging with one explicit
path and one surface reference per terminal. These must be evaluated separately.
Do not claim bitwise preservation of the former contour-averaged observable.

GetDP evaluates the field at quadrature coordinates; it does not snap the path
to triangle edges or promise exact integration across every triangle crossing.
Four-point Gauss integration is not automatically exact when a line element
crosses discontinuities in the piecewise field. Refining the measurement mesh
controls that integration error. Export exposes `VoltageRefinements` for this
purpose; the current Julia public options do not expose that exporter control.
Validation can refine only the measurement curves in a test-owned mesh without
introducing another production preparation stage or a new public API.

The recorded [native measurement study](fem-native-measurements.md) compared
one-path native extraction with a matched-endpoint triangle-partition reference
on a fixed coarse mixed-wire solution. Its maximum relative P difference was
1.122896% at the original line grading and 0.018727% at 32 subdivisions relative
to that grading, with nonmonotone intermediate changes. These are previously
recorded results, not measurements repeated while writing this plan, and not
physical accuracy bounds for arbitrary cable systems.

The removal eliminates Julia's extra mesh import, triangle arrays, clipping
passes and generated quadrature files. Native work remains: one-dimensional
meshing and, for each quasi-full excitation, building a field representation,
native spatial lookup and quadrature. No timing or memory speedup is claimed
from source inspection. Any future performance claim must report first-use
and warmed Julia timings, plus GetDP extraction time and peak memory on the
same mesh, source columns, field-map setting and measurement convention.

## 5. Pre-execution cleanup inventory (resolved)

| File | Verified stale dependency | Required disposition |
|---|---|---|
| `test/manual/calculations/run_fem_voltage_reference.jl` | Calls `_prepare_voltage_paths!`, `_voltage_path_file`, `_voltage_path_points`; instruments generated arrays | Retire the legacy preparation-based solve modes. Preserve saved-result reading/plotting and historical numerical records. Provide any current extraction/convergence study through native measurement fixtures, without restoring private helpers |
| `docs/src/fem.md`, manual voltage-reference paragraph | Advertises `run_fem_voltage_reference.jl --full` as a current workflow | Replace the stale invocation with the maintained native measurement checks and their documented convergence procedure |
| `test/manual/calculations/fem_static_boundary/mesh.jl` | Five-argument branch calls `_write_voltage_paths` | Remove that obsolete preparation mode; retain ordinary native geometry meshing. Mark the historical dependent experiment as requiring its recorded source revision |
| `test/manual/calculations/run_fem_terminal_extraction.py` | Imports removed exporter `measurements` and `matrices` modules and supplies `PathDataPath` | Retire the obsolete solve entrypoint; preserve its existing data and description as historical evidence, without recreating Python modules |
| `test/manual/calculations/fem_static_boundary/run.py` | Supplies `PathDataPath` and relies on the obsolete meshing mode | Retire its preparation-dependent invocation from current instructions; preserve the historical experiment and results |
| `test/manual/calculations/fem_cartesian_compactification/run_cable.py` | Supplies old generated path files from retained runs | Treat as a historical source-snapshot experiment, not a current backend example; remove obsolete current-backend invocation claims |
| Historical research Markdown | Describes clipping, contour averages and generated paths as the then-current implementation | Keep numerical history, clearly identify the old source/measurement convention and link to the current native contract |

Retiring means removing the executable legacy preparation route and its current
usage instructions, while keeping scientific records and independent
saved-result inspection. It does not mean adding a runtime fallback or copying
the deleted helpers into a “test utility”. Porting unrelated boundary-condition
or compactification research to a new formulation is outside this completion
task. Do not delete unrelated Python plots, experiments, retained runs or
manually edited exported projects.

## 6. Ordered execution plan

1. **Establish the exact target revision and ownership.** Record `git status`,
   the task file manifest and the active package path. Inspect the call chain
   above in a fresh Julia process when execution is authorized. Preserve all
   unrelated work. If this checkout is the target, record the production
   removal as already present. An older target must receive the same native
   geometry and shared formulation changes before its old path is deleted.
   Exit: one verified production call chain and a bounded change manifest.

2. **Lock the measurement contract.** Retain reference-to-conductor orientation,
   one CAD endpoint per terminal, the host-dependent reference, source units,
   metal exclusion and PML pullback. Use the quasi-TEM and quasi-full expressions
   above. Preserve field equations, boundary/gauge constraints, physical ordering
   and result shape. Record the deliberate departure from contour averaging.
   Exit: equations, geometry and native postoperations agree on the same P.

3. **Close production removal, without creating another stage.** Verify absence
   of the deleted module, extension include, adapter digest, call to
   `_prepare_voltage_paths!`, `PathDataPath`, generated path arrays and file
   checksums. Verify workers directly consume ordinary `.msh` and model data.
   Keep `_prepare_run_inputs!`, `bases.pro`, input snapshots and result validation.
   Keep one shared pair of GetDP formulations. Preserve the existing mesh and
   run-identity invalidation; only advance versions again if a further semantic
   change actually requires it. Exit: no production measurement preparation
   dependency or compatibility path remains.

4. **Remove stale executable consumers.** Apply the per-file dispositions in
   section 5. Preserve plotting of retained historical outputs where it is
   independent of the deleted runtime. Current study instructions must point
   to native fixtures; retired study results must retain their original
   reference/averaging labels. Do not silently compare an old averaged column
   with a new single-path column as if they measured the same functional.
   Exit: no runnable Julia or advertised current workflow calls the removed
   helpers; no current workflow imports deleted exporter Python modules.

5. **Fill meaningful validation gaps.** Reuse the existing native integration,
   manufactured extraction, geometry, worker and exporter tests. Add only the
   missing orientation, field-support and independent-refinement checks below.
   Keep convergence and performance experiments outside production. If a new
   check exposes a defect, correct the owning geometry/postoperation; do not
   restore a clipper, loosen a tolerance or insert automatic retries.
   Exit: each acceptance criterion has reproducible evidence.

6. **Verify execution and report completion.** Run the bounded commands below,
   check a fresh run's inputs/worker arguments and retained raw output, and
   repeat the direct/export shared-mesh comparison. Publish the changed-file
   manifest, exact commands, binaries and outcomes, separating mathematical
   checks, physical convergence and performance. Do not overwrite historical
   studies or external exports to make them appear current. Exit: all required
   gates pass or the remaining scientific/implementation issue is explicitly
   reported; no unsupported “fully validated” claim.

## 7. Acceptance and validation criteria

| Gate | Validation | Pass criterion |
|---|---|---|
| V01 Production removal | Inspect extension registration, input capture and worker commands; inspect a fresh retained run | No preparation helper, `PathDataPath`, generated `paths*.pro` or serialized path/reference weight arrays; ordinary input snapshots remain |
| V02 Geometry contract | Overhead, buried, mixed, insulated, noncircular and disconnected-terminal cases | Correct groups, host reference, deterministic CAD endpoint and reference-to-conductor orientation; buried PML split; no forced volume embedding |
| V03 Native complex integral | Existing analytic edge-field fixture on a nonconforming open path and closed loop | For `(1+2j)*(1-y,x,0)`, retain exact-reference checks `0.24+0.48j` and `2+4j` at absolute tolerance `1e-13`; off-mesh field check unchanged |
| V04 Orientation and field support | Reverse a native test path; include a path through a known zero-support metal interval and a separate malformed outside-domain case | Reversal changes the integral sign; metal interval contributes zero; the malformed geometry is identified as invalid, not accepted as metal. Test-level geometry checks suffice; no new preparation pipeline |
| V05 Production extraction | Existing basis-dependent manufactured scalar/vector values, both physics choices, non-unit source amplitudes and maps disabled | Existing `rtol=1e-11` expectations pass for every response/source pair; correct reference and j*omega sign; different bases refresh the native field; no `.pos` dependency |
| V06 Route preservation | Direct Julia and detached GetDP on the same mesh, frequency, source and measurement convention | Retain componentwise bound `abs(a-b) <= 2e-9*abs(b) + 100eps(Float64)*maximum(abs,b)`, separately for real and imaginary primitive P/Z and resulting Y/Z; no tolerance relaxation |
| V07 Independent line refinement | Hold solved volume field/geometry fixed, refine only measurement curves | Volume coordinates/connectivity are unchanged; report P and Y sensitivity separately. For a manufactured piecewise field, include a final line mesh aligned with its known discontinuity and require its analytic integral within `1e-13` absolute error. Record the coarser nonconforming errors without requiring monotonicity. No universal engineering tolerance is inferred from the toy |
| V08 Mesh/restart validity | Cached mesh, supplied mesh, old run identity, partial and complete restart | Current groups accepted; old meshes lacking groups rejected/remeshed through existing policies; incompatible historical results not adopted; valid columns remain reusable |
| V09 Worker ownership | Multiple frequencies/workers, factorization reuse, maps on/off | Same numerical columns and ordering; independent process fields; no new preparation subprocess or generated measurement dependency |
| V10 Complete cleanup | Static audit of production, maintained tests and advertised manual entrypoints | No live old helper calls or deleted-module imports; remaining text matches are explicit historical descriptions or negative tests |
| V11 Scientific scope | Compare matched endpoints separately from historical contour-averaged results | No claim that the changed observable is identical; scalar diagnostic never substitutes for full quasi-full P; physical/PML convergence is reported independently |

The fixture in V03 establishes native interpolation, orientation and complex
arithmetic. V05 establishes the actual production extraction and normalization.
V06 establishes transport/route preservation. None alone establishes physical
accuracy of arbitrary production meshes. For an engineering accuracy claim,
use the user's study-specific accuracy criterion fixed before inspecting the
new results; this plan does not invent one.

Use the existing small mixed-wire toy: 5 mm-radius wires at (0, 0.1) and
(0.2, -0.1) m, earth rho = 100 ohm m and relative permittivity 10, frequencies
50 Hz and 10 kHz, eight PML layers, no bundle/Kron/transposition reduction.
Exercise both formulations and both source columns. The small PML count is
appropriate for execution and preservation evidence; it is not an engineering
reference mesh. Supplement it with geometry-only insulated/polygon/disconnected
cases and the independent manufactured fields above.

## 8. Execution commands

These are the original execution instructions; the linked execution record
identifies the commands actually run and their outcomes. Run from the
repository root with the supported test project. Use the repository's
normal dependency preparation once; on this host the already established
writable depot can be selected with
`JULIA_DEPOT_PATH=/tmp/lcm-onelab-depot:/home/amartins/.julia`.

```sh
# Audit exact legacy dependencies; classify negative tests and historical prose.
rg -n '_prepare_voltage_paths!|_write_voltage_paths|_voltage_path_points|_voltage_receiver_quadrature|_voltage_path_file|PathDataPath' ext src test docs
rg -n 'from (measurements|matrices)|include\("voltage_paths.jl"\)' ext src test

# Confirm this process is using this checkout, rather than another installation.
julia --project=test --startup-file=no -e 'using LineCableModels, Gmsh; println(pkgdir(LineCableModels)); println(pathof(Base.get_extension(LineCableModels, :LineCableModelsGmshExt)))'

# Focused existing scientific and execution checks, extended only where needed.
julia --project=test --startup-file=no test/runtests.jl \
  fem_native_measurements fem_quasi_full fem_mesh_grading \
  fem_workers fem_resume fem_export

# Emit a new disposable native bundle; do not overwrite a historical study.
julia --project=test --startup-file=no test/manual/fem/export_onelab_toy.jl \
  /tmp/lcm-native-voltage-acceptance

# Optional inspection of the emitted native model after numerical checks.
gmsh /tmp/lcm-native-voltage-acceptance/study.pro
```

Execute numerical toy solves through the maintained tests and the exported
native Run controls or direct Gmsh/GetDP commands documented in the bundle.
Record per-frequency meshes and selected controls for matched comparisons;
the native ONELAB scan overwrites its working mesh at each frequency.
Manual research files must remain excluded from automatic test discovery.

Completion means that current Julia-driven FEM computation has one native
measurement implementation, no externally prepared voltage-path data, no
stale executable callers presented as supported, and evidence for the stated
numerical contract. The linked execution record closes this completion plan.
