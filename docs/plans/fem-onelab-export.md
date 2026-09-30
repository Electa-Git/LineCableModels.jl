# Native ONELAB export: execution plan

Date: 2026-09-28. Status: implemented and validated; A01-A14 satisfied.
Execution evidence and remaining numerical/platform limits are recorded in
[fem-onelab-export-validation.md](fem-onelab-export-validation.md).

**Goal:** `export_data(:onelab, problem_or_system, formulation; ...)` produces
an editable, relocatable project that meshes, solves, integrates voltage paths,
assembles and reduces matrices, computes admittance, and displays results using
only Gmsh, GetDP and ONELAB. Julia generates the project and has no runtime role.
No Python interpreter, Python package, driver, shell launcher, Julia process,
LineCableModels installation or repository checkout is needed to use the export.
The applications' ordinary bundled/native libraries are allowed dependencies.

This is completion of the intended feature. The existing driver was a transient
prototype. It imposes no compatibility contract: delete its runtime, interfaces,
installation instructions and tests after replacing the useful behavior natively.
Do not retain a second execution path, fallback, deprecation layer or migration
framework. Earlier prototype validation is evidence about individual components,
not acceptance of this goal.

The endpoint is the implemented and validated exporter, including actual ONELAB
interaction and a freshly generated usable toy. A plan, algebra fixture, working
CLI alone, or partial driver removal cannot close the implementation goal.

**User-visible contract**

- Keep the existing Julia export dispatch for `LineParametersProblem` and the
  `LineCableSystem` convenience method. Export itself never runs a solver or GUI.
- `gmsh study.pro` opens the model with readable geometry and working native
  Check, mesh, Run and Stop behavior. GetDP appears under Gmsh's normal solver
  registration, using its configured/discovered executable. The project creates
  no driver registration and no additional executable-path parameter.
- Geometry and evaluated material frequency cases are fixed by Julia at export.
  The user can select an exported frequency, either supported formulation,
  one excitation or the full matrix, mesh refinement, field outputs, solver
  settings and per-metre/total output. Named material/source/connection values
  and constraint assignments remain directly editable in native source files.
  Arbitrary frequency-dependent Julia material functions are not exported as
  executable functions; the frequency table is explicit.
- Points, curves, surfaces, physical groups, materials, regions/domains,
  terminals, references, constraints, equations, integration rules and matrix
  operations are visible in ordinary `.geo`/`.pro` syntax. Comments retain the
  connection to cable/layer/terminal names. No JSON or generated runtime code
  becomes an alternative model authority. Runs never rewrite model sources.
- Both `quasi_tem` and `quasi_fw`, primitive and reduced matrices, bundle/Kron
  options, ideal transposition, terminal ordering, complex phasor convention,
  units and normalization remain supported. A selected excitation produces
  diagnostic columns; it cannot produce a purported full inverse.
- One exported frequency case is selected per native run, matching the present
  interaction. Every exported case is runnable. A batch scheduler is unnecessary.
- Existing headless Julia computation and separate Makie mesh/field inspection
  remain supported. The exported project also works without an ONELAB server via
  direct Gmsh/GetDP commands.

**Native ownership and files**

| Responsibility | Owner in the completed feature |
|---|---|
| Problem dispatch, fixed geometry, frequency/material tables, file emission | Julia exporter in `ext/LineCableModelsGmshExt/export.jl` |
| Primitive geometry, physical groups, mesh size/grading/refinement | `study.geo`, `geometry/case-*.geo`, native Gmsh commands |
| Named coefficients, connections, excitation and ONELAB controls | `study_data.pro`, consumed directly by Gmsh/GetDP |
| Entry, effective case selection, public resolution | `study.pro`, directly including native files |
| FEM regions, materials, constraints, spaces and equations | The maintained `getdp/*.pro` sources copied through `FEM_GETDP_SOURCES` |
| Voltage reference and field circulation | Shared native GetDP post-processing on physical points/curves |
| Matrix collection, reduction, algebraic solves and output | One owned native `getdp/line-parameters.pro` implementation, included in the exported project |
| Executable selection and process start/stop | Standard Gmsh/GetDP ONELAB integration |
| Results, matrix entries, units and maps | Native GetDP `PostOperation`; native ONELAB publication and Gmsh views |
| Export overwrite ownership | Existing exporter file manifest, with removal of obsolete owned files |

Do not fork the FEM equations for export. Expose the existing solve/measurement
operations once and invoke them from the detached native resolution. Its small
algebraic systems belong alongside the field system. Preserve the direct Julia
worker resolution and its raw-output contract; do not make Julia compute use
ONELAB. Extract a shared macro only where the two resolutions actually need the
same operation sequence, without callback registries or generic driver machinery.

The main entry defines its own data path, selected frequency/PML values,
requested basis list and output directory. Direct invocation must not require
knowledge of `DetachedSolve`, run snapshots, `bases.pro`, worker internals or
long lists of injected coefficients. Remove the export's recursive self-launch.

Target ordinary usage, to verify verbatim during implementation:

```sh
gmsh study.pro
```

```sh
gmsh study.geo -setnumber BuildMesh 1 -0
getdp study.pro -msh study.msh -solve LineCableModelsFEM
```

Native command-line parameters select nondefault exported cases. Document the
exact verified flags alongside the resulting GUI controls. Opening/checking
must not mesh or solve. Use the documented Gmsh automatic meshing/action
mechanism and native refinement operations; verify that the requested uniform
refinement is applied once in both entry modes. Do not recreate a mesh cache in
another language. Gmsh owns remeshing after geometry-affecting changes.
[Gmsh solver options](https://gmsh.info/doc/texinfo/gmsh.html#Solver-options)

**Scientific calculation sequence**

1. Read the selected case and initialize the complex field system. Validate
   ordinary input invariants: frequency index, material-table lengths, terminal
   map, nonzero normalization and mesh physical groups.
2. For each requested basis, update the constraints, solve, and evaluate the
   existing normalized primitive `Z` and `P` entries. Retain field factorization
   reuse where already supported. Store entries in named GetDP runtime variables
   indexed by response and excitation; raw tables are outputs, not an external
   parser protocol required to continue solving.
3. For a complete matrix, apply the same terminal permutation and bundle row/
   column transformations as `Engine.reduce_primitive_matrices`. Then eliminate
   the prescribed rows/columns by solving `M_ee X = M_ek` and forming
   `M_reduced = M_kk - M_ke X`, separately for `Z` and `P`. Preserve the existing
   special treatment of grounded conductors when bundle reduction is enabled
   without Kron reduction. Apply ideal transposition after reduction. Compute
   these maps from the visible connection data whenever GetDP parses a run so
   manual connection edits take effect.
4. Compute `P_reduced Y = I` through a native complex algebraic formulation with
   one global unknown per row and `GlobalTerm` coefficients, followed by native
   `Generate`/`Solve` operations for its right-hand sides. Use this same arbitrary
   size route for every terminal count. Do not introduce a 2-by-2/3-by-3 branch,
   tensor-size ceiling, handwritten elimination or external inversion.
5. Write primitive `Z`, primitive `P`, reduced `Z`, reduced `P` and `Y`, with named
   rows/columns, frequency, real/imaginary parts and units. Primitive `P` and
   reduced `P` are inverse-admittance coefficients in ohm m; no extra `j omega`
   factor enters `Y`. Total output scales `Z` and `Y` by line length and leaves
   the per-unit-length `P` table explicitly labelled.
6. Publish selected scalar/matrix outputs and native `.pos` maps to ONELAB/Gmsh.
   Mark completion only after all requested numerical operations and output
   operations succeed. Keep partial columns visibly distinct from full results.

The arbitrary-size algebra route is supported by the native
[GetDP resolution operations](https://getdp.info/doc/texinfo/getdp.html#Types-for-Resolution).
The existing maintained 4-by-4 fixture verifies the algebra independently.
A further planning experiment fed mutable runtime coefficients into that system
and solved first `P`, then `2P`, in one process: maximum identity residuals were
`4.440892098500626e-16` and `4.441027621704298e-16` on GetDP 3.5. Initialize the
complex system before constructing complex runtime values: doing so before
initialization discarded the imaginary part in the first pass. Preserve this
case as a production integration test. Scratch evidence is under
`/tmp/lcm-native-export-plan/`; it is not an exported dependency.

Do not reproduce NumPy's SVD condition-number machinery just to preserve the
prototype's `inversion.txt`. Report the actual residual and solver diagnostics.
If a condition estimate is reported, identify its norm and estimator. Do not
label a different estimate as the old 2-norm condition number. Singular/nonfinite
results cannot be published as successful; warning/accuracy policy must not be
replaced with speculative retry or acceptance machinery.

**Voltage integration is part of the native model**

Keep the native replacement already implemented: explicit reference points and
oriented measurement curves, in-memory `StoreInField`/`ComplexVectorField`, and
GetDP `Integral` with the appropriate line Jacobian and quadrature. No external
triangle intersection, Python/Julia path preparation, generated weights or
`paths.pro` may reappear. No field-map file round trip is needed for integration.

Retain the actual formulation's scalar-reference plus vector-potential
circulation expression. In quasi-full-wave, a zero scalar potential at the outer
boundary does not by itself remove the inductive line contribution. Overhead
paths start at the air/earth interface; buried paths start at the prescribed
outer bottom reference. Preserve orientation, metal exclusion, terminal/source
normalization and the pulled-back PML convention. The single explicit path
convention already adopted remains visible; do not silently reintroduce contour
averaging.

Native interpolation across volume-element boundaries is numerical quadrature,
not exact triangle clipping. Include an explicit measurement-curve refinement
control implemented in `.geo`, preserving grading endpoints and PML splits. Test
its influence separately from volume-mesh refinement. The previous coarse toy
showed 1.12% extraction difference, decreasing to 0.019% with line refinement;
this is evidence requiring disclosure, not a new default engineering tolerance.
Acceptance of runtime independence does not certify physical mesh convergence.
No existing numerical tolerances may be loosened to obtain a pass.

**Native execution and result lifecycle**

- Share `DefineConstant` controls between the entry files, allowing ordinary
  ONELAB selections while keeping named model coefficients file-authoritative.
  Frequency changes must update geometry, materials and PML together.
- Register a normal GetDP resolution and postoperations. Standard Solver/GetDP
  is the only executable authority; use native thread/options controls without
  promising BLAS-thread semantics that GetDP cannot enforce. Name options
  according to what they actually control and test their effect.
- Use native GetDP runtime variables for matrix data. ONELAB exchanges controls
  and displayed outputs, not intermediate matrices needed for arithmetic.
- Use simple derived output directories distinguished by frequency, formulation
  and full-matrix/selected-basis mode. Native file operations clear a current-run
  completion marker before computation and create it last. Repeated runs replace
  their derived outputs. Do not reproduce Python's snapshot/hash/cache framework.
- Ensure changing mode, mesh or inputs cannot leave old results presented as
  current: clear/invalidate the relevant ONELAB result entries and derived
  publication state on Check/change and before Run. Failure/Stop must leave no
  valid current completion marker. Partial files are permissible but labelled.
- Native files and `.pos` output remain inspectable without Julia. Use standard
  Gmsh result merging; the existing Makie readers continue to consume compatible
  mesh and field artifacts separately.
- Opening and execution use paths relative to the project, not its original
  export location or current shell directory. A stale project `.db` is derived
  ONELAB session state, not a runtime dependency. Test clean and reused sessions.

GetDP provides native output publication through `SendToServer` and runtime
ONELAB parameter operations; use these directly.
[GetDP postoperations](https://getdp.info/doc/texinfo/getdp.html#Types-for-PostOperation)

**Ordered implementation work and completion gates**

| Step | Work | Evidence required before closing the step |
|---|---|---|
| 1. Establish the native entry and data flow | Replace the launcher in a generated toy; initialize the shared FEM system, retain computed matrix entries, and feed them to the algebraic formulation in the same process | A real exported FEM run computes `Y` from its own measured `P`; the runtime-coefficient complex test passes; no external result parser participates |
| 2. Complete native matrix operations | Implement permutation, bundle transformations, native Schur solves, transposition, inversion, units and labels in the owned `.pro` source | Nontrivial complex 1-, 2- and 4-terminal algebra fixtures; all supported reduction flag combinations; selected-basis semantics; comparison against the engine-owned Julia reference |
| 3. Complete Gmsh/ONELAB execution | Native case/mesh/refinement controls, one GetDP registration, Check/Run/Stop, read-only results and field merging; relative output paths | Fresh GUI and CLI runs both work; case changes remesh correctly; native edits affect the next solve; failures and stopped runs cannot masquerade as completed |
| 4. Finish native measurements in export | Ensure every emitted geometry has readable references/paths and native refinement; remove remaining preparation references in maintained callers | Overhead, buried, mixed, polygon and disconnected-terminal coverage; complex line orientation and basis-refresh tests; independently reported line convergence |
| 5. Remove the prototype completely | Delete runtime assets, obsolete controls, Python tests/CI/dependency setup, notices and old instructions; replace them with native tests/documentation | Export-related dependency/reference audit is clean; fresh export has no scripts/interpreter requirement; overwrite removes obsolete owned assets without touching unrelated files |
| 6. Validate and deliver the feature | Run the acceptance matrix below, rebuild the user's representative toy, write final native-only usage and execution evidence | Every required gate passes; a person can open the fresh `.pro`, mesh, run, inspect matrices and fields without Julia/Python or repository access |

Integrate these steps as one feature completion. Intermediate prototypes are
allowed for development, but cannot be reported as the completed exporter.

**Exact removal and update scope**

Delete `onelab_export/driver.py`, `matrices.py`, `onelab.py`,
`requirements.txt`, `ONELAB-LICENSE.txt` and the shell-handoff `launch.pro`.
The already removed `measurements.py` and `voltage_paths.jl` stay absent. Remove
`DetachedSolve`, `Launcher*`, `GetDPExecutable`, `Solver.PythonInterpreter`,
Python CLI/environment settings, process polling, JSON runtime manifests,
NumPy imports and the corresponding file-copy inventory from the exporter.
Remove `entry.txt` if no native consumer needs it; the `.pro` itself is the entry.

Replace `test/fem_onelab/` helper tests with Julia-maintained tests that invoke
the real Gmsh/GetDP executables. Remove ONELAB-specific `setup-python`, pip and
unittest steps from `.github/workflows/CI.yml`; retain native system libraries
needed by Gmsh. Update `test/extensions/fem_export.jl`, the manual export and
validation runners, `test/README.md`, `docs/src/fem.md`, the bundle README,
public export docstrings, and `THIRD_PARTY_NOTICES.md` accordingly. Remove only
the notice for the client no longer shipped; preserve Gmsh/GetDP notices.

Replace stale prototype completion claims in plan/validation documents with the
final native contract and actual evidence. Final user documentation contains
native usage only. Audit maintained FEM manual callers for deleted helper usage;
update in-scope callers without retaining compatibility functions. Unrelated
PSCAD Python integration and independent research scripts are outside this
feature and must not be broadly deleted on a keyword match.

For `overwrite=true`, use the existing ownership manifest generically: validate
its relative paths stay inside the destination, remove previously owned files
absent from the new inventory, and replace current owned files. Preserve
unrelated files and user-owned virtual environments/research outputs. Do not add
special handling named after the prototype. The new export must not reference
any retained unrelated files.

Do not rewrite the user's global `.gmshrc` as part of an exporter run. Prove the
feature with a clean ONELAB configuration and with standard GetDP configured by
absolute path while `getdp` is absent from PATH. A pre-existing global driver
entry may still need user-side removal; it must never be emitted or consulted by
this feature. Verify that distinction explicitly in GUI evidence.

**Acceptance and validation matrix**

| Gate | Required validation | Pass condition |
|---|---|---|
| A01 Public export | Problem and system dispatch; caller-owned Gmsh state; export without solver/display | Same public return contract; native files emitted without mesh/solve/UI; caller state preserved |
| A02 Bundle independence | Copy bundle to a path with spaces, outside the repo; run from another cwd with interpreter executables unavailable, environment cleared and no Julia/LineCableModels resources | Gmsh meshes and GetDP produces complete matrices and maps; process-execution evidence shows no Python/Julia or custom launcher; absence from PATH alone is insufficient proof |
| A03 Native GUI | Open fresh `study.pro`; Check, mesh, full Run, selected basis, Stop, rerun, case change; absolute standard GetDP path with solver absent from PATH | Working geometry/mesh, inputs, named outputs and fields; one normal GetDP executable authority; no generated driver entry or manual registration step |
| A04 Primitive parity | Shared mesh, 50 Hz and 10 kHz, both formulations; compare detached primitive `Z/P` with current Julia backend | Existing per-component bound `abs(a-b) <= 2e-9*abs(b) + 100eps(Float64)*maximum(abs,b)` retained separately for real and imaginary components |
| A05 Reduced parity | Repeated phase IDs, grounded terminal, retained grounds without Kron, reordered IDs, transposition on/off; compare reduced `Z/P/Y` with engine reference | All supported flag combinations agree within established component tolerances; names/order/units are exact |
| A06 Arbitrary matrix size | Complex nonsymmetric well-conditioned 1-by-1, 2-by-2 and 4-by-4 fixtures; runtime coefficient changes in one process | Native inverse agrees with independent Julia solve at existing `1e-13` absolute fixture tolerance; `max(abs(P*Y-I)) <= 1e-13`; no size ceiling or stale imaginary data |
| A07 Matrix failure and partial basis | Singular matrix, invalid normalization, failed solve, Stop; full run followed by a single basis and a failed run | No successful current full `Y` for incomplete/failed work; explicit diagnostic columns; logs/native errors retained; no retries that change the requested problem |
| A08 Native measurements | Analytic complex open/closed paths and off-mesh field; reference/metal/PML paths; basis-dependent field; maps disabled | Existing analytic `1e-13` and manufactured extraction `1e-11` checks retained; correct orientation/reference; no preparation files or hidden field-map dependency |
| A09 Line refinement | Fixed volume solution with refined measurement lines; separately refined volume mesh; native curve counts inspected | Refinement control actually changes measurement resolution; errors and sequence recorded separately from model error; no silent accuracy claim from parity alone |
| A10 Manual edits | Change named conductivity, source amplitude, a constraint and connection assignment; mesh scale/refinement; switch frequency | Changed sources read directly and remain untouched by execution; expected field/matrix response; amplitude-normalized matrices invariant; no stale mesh/case data |
| A11 Output/inspection | Primitive/reduced matrices, per-metre/total, `.pos` maps, geometry/mesh overlays and existing Makie imports | Exact units and ordering; no extra `j omega`; native Gmsh viewing and separate artifact readers work |
| A12 Repeatability and lifecycle | Repeat Run, failed Run after success, selected/full switches, reopened and relocated session | Correct completion state and outputs; no accumulation of old matrix rows, duplicated parameters or stale fields; no orphan solver process after Stop |
| A13 Removal and overwrite | Fresh bundle and overwrite of a disposable previous export; inspect owned inventory and maintained code/docs/CI | No exporter Python files/dependencies/launcher/preparation remains; stale owned assets removed; unrelated files unchanged; no compatibility path |
| A14 Shared backend preservation | Focused FEM export, native measurements, both formulations, mesh grading, workers/resume and artifact import tests | Required checks pass without changing existing tolerances or altering native headless execution/Makie contracts; unrelated baseline failures reported separately |

For A02, prefer a disposable native-runtime-only environment; otherwise use
process tracing plus deliberate rejection of interpreter execution and record
what was actually isolated. Do not invoke a Python-based Gmsh wrapper: test the
native Gmsh binary. Julia may orchestrate repository tests and independently
check numbers outside the detached child runtime; it may not supply the child's
matrix results or generated runtime inputs after export.

Use the small mixed overhead/buried two-wire study for GUI/CLI behavior and
shared-mesh parity. Use algebra fixtures for the reduction/inversion flag matrix
without multiplying expensive FEM runs. Include one physical bundled/grounded
case to validate the connection between the two stages. Execute both the
backend's GetDP 3.5 artifact and the user's installed ONELAB GetDP where available;
record versions. The initial demonstrated platform is Linux. Removing shell
handoffs is required, but do not claim tested Windows/macOS GUI behavior without
running it there.

Before updating `.linecablemodels/fem/onelab-two-bare-wires`, recover and compare
its actual geometry, materials, frequencies and terminal ordering with its source
problem. Validate a fresh native bundle first. Do not replace this user's study
with the generic toy or overwrite unrecorded manual edits. Deliver a clearly
identified native example ready to open, and record the status of the user's
existing study explicitly.

**Completion report**

Provide the native bundle path, exact verified opening/CLI commands, supported
controls, test commands/results, numerical residuals and line-refinement limits,
plus GUI evidence. Identify any unverified platform or study. Record a scoped
changed-file manifest and preserve the existing unrelated dirty work; no commit
or publication is included. Mark the feature complete only when A01–A14 are
satisfied. If any required gate remains unverified, report the feature as
incomplete, even if individual numerical experiments pass.
