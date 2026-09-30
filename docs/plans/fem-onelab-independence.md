# Headless FEM and independent Makie inspection

Status: implementation and verification complete; preexisting visual failure
and environment limits are recorded below.
Date: 2026-09-26.

The user authorized this plan, implementation, and verification as one active
goal. Computation is headless. Mesh and field inspection are separate Makie
operations. The former `ui=true` workflow is intentionally removed rather than
ported. No commit, push, publication, or unrelated scientific redesign is included.

## Required outcome

1. Owned Julia geometry constructs Gmsh primitives and meshes them without
   ONELAB or FLTK. The existing GetDP subprocess path receives explicit files
   and arguments and returns validated output files.
2. A saved mesh can be loaded and inspected alone or underneath the existing
   geometry preview without meshing, resolving GetDP, or solving a problem.
3. A saved field map can be loaded and plotted through the existing PlotBuilder
   shell without GetDP, ONELAB, FLTK, or a live solver session.
4. Mesh and field data are detached Julia values. Gmsh session ownership,
   physical coordinates, element identities, units, and field discontinuities
   remain explicit. Makie retains presentation ownership.
5. Current numerical behavior, worker isolation, partial recovery, cancellation,
   artifacts, and unrelated working-tree changes are preserved.

Independence means no owned ONELAB/FLTK execution or parameter contract. Gmsh
remains the native meshing and file-reading library; no custom Gmsh build is
required merely because that library also provides optional ONELAB functionality.

## Starting evidence

The tracked HEAD is `94e4b058706031e937954802cec3e2921c4dd751`, but HEAD is not the
implementation baseline: this checkout contains ongoing FEM/PML changes.
The initial tracked and nonignored untracked files, SHA-256 manifest, and Git
status were captured in `/tmp/lcm-onelab-baseline-lg30kjop` (732 files).
That snapshot is preservation evidence and is not a source of replacement files
for unrelated user changes.

Verified source boundaries:

- `ext/LineCableModelsGmshExt/geometry.jl:_build_geometry!` already constructs
  native Gmsh primitives.
- `mesh.jl` already invokes `gmsh.model.mesh.generate(2)` and writes MSH 4.1.
- `workers.jl:_getdp_command` and `_start_worker!` already run GetDP directly;
  there is no `-onelab` worker transport to replace.
- `onelab.jl`, Gmsh-session snapshots, the headless publication calls, and
  `_ui_solve!` retain the unwanted coupling.
- `getdp/model.pro` attaches ONELAB options to its Physics constant.
- `results.jl:_merge_maps!` and the UI mesh-loading branch exist for FLTK display.
- `PlotBuilder.preview`, `plot`, `import_data`, the Makie shell, and `UIPlot`
  are existing public owners to extend.
- Retained native POS files contain element-local scalar/vector values with
  separate real and imaginary output steps. Saved examples also carry physical
  units, basis/frequency labels, material tags, and PML masks.

References: [GetDP CLI](https://getdp.info/doc/texinfo/getdp.html#Running-GetDP)
and [Gmsh API](https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-application-programming-interface).

## Ownership and public interaction

Use existing owners and generics; do not introduce a transport service, a second
solver framework, a second plotting shell, or a registry that repeats dispatch.

- Engine owns small passive FEM mesh/field values and their scientific invariants.
  They are not line-parameter matrices and do not require a fake calculation or
  reconstruction of a LineParametersProblem.
- ImportExport owns the `import_data` public format boundary. The Gmsh extension
  implements native mesh/POS decoding through that boundary.
- The Gmsh extension owns native session initialization, locking, scratch models,
  extraction, and cleanup. Importing a file must not call `_getdp_selection`,
  `_fem_input_record`, mesh generation, or computation-run preparation.
- The Makie extension consumes detached values and draws native mesh/line plots.
  File conveniences acquire through `import_data` and use the same rendering
  methods as already-loaded values. There is no private extension-to-extension API.
- Existing `preview(system; ...)` gains `mesh=...`; `plot` accepts saved mesh and
  field files and imported values. Existing line-result plotting stays intact.

Target interactions:

```julia
mesh = import_data(:msh, "frequency_0001.msh")
preview(system; mesh)
plot(mesh)
plot("frequency_0001.msh")

field = import_data(:pos, "e_f0001_b0001.pos")
plot(field; component=2, part=:real)
plot("e_f0001_b0001.pos"; component=2, part=:real)
```

Native file decoding must retain raw step metadata. An external POS view must
not become complex merely because it contains two steps. Owned harmonic output
is recognized from its recorded convention; otherwise the reader requires an
explicit step/representation selection. Preserve multiple views explicitly.

## Execution sequence and gates

### 1. Establish the effective baseline

- Record starting source/dependency identities and relevant test outcomes.
- Use small existing transport, workers, resume, options, and geometry checks.
- Capture native matrix and field output from representative aerial/buried
  cases, both physics modes, and more than one frequency/source. Keep original
  output and comparison settings. Existing retained artifacts can supplement,
  but cannot replace, a known-current preservation run.
- Reproduce any encountered failure on the saved starting implementation before
  attributing it to this change. Do not change unrelated PML equations, accuracy
  policy, fixtures, or tolerances to obtain a pass.

Gate: the starting implementation and each comparison's limits are recorded.

### 2. Remove ONELAB and interactive computation

- Remove `onelab.jl`, its includes and fingerprints, publication and completion
  calls, and ONELAB session state.
- Remove `_ui_solve!`, FLTK polling, UI-only mesh/map loading, and `ui` as a
  supported computation option. Unsupported options must remain explicit errors.
- Keep `plot_field_maps` as the existing artifact-generation control. Keep
  `keep_run_directory` behavior explicit; inspection cannot resurrect deleted files.
- Replace ONELAB Physics metadata with ordinary CLI-controlled constants.
- Preserve the direct worker command, frequency batching, factorization reuse,
  file validation, progress logging, process cleanup, and checkpoint ownership.
- Preserve interruption/cancellation using the existing execution owner, without
  retaining a GUI pump solely for testing or adding a new callback framework.
- Update all live callers, docs, and tests. Historical execution records remain
  historical evidence; do not rewrite them to describe the new implementation.

Gate: native headless solves and worker lifecycle checks pass without owned
ONELAB/FLTK calls. The changed source fingerprints reject incompatible resume;
historical files remain readable for inspection.

### 3. Introduce detached native-file imports

- Retain coordinates in Float64, sparse/nonconsecutive node tags, element blocks
  and order, entity membership, and all physical-group memberships.
- Read in uniquely named scratch models/views using the same Gmsh ownership and
  lock discipline. Restore caller-owned state on success and failure; finalize
  only an importer-owned session. Release resources before creating Makie objects.
- Use native Gmsh extraction APIs instead of a handwritten MSH/POS text parser.
- Preserve list-based POS element-side samples, components, views, and native
  output-step times. These are the field maps emitted by the existing solver.
  Model-based views require explicit export to list-based POS; they must not be
  guessed from an exception or silently converted into another association.
- Separate file validity from solver admissibility. An ordinary valid planar
  mesh need not contain LCM-specific physical tags to be displayed.
- Cover the first-order planar triangles emitted by this stack, lower-dimensional
  boundaries, and supported quadrangles. Retain unsupported blocks and reject an
  unsupported rendering request precisely rather than dropping data or flattening
  a nonplanar/volume mesh. Higher-order display must be explicit tessellation,
  never silent removal of midside nodes.

Gate: fixtures check sparse tags, physical groups, scalar/vector data, complex
step interpretation, discontinuities, malformed files, and session cleanup.
Import succeeds with an invalid GetDP override and no display.

### 4. Render saved meshes in the existing shell

- Draw standalone mesh edges with native Makie primitives and equal aspect.
- Extend geometry preview with a mesh layer below existing geometry, with mesh
  visibility/appearance controls and retained geometry/material behavior.
- Map solver horizontal x / vertical y to preview horizontal y / vertical z
  explicitly; apply the same convention to vector components and labels.
- Preserve preview limits, reset behavior, controls, guides, exports, and live
  `UIPlot` handles. Offer full-mesh limits through ordinary axis controls.
- Batch segments; do not create one Makie plot object per element.

Gate: rendered mesh and geometry align at known coordinates; preview without a
mesh retains existing behavior; detached meshes render after Gmsh finalization.

### 5. Render saved field maps

- Use the same PlotBuilder shell, native colorbar/axis controls, and export path.
- Support real/imaginary parts, magnitude, and component phase; vector phase
  requires a component. Support explicit vector display with documented selection.
- Preserve element-local values: duplicate render vertices where needed instead
  of welding discontinuous samples or averaging across material interfaces.
- Retain units, source normalization, basis/frequency identity, scaled quasi-fw
  quantities, material tags, and PML-continuation labels. Do not infer scientific
  meaning solely from filenames or hide PML as if it were physical material.
- Support geometry outlines and optional mesh overlays without opaque geometry
  hiding field values. Load only requested views/maps; a frequency change uses
  its own retained coordinates/connectivity.
- Apply logarithmic color limits only to valid positive support; report excluded
  or undefined samples without modifying retained scientific values.

Gate: known scalar/vector fixtures verify values, components, labels and jumps;
real retained maps from both physics modes render with correct units and overlays.

### 6. Complete validation and documentation

- Re-run affected numerical/lifecycle checks and current-baseline comparisons.
- Run ordinary, core-boundary, quality, and relevant visual checks; use existing
  backend activation paths for Cairo/GL/WGL and report which interaction was
  actually exercised.
- Inspect rendered images and exports, not only object types and file existence.
- Update FEM, plotting, import API documentation, examples, and test commands.
- Run docstrings/docs checks and task-scope diff/format checks. Preserve unrelated
  formatting and edits; do not run an unbounded repository rewrite.
- If claiming performance, measure first-use and warmed execution separately,
  distinguishing mesh creation, worker startup, assembly/solve, map output,
  import, and rendering. ONELAB removal alone establishes no speedup.
- Record commands, outcomes, artifacts, baseline failures, and unverified scope.

Gate: all required implemented workflows are exercised end to end. Never label
the goal complete while required implementation or verification remains.

## Validation commands

Run from the checkout root with Julia 1.12 and the existing manifests. First
inspect the selected items with `--list` when selection scope is uncertain.
The repository runner combines filename alternatives with tag filtering.

```sh
julia --project=test test/runtests.jl extensions/fem_transport.jl extensions/fem_workers.jl extensions/fem_resume.jl unit/engine/options.jl
DISPLAY= julia --project=test --compiled-modules=no test/runtests.jl tag:fem_numerical
julia --project=test test/runtests.jl
julia --project=test/core test/runtests.jl tag:core_only
julia --project=test test/runtests.jl tag:quality
LINECABLEMODELS_TEST_PLOTTING=true julia --project=test/visual --compiled-modules=no test/runtests.jl tag:visual
julia --project=docs docs/doctest.jl
julia --project=docs docs/make.jl
git diff --check
```

Add focused selections for new import/render tests. Use an isolated writable
depot before diagnosing managed-host cache permissions as a code failure.
Prepare changed dependency environments once before testing; do not overlap cold
preparations. Do not launch PSCAD native tests or the existing long PML research
campaign as part of this task. Numerical preservation is implementation evidence,
not a new claim of physical convergence or accuracy.

## Execution record

- [x] User selected separate headless computation and Makie inspection.
- [x] User authorized plan, implementation, and verification in goal-pursuit mode.
- [x] Active goal created; no token budget was requested.
- [x] Initial dirty source snapshot captured.
- [x] Baseline checks and numerical preservation artifacts captured.
- [x] ONELAB/FLTK removed from computation.
- [x] Detached native-file imports implemented and verified.
- [x] Mesh inspection and preview overlays implemented and verified.
- [x] Field inspection implemented and verified.
- [x] Required composed workflows, documentation, and final audit verified.

### Final implementation

The owned computation path contains no ONELAB or FLTK calls. Gmsh constructs
and meshes the existing geometry headlessly; GetDP receives the existing explicit
command arguments and files. The `ui` option and its display-only branches are
removed. Worker cancellation, retained runs, source fingerprints, reuse, and
partial recovery remain in the existing execution owner.

`import_data(:msh, path)` returns a detached `FEMMesh`; `import_data(:pos, path)`
returns detached `FEMFieldMap` values. Native Gmsh APIs own decoding, including
binary files, sparse tags, physical memberships, element-side samples, view
selection, and output-step times. Only explicitly marked harmonic output or
`representation=:complex` combines two steps into a phasor.

The existing Makie shell now provides `plot(mesh)`, `plot(field)`, filename
conveniences, and `preview(system; mesh)`. Field plots retain discontinuities and
metadata and support components, real/imaginary parts, magnitude, phase, native
color scales, geometry/mesh overlays, vector arrows, and SVG export. Arrow
normalization and length scaling use native Makie `arrow_attributes`. The
rendering scope is first-order planar meshes and triangular/quadrangular list
POS fields; unsupported geometry is rejected explicitly without altering imported
data. See `docs/src/fem.md` for complete examples and representation conventions.

### Measured checks

All commands used Julia 1.12.7 and the writable depot overlay
`JULIA_DEPOT_PATH=/tmp/lcm-onelab-depot:/home/amartins/.julia`. Results below are
from the final source or the explicitly identified starting snapshot. Counts
from separate invocations are not presented as one clean full-suite run.

| Check | Result | Evidence |
| --- | --- | --- |
| Starting transport/workers/resume/options | 276/276 passed | `/tmp/lcm-onelab-baseline-tests.log` |
| Final import/transport/workers/resume/options | 323/323 passed; includes 48 native import assertions | `/tmp/lcm-onelab-targeted-final.log` |
| Core without optional Gmsh | 35/35 passed | `/tmp/lcm-onelab-core-final.log` |
| Full native FEM selection | 191 passed; two invocation errors resolved by the reruns below | `/tmp/lcm-onelab-native.log` |
| Fresh GetDP artifact, compiled and uncompiled modes | 4/4 passed with network access | `/tmp/lcm-onelab-artifact-network.log` |
| Complete multi-frequency/reuse/recovery case | 100/100 passed with source fixed for the duration | `/tmp/lcm-onelab-native-reuse-final.log` |
| Ordinary suite | 23,716 passed, one sandbox X11 error; activation rerun below passed | `/tmp/lcm-onelab-ordinary.log` |
| Quality | 2,389 passed, 18 new text-display failures corrected; focused rerun below passed | `/tmp/lcm-onelab-quality.log` |
| Owner-local text display after correction | 512/512 passed; starting source 494/494 | `/tmp/lcm-onelab-text-display-final.log` |
| Full workspace visual selection | 3,611 passed, one preexisting undefined-phase/log-axis error reproduced on the starting snapshot | `/tmp/lcm-onelab-workspace-visual.log`, `/tmp/lcm-onelab-baseline-resolution.log` |
| Affected pinned preview/material/control/FEM selection | 264/264 passed | `/tmp/lcm-onelab-preview-final.log` |
| Final FEM inspection including arrow attributes | 27/27 passed | `/tmp/lcm-onelab-arrows-final.log` |
| Existing Cairo entry point | Starting and final sources each 26/26, with both workspace and pinned Makie versions | `/tmp/lcm-onelab-baseline-cairo.log`, `/tmp/lcm-onelab-cairo-final.log`, `/tmp/lcm-onelab-baseline-pinned-cairo.log`, `/tmp/lcm-onelab-current-pinned-cairo.log` |
| Cairo/GL activation | 16/16 passed with access to the existing display | `/tmp/lcm-onelab-gl-activation.log` |
| WGL activation | 8/8 passed; detached mesh/field figures also constructed | `/tmp/lcm-onelab-wgl.log` |
| Doctests | Passed | `/tmp/lcm-onelab-doctest.log` |
| Final documentation build | Passed with existing citation/size warnings | `/tmp/lcm-onelab-docs-final.log` |

The native preservation script `/tmp/lcm-onelab-preservation.jl` exercised two
terminals (one aerial and one buried), 50 Hz and 10 kHz, and both physics modes
using the same fixed eight-layer PML discretization. Z, Y, and primitive matrices
compare exactly: maximum absolute differences are zero. All 56 quasi-TEM and
64 quasi-fw map payloads compare byte-for-byte after the first view-label line.
The intentional label change records `phasor=real,imag`; sampled values did not
change. Before/after artifacts are `/tmp/lcm-onelab-preservation-before` and
`/tmp/lcm-onelab-preservation-after`; comparison output is
`/tmp/lcm-onelab-preservation-comparison.log`. This is preservation evidence,
not a new physical-convergence or performance claim.

Native saved-file inspection rendered retained maps from both formulations,
meshes, and geometry overlays with an invalid GetDP override and no display.
Final Cairo images and SVG exports are in `/tmp/lcm-onelab-visual-final`.
Both formulations were visually inspected at original resolution for coordinate
alignment, field support, units, scaled quasi-fw labels, and PML qualifications.
The GL backend also rendered mesh and discontinuous field images; evidence is in
`/tmp/lcm-onelab-gl-spatial.log` and `/tmp/lcm-onelab-visual/gl-mesh.png` /
`gl-field.png`. WGL browser interaction was not exercised.

### Invocation limits and scope audit

- The first fresh-artifact invocation could not resolve download hosts inside
  the sandbox. The network-enabled rerun passed in both Julia compilation modes.
- One native parent process retained the source identity from before a local
  mesh-variable rename, while its recovery subprocess saw the edited source.
  Recovery correctly rejected the mismatch. The entire case passed after
  holding the source fixed; no numerical tolerance or acceptance rule changed.
- The ordinary suite's only error was GLFW opening the restricted X11 display.
  The same Cairo/GL activation item passed with display access.
- The full visual suite's only error is construction of an all-undefined phase
  plot in `test/integration/observable_resolution_plots.jl:31`: Makie rejects
  x limits `(0.0, 10.0)` for `log10` at `shell.jl:332`. The identical item on
  the initial source snapshot produced the same error after seven passing
  assertions. Both the test and shell are byte-identical to the starting tree.
  This preexisting issue is recorded without changing unrelated plotting code.
- Quality failures concerned display ownership of the four new detached data
  types. Existing `TextDisplay.@showfields` definitions resolved them. All other
  quality items passed in the full invocation, including explicit imports.
- The workspace uses Makie 0.24.15 / CairoMakie 0.15.15; `test/visual` pins
  Makie 0.24.13 / CairoMakie 0.15.13. Full pinned visual invocations terminated
  twice with signal 15 before a test result. No assertion outcome is inferred
  from them; the pinned focused selections and full workspace run are reported
  separately above.
- All 32 changed Julia files parsed without executing manual studies.
  `git diff --check` passed. The explicit task manifest is
  `/tmp/lcm-onelab-task-paths.json`: 45 changed paths against the initial dirty
  tree, with no changes outside that manifest.
- Geometry, voltage paths, PML coefficients, materials, integration, and Jacobian
  sources remain byte-identical to the initial snapshot. The preexisting Makie
  shell and unrelated user edits are also unchanged. No commit, push, or
  publication was performed.
