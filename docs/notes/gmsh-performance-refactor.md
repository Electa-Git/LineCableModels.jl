# Gmsh/GetDP scan performance: diagnosis and refactor design

2026-09-08. This records the initial investigation and executable diagnostic
patch. Production integration and full-campaign validation subsequently
completed; see the [acceptance report](gmsh-performance-validation.md).
The investigation left production code, meshes, retained results, and
campaign state untouched. The supplied documents were read in full
as technical references; their embedded examples and instructions were not
treated as additional user requests.

The [implementation plan](gmsh-performance-implementation-plan.md) orders the
production changes and defines the acceptance gate for each stage through a
complete campaign benchmark.

## Recommendation and measured result

Make **one GetDP process own one frequency, its mesh, and all requested terminal
excitations**. Assemble and factor the coupled system once. For subsequent
terminals, update the prescribed currents, regenerate only the RHS, and solve
using the existing factors. Julia continues to own the frequency scan and
precomputed meshes.

This is supported by an experiment on copies of the completed 18 kV campaign:

| Model at 316.227766 Hz | Separate processes, current operations | One process, current operations | One process, RHS/factor reuse |
|---|---:|---:|---:|
| Detailed, 189,505 GetDP DOFs | 91.13 s | Not measured | **15.60 s** |
| Homogenized, 41,762 GetDP DOFs | 24.13 s | 20.66 s | **4.17 s** |
| Processes per frequency | 9 | 1 | 1 |
| Symbolic / numeric LU factorizations | 9 / 9 | 9 / 9 | **1 / 1** |
| Linear solves | 9 | 9 | 9 |

All timings cover nine excitations and both primitive matrices. Merely moving
the loop saves comparatively little. Reusing assembly and factors produced
about **5.8×** acceleration in both controlled comparisons.

The batch probe also passed at 0.1 Hz and 1 MHz for both models, and at the
middle frequency for the second detailed and homogenized material profiles:
eight complete 9-by-9 Z/P pairs, 1,296 complex entries in total. Every emitted
complex entry matched the corresponding retained result exactly after parsing
the TSV values. Every batch recorded one symbolic factorization, one numeric
factorization, and nine solves. The largest reported `GetResidual` norm across
these batches was approximately 2.31e-10; RHS norms were 1.

These are single-sample solver measurements on this machine, with existing
meshes, inherited thread settings, and field maps disabled as in the campaign.
They include startup, preprocessing, residual checks, and raw output; they
exclude mesh copying, meshing, Julia, and ONELAB. They establish an optimization
opportunity, not a measured speedup for a complete 101-frequency campaign.
They establish output equivalence on these samples, not independent physical
validation or equivalence of every field value.

Machine-readable results, executable identity, source and mesh hashes, memory
observations, and PETSc event counts are in
[gmsh-performance-evidence.json](gmsh-performance-evidence.json).

## Where the time goes

In [getdp.jl](../../ext/LineCableModelsGmshExt/getdp.jl),
`_run_getdp_unlocked!` nests a terminal loop inside the frequency loop and calls
`gmsh.onelab.run` for every pair. `FEMSolveBasis` in
[quasi_tem.pro](../../ext/LineCableModelsGmshExt/getdp/quasi_tem.pro) performs
`UpdateConstraint`, `Generate`, and `Solve` on every invocation. Precomputing
the meshes already avoids meshing each excitation, but each GetDP process
still rereads its mesh, constructs its DOFs, assembles the same matrix, and
factors it again.

The `dielectric-baseline-fem-common-band` campaign records 17,097.60 seconds
(4.75 hours) for the two detailed 18 kV profiles and 4,989.23 seconds for the
two homogenized profiles. The detailed runs `run-Pbj2uz` and `run-orkY1H`, and
homogenized runs `run-YiqIRh` and `run-aZVWtq`, each record 909 launches:
101 frequencies times nine terminals. These are terminal **columns** in the
current Julia convention, with response terminal indexing rows. Some legacy
Python code uses the opposite orientation.

The retained `logs/getdp.log` files are empty. Gmsh logger capture at the
campaign's zero verbosity supplies no phase breakdown. Live invocation counts
also lag execution because state is persisted at lifecycle transitions. Raw
file timestamps can indicate throughput, but cannot identify solver phases.

A profiled first solve of the detailed middle-frequency mesh showed:

| Stage | Approximate wall time |
|---|---:|
| Input parsing, mesh loading, preprocessing and initialization | 1.39 s |
| Updating constraints | 0.58 s |
| Full assembly | 4.45 s |
| Initial linear solve, including factorization/setup | 3.43 s |
| Each later excitation with group RHS assembly and factor reuse | 0.66–0.67 s |

For a representative later excitation, constraint updating consumed about
0.57 s, RHS generation 0.03 s, and the reused solve about 0.06 s. Scalar
`OnRegion` Z/P output was negligible. Optional whole-mesh field output was
disabled and cannot explain this campaign's cost.

The actual executable is GetDP `3.6.0-git-1cf7fa06`, built 2024-08-31, with
complex PETSc 3.14.4. Its observed configuration is `KSP=preonly`, `PC=lu`,
factor package `mumps`. Thus factorization means direct sparse LU, not an
unobserved iterative convergence problem. The maintained four `.pro` files
match the four sampled runs' retained solver snapshots byte for byte.

Parsing is secondary: detailed `model_data.pro` contains about 1.63 MB of
whole-scan material arrays for 366 material regions, and parsing/loading the
problem before mesh loading took about 0.20 s in the probe. Eliminating that
repetition helps, but assembly and factorization dominate.

Julia also repeatedly parses and appends column files, recounts both growing
aggregate tables after every column, and rewrites cumulative logs. These are
avoidable costs that grow quadratically with job count; they are secondary to
the measured solver work.

## What must remain invariant

For a fixed frequency and mesh, the implemented problem has the form

```math
A(f,\mathcal M,\sigma,\epsilon,\mu,\gamma)\,x_k=b_k.
```

The coupled `a`, `ur`, and `phi` equations are linear with fixed propagation
constant. Selecting terminal k only changes prescribed global currents through
`$FEM_I~{t}`. It does not change material coefficients, function spaces, support
groups, or the constrained DOF set. `GlobalTerm [Dof{I}, {U}]` supplies the
terminal current source. Consequently all columns share the same A.

The refactor must preserve:

- The coupled formulation, its conductive and displacement terms, gamma and
  gamma-squared terms, integration rules, and infinity transformation.
- All terminal current constraints: the driven terminal gets one ampere and
  all other terminal currents are reset to zero on every excitation. They are
  not changed to zero-voltage constraints.
- Magnetic outer Dirichlet conditions and the existing scalar-potential
  reference on the outer earth boundary. The outer air boundary retains its
  natural condition. The legacy magnetic/electric variants differ here.
- Resolved frequency-dependent admittivity, including dielectric losses and
  temperature corrections. GetDP must not add dielectric loss a second time.
- Current extraction: `Z = -U / UnitSource`,
  `P = Phi / (gamma * UnitSource)`, indexed `[response, basis, frequency]`.
  A single excitation yields a column of **both** matrices.
- Existing reduction order and the condition-checked solution of `P * Y = I`
  in Julia. No extra `jω`, transpose, symmetrization, or regularization.

Factor reuse is invalidated by any change to the mesh or DOF layout, frequency,
material coefficients, gamma, shell radii/Jacobian, constraint support, or
operator configuration. A changed RHS alone is the permitted case. With the
observed direct `preonly` solver, stale LU factors after an operator change
would solve the wrong problem. PETSc's reuse setting deliberately permits
keeping a preconditioner even when matrix values change; the application must
enforce the invariant. See the official
[PETSc reuse contract](https://petsc.org/release/manualpages/KSP/KSPSetReusePreconditioner/).

## GetDP operation sequence and DSL contracts

The tested sequence is:

```text
InitSolution[Sys_FEM];
// Existing frequency, index, SetFrequency and output-time initialization.
Evaluate[$FirstSolve = 1];
For basis In {1:NumTerminals}
  Evaluate[$FEMBasisTerminal = basis];
  Call FEMSetBasisCurrent;
  UpdateConstraint[Sys_FEM];
  Test[$FirstSolve]{
    Generate[Sys_FEM];
    Solve[Sys_FEM];
    Evaluate[$FirstSolve = 0];
  }{
    GenerateRHSGroup[Sys_FEM, Terminals];
    SolveAgain[Sys_FEM];
  }
  PostOperation[FEMAppendRaw];
EndFor
```

`For` here expands operation statements during parsing; `Evaluate` changes
runtime variables when those operations execute. A parser assignment to a
constant is not a substitute for setting the runtime current variables.
`UpdateConstraint` is necessary because assigned constraint values otherwise
retain preprocessing values. The first solve is the **first requested**
terminal, which need not be terminal 1 when resuming.

`GenerateRHS` plus `SolveAgain` is demonstrated in the supplied
`models-getdp/ElectromagneticScattering/scattererTmatrix.pro` at lines 176–191.
`models-getdp/GetDDM/SchwarzMacros.pro` at lines 132–170 also demonstrates
`GenerateRHSGroup` and scoped constraint updates. Both RHS operations parse
and run on the installed executable. The supplied GetDP 3.5 manual does not
list these RHS operations, so a production change must declare or probe its
supported executable capabilities rather than assume every GetDP version
supports the same grammar.

The group-restricted RHS operation is appropriate here because the changing
source is the terminal `GlobalTerm`; all fixed prescribed boundary values are
zero. It must be reconsidered if nonzero boundary sources, distributed sources,
or changing support groups are added. Full `GenerateRHS` is the conservative
alternative when RHS support cannot be restricted. Switching to
`UpdateConstraint[Sys_FEM, Terminals, Assign]` also preserved the sampled
results, but gave no material speed improvement: 15.40 s versus 15.60 s for
the detailed case. Do not count on that syntax to eliminate a mesh traversal.

Output handling needs explicit changes alongside the numerical loop:

- `AppendToExistingFile 1` has within-PostOperation behavior; it does not mean
  append safely across successive calls. The diagnostic uses value **2** with
  fresh, isolated files. GetDP then emits additional blank lines, which the
  production parser already ignores.
- Production should keep independently validated column checkpoints. Use
  basis-specific PostOperations expanded by the maintained `.pro` template,
  or verified runtime filename suffixes, with separate attempt directories.
  A single repeatedly appended frequency file is insufficient for the
  existing resume contract without explicit framing and validation.
- Field-map filenames and labels currently use parser `BasisTerminal`.
  Moving only the solve loop would reuse names and overwrite maps. Move map
  emission inside the loop and bind its names to the actual basis. Preserve
  all nine optional quantities and `LastTimeStepOnly` behavior.
- Keep `SetFrequency` responsible for physics and `SetTime[FrequencyIndex]`
  responsible for output identity. Do not replace the latter with Hz or
  accumulate an unnecessary solution history.

The runnable [diagnostic patch](gmsh-terminal-batch-probe.patch) demonstrates
the numerical sequence and append behavior. It intentionally omits production
resume integration and supports **field maps off only**.

## Why one process for the entire frequency scan is a later option

GetDP can describe several systems with distinct `NameOfMesh` declarations;
the supplied GetDDM examples do this for separate domains. A generated
multi-system resolution could sequence an entire scan inside one process.
That is a possible design, but it requires explicit systems/meshes and careful
memory lifecycle management. The existing system cannot be remeshed merely by
calling `SetFrequency` or `GmshRead`.

There is also a parser/runtime dependency: `materials.pro` indexes material
arrays using parser `FrequencyIndex`, and the shell Jacobian uses the selected
radii. A runtime loop that changes only frequency would retain stale material
bindings or domain parameters. Each mesh needs its matching system and data.

With genuinely different frequency meshes, the matrix dimensions/numbering
and values can change. Numerical factors therefore cannot be shared across
frequencies. Even on an identical mesh, changing omega or dispersive material
values requires new numerical factors. Reusing graph ordering or element
geometry is a separate optimization with narrower invariants.

One process per frequency already changes the campaign from 909 launches,
assemblies and factorizations to **101 of each**, while retaining 909 solves.
That captures the demonstrated reuse with bounded memory and a natural retry
boundary. A full-scan process might subsequently save some startup/parsing,
but it does not provide another ninefold factorization reduction.

`GenerateSeparate`/`Update` examples concern linear constant-coefficient time
integration; they do not automatically cache this remeshed harmonic scan.
Likewise, GetDP's `-cache` option is documented for network computations, not
general LU caching. ONELAB parameter transport is not a factor cache.

## Mesh findings and design continuity

The original characteristic-length and field approach is already recognizable
in the Julia implementation. `geometry.jl` shares interface curves and points,
builds conformal surfaces, and checks material adjacency. `model.jl` computes
local solid/annulus/strand/foil sizes. `_configure_mesh!` in `mesh.jl` enables
point sizes, disables extension from boundaries, and combines cable-exterior
`Distance`/`Threshold` fields using `Min`, sampling 100 and a growth parameter
of 1.2. This follows the useful legacy `.geo` design rather than assigning the
smallest internal feature size to the entire earth domain.

| Frequency | Detailed nodes | Detailed cable triangles | Detailed exterior triangles | Homogenized nodes |
|---|---:|---:|---:|---:|
| 0.1 Hz | 117,804 | 221,612 | 13,762 | 33,184 |
| 316.227766 Hz | 116,971 | 221,612 | 12,096 | 32,314 |
| 1 MHz | 116,123 | 221,614 | 10,398 | 31,497 |

Exterior includes finite air/earth and both infinity-shell halves. These are
actual 2D element counts, not the MSH count that also includes boundary
elements. At the middle frequency, about 94.8% of detailed triangles are
inside the cables. The finite radius changes from approximately 15.9 km at
0.1 Hz to 5.03 m at 1 MHz, while total node count changes by only about 1.4%.
The detailed case has 366 material regions versus 18 in the homogenized case.

The mesh may still have opportunities for improvement, particularly around
the detailed internal features. Element counts alone cannot establish that
those elements are unnecessary. Preserve the meshes for the orchestration
refactor. A later mesh study should report counts and quality per material and
compare Z/P/Y convergence when changing local resolution. Conductor interior
sizes currently follow geometry, rather than conductor skin depth, so a
coarsening proposal also needs high-frequency skin/proximity validation.
Homogenizing strands, removing thin regions, or replacing boundaries with
impedance conditions is a separate modeling decision.

## Refactor boundaries and delivery order

| Change | Main locations | Reviewable outcome |
|---|---|---|
| Frequency batch protocol | `getdp.jl`, maintained `quasi_tem.pro` | One process per mesh; list of requested bases; first full solve then RHS/reused solves; existing equations unchanged |
| Completion/resume | `getdp.jl`, `results.jl`, run state | Independently validated column artifacts; ordered aggregation; failed attempts cannot declare a complete scan |
| Solver diagnostics | `getdp.jl`, run metadata | Per-frequency stdout/stderr, exit status, elapsed stages, mesh/DOF counts, invocation and completed-column counts persisted separately |
| Headless transport | `getdp.jl`, `model.pro`, `compute.jl` | Julia launches a `Cmd` with immutable paths and explicit parameters; ONELAB remains available for the UI |
| Maps and public documentation | `quasi_tem.pro`, `docs/src/fem.md`, extension tests | Correct frequency/basis names and labels; documented invocation semantics |

First integrate batching and checkpoint handling using the current transport;
this isolates the measured optimization. Then separate headless execution from
ONELAB. The probes already demonstrate that the physics runs standalone when
`ModelDataPath` is supplied directly. Keep data-only input serialization:
Julia supplies parameters and a requested-basis list; the maintained `.pro`
owns loops, equations and outputs. There is no need to restore a Python driver
or emit a second solver program into `model_data.pro`.

For a frequency attempt, validate and skip existing compatible columns,
submit only missing bases, and factor on the first missing basis. Write new
outputs into an exclusive attempt directory. Adopt a column only after its Z
and P contain exactly one finite value for every response, with the expected
frequency/index/basis, and required maps are complete. Completion framing must
distinguish a fully written column from a killed write. On failure, preserve
valid columns, diagnostics and partial artifacts; a later attempt rebuilds
the factorization and solves only missing columns. Aggregate validated columns
once in deterministic order, then publish the existing scan completion marker
and checksums. Persist process launches separately from completed columns.

Version the execution protocol and record solver/options/source identities.
Keep the current incompatibility checks for resume; do not reinterpret or
overwrite already completed campaigns under the new implementation. Keep a
single-column/full-assembly diagnostic path for comparisons.

For standalone execution, use Julia `Cmd` argument arrays and direct file
capture, avoiding shell quoting and dependency on shared ONELAB values. The
existing run-specific sockets and session restoration remain relevant to UI
execution. ONELAB may return existing parameter values unless the correct
read-only/definition rules are used; making the numeric job self-contained
removes that ambiguity. Capture actual GetDP output rather than assuming
`gmsh.logger` contains the external solver's timings. Distinguish a solver
failure from Gmsh's known post-output/automatic-merge failures.

Only after batching is stable, consider a bounded pool of independent
frequency processes. Gmsh/ONELAB state is shared within a Julia process, so
threading `gmsh.onelab.run` is not the concurrency design. Mesh preparation
stays under its current session ownership. Limit solver concurrency using
measured memory and solver thread settings: the middle-frequency probes
reported approximately 895 MB for a detailed worker and 307 MB for a
homogenized worker. Larger cables need their own measurements. Parallelizing
terminals would replicate precisely the matrices and factors we can share.

The next residual cost after reuse is constraint updating. Investigating a
more direct RHS representation could be worthwhile, but changes how sources
are assembled and needs a separate equivalence study. Simply scoping
`UpdateConstraint` did not resolve it. Multi-RHS solver bindings, exported
matrices, iterative solver tuning, frequency-mesh buckets, and separate
near/far submeshes should follow evidence of a remaining need.

## Frequency parallelism: tested design

Two isolated GetDP processes can safely solve different frequency batches at
the same time. Julia only needs to supervise external processes; multiple
Julia execution threads are not required for the solver processes to use
different CPU cores. Use a bounded queue of frequency jobs, supervised through
Julia tasks and `Cmd`, with a single coordinator handling shared run state.
Julia documents external process and per-command environment handling in its
[process manual](https://docs.julialang.org/en/v1/manual/running-external-programs/).

A further diagnostic ran the detailed 18 kV model at 0.1 Hz and 1 MHz, with
all nine terminals batched at each frequency. Each process had
`OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` in its own environment:

| Execution of the same two batches | Total wall time |
|---|---:|
| Sequential | 29.51 s |
| Two concurrent processes | **14.90 s** |

This is **1.98×** acceleration on top of terminal factor reuse. All 324 parsed
complex Z/P entries were identical between sequential and concurrent runs.
Relative Frobenius differences from the original campaign were at most
1.88e-12 after changing the inherited thread configuration. Each of the four
process executions recorded one symbolic factorization, one numeric
factorization and nine solves. These measurements are in the accompanying
evidence JSON. Higher worker counts, field maps, and production recovery
under concurrent failures have not yet been tested.

The existing runner cannot simply be wrapped in `Threads.@threads`: it shares
`Solver.SocketName`, Gmsh logger state, ONELAB parameters, run counters and
aggregate files. The worker boundary should have these properties:

- Workers launch standalone GetDP, with a fixed mesh and immutable material
  data, and never call the shared Gmsh/ONELAB API.
- Every frequency attempt owns its working directory, `-name` prefix, raw
  output, maps and logs. Workers do not append to the same aggregate table.
- One coordinator validates completed columns, persists run metadata and
  counters, and assembles the final frequency-indexed arrays. Worker finish
  order never becomes matrix frequency order. Gmsh/UI updates stay with the
  session owner.
- The queue limits active solver processes. Begin with the tested two-worker
  configuration, then benchmark four and eight against CPU and memory limits.
  The previous middle-frequency detailed probes reported about 0.9–1.0 GB per
  process; reserve additional memory for Julia/Gmsh and larger cases.
- Set thread limits on each child command, not by concurrently mutating Julia's
  global `ENV`. Multiple frequency workers each starting a large BLAS/OpenMP
  team would oversubscribe the host. Keep one solver thread per worker for the
  initial comparison, then measure other worker/thread combinations. The
  applicable BLAS controls depend on its build; see
  [OpenBLAS runtime variables](https://www.openmathlib.org/OpenBLAS/docs/runtime_variables/).
- Cancellation and failure handling explicitly reap launched processes and
  retain completed column checkpoints. A failed worker must not leave an
  untracked solver writing into an attempt that is being retried.

## Acceptance checks for integration

1. Compare full primitive Z/P columns using the same retained meshes, inputs
   and executable. Cover low/middle/high frequency, detailed and homogenized
   cables, both material profiles, and a smaller multi-terminal fixture.
   Use a tight norm-based roundoff tolerance plus an absolute tolerance near
   zero, calibrated to the unchanged baseline. The existing 10% legacy
   cross-mesh reference threshold is much too loose for this change.
2. Check one LU setup/factorization per frequency and N solves with PETSc event
   counters. Record phase timings separately from optional field-map I/O.
3. Check source resetting, nonconsecutive/missing bases, reordered requested
   bases, failure after a complete column, interrupted writes, duplicate rows,
   wrong indices/Hz, and reuse of an already completed scan. Never infer
   success from process return or row count alone.
4. Compare reduced Z/P and final Y through the existing Engine path, including
   its condition-aware inversion residual. Check existing reciprocity and
   passivity diagnostics without modifying entries to force them to pass.
5. Exercise maps on/off, all expected names and labels, headless/UI session
   restoration, transport errors, and executable capability rejection.

The probes in this note have completed the retained-input numerical and PETSc
counter checks for eight frequency/profile combinations. They have not
integrated or tested the production resume/UI/map behavior, nor rerun a full
101-frequency campaign or the package test suite.

## Reproducing the diagnostic

The original scratch directory is `/tmp/lcm-fem-performance`, with per-probe
commands, logs, raw files and summaries. The durable patch can be applied to
a fresh copy of a retained run's `input/getdp` directory using `patch -p1` from
the directory containing that copy. It adds a guarded `ModelDataPath` input
for standalone execution; supply the path explicitly with `-setstring`.

For `run-Pbj2uz` frequency index 51, copy `input/model_data.pro` and
`mesh/frequency_0051.msh` too, and create fresh `raw/jobs` and `maps`
directories. Invoke the recorded executable with these arguments, replacing
`SCRATCH` with the absolute fresh directory:

```text
SCRATCH/getdp/model.pro
-solve LineCableModelsFEMScan
-msh SCRATCH/mesh.msh
-name SCRATCH/solver
-v 5 -cpu -ksp_view
-log_view :SCRATCH/petsc.log
-setstring ModelDataPath SCRATCH/model_data.pro
-setstring RunDirectory SCRATCH
-setstring RawOutputStem probe
-setnumber FrequencyIndex 51
-setnumber FrequencyHz 316.22776601683796
-setnumber Val_Rint 283.02195830623396
-setnumber Val_Rext 353.77744788279244
-setnumber BasisTerminal 1
-setnumber PlotFieldMaps 0
```

Capture stdout and stderr to the scratch log. Compare all 81 nonblank rows in
each `probe-Z.tsv` and `probe-P.tsv` with the nine retained
`getdp-f0051-bNNNN-{Z,P}.tsv` files by `(frequency, response, basis)`; validate
Hz, uniqueness, completeness and finiteness as well as numerical differences.
Frequency index 101 uses `mesh/model.msh` and `model.json` in these runs, not
`frequency_0101.msh`. Use each selected mesh's metadata for its exact Hz and
radii. Always use fresh outputs: the probe deliberately appends across calls.

## Reference review

All seven supplied documents under `LineCableLab/docs` were consumed:
`onelab-parsing.md`, `onelab-nonlinearsolvers.md`, `onelab-linearsolvers.md`,
`onelab-integration-postprocessing.md`, `onelab-syntax.txt`,
`onelab-json.txt`, and the complete `getdp-manual.txt` (GetDP 3.5 manual).
The parser/runtime split, assigned-constraint updates, file append semantics,
ONELAB value precedence, and system mesh declarations directly inform the
design above. The official online
[GetDP resolution reference](https://getdp.info/doc/texinfo/getdp.html#Types-for-Resolution)
also documents `SolveAgain`, `UpdateConstraint`, and separate-generation
operations; executable-specific behavior was checked experimentally.

Relevant implementation sources reviewed:

| Source beneath `/home/amartins/Documents/KUL` | Applicable finding |
|---|---|
| `linecablemodels-fem/linecablemodels.py`, `linecablemodels.geo`, `lib/electrodynamic_full_wave.pro`, magnetic/electric variants | Julia's current formulation and conformal mesh design have recognizable predecessors; the legacy solve macro also regenerates and solves for each basis |
| `LineCableLab/driver/linecablelab.py`, `gmsh/model.geo`, `getdp/{model,fullwave,magnetic,electric,jacobian_integration}.pro` | Standalone and ONELAB launch paths, stable client/launcher handling, isolated meshes and outputs; the driven scan still launches per frequency and terminal |
| `MultiLayerZY/fem/templates/cable.jl`, `cable.geo` | Direct Gmsh Julia API geometry/physical groups and characteristic lengths; alternative conformal fragmentation examples |
| `MultiLayerZY/fem/templates/thin-wire-mutual-eig/electrodynamic_hybrid.pro` | Runtime currents and constraint updates; gamma iteration is a distinct nonlinear mode that cannot inherit fixed-operator factor reuse blindly |
| `MultiLayerZY/fem/templates/models-getdp/ElectromagneticScattering/scattererTmatrix.pro` | First assembly/solve followed by `GenerateRHS` and `SolveAgain` for changing excitations |
| `MultiLayerZY/fem/templates/models-getdp/GetDDM/{Schwarz,SchwarzMacros}.pro` | Explicit per-system mesh ownership, restricted RHS assembly, reusable factors |
| `MultiLayerZY/fem/templates/models-getdp/ThermalConduction/Thermal.pro` | `GenerateSeparate`/`Update` reuse for its constant-coefficient transient problem; different assumptions from this scan |

The template trees were scanned for relevant solve, mesh, RHS and frequency
patterns. Gmsh Julia tutorial examples in the installed ONELAB distribution
also corroborate the API field construction. These examples supply particular
mechanisms; their alternative formulations and boundary conditions are not
substitutes for the working coupled model.
