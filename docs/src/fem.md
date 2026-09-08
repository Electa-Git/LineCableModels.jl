# Gmsh/GetDP finite-element backend

[`LineCableModelsFEM`](@ref) is the Julia-native, coupled quasi-TEM finite-element
backend for `LineParametersProblem`. Gmsh is a weak dependency: the public
formulation and option types are always available, while the `compute` method
is activated by loading Gmsh.

```julia
using LineCableModels
using Gmsh

fem = Formulation(
    :LineCableModelsFEM;
    options = (
        reduce_bundle = true,
        kron_reduction = true,
        ideal_transposition = false,
        temperature_correction = true,
    ),
    fem_options = (
        mesh_policy = :reuse,
        gmsh_verbosity = 2,
        getdp_verbosity = 2,
        frequency_workers = 2,
        solver_threads = 1,
    ),
)

parameters = compute(problem, fem)
```

Julia-side milestones and an optional Julia log file remain computation
options, independent of the two external verbosity controls:

```julia
parameters = compute(
    problem,
    fem;
    options = (
        verbosity = (default = 1,),
        log_file = "fem-julia.log",
    ),
)
```

GetDP remains an external process, but no separate installation is required on
the supported artifact platforms. The first FEM calculation there downloads
the package's lazy, hash-verified GetDP 3.5.0 complex-PETSc artifact. Loading
LineCableModels or Gmsh alone does not download it. No Python runtime or
`GetDP.jl` problem generator is used.

The adapter writes one run-local `input/model_data.pro` containing resolved
region tags, terminal names, material coefficients and domain dimensions.
Maintained GetDP files own material-domain bindings, the equations and the
parameterized field-map operation. Material conductivity is the effective
real part of the selected complex admittivity: dielectric losses are already
included, so GetDP does not add another loss-tangent contribution.

Each solver job selects its material coefficients by frequency index. When
`plot_field_maps=true`, basis-specific output operations write the nine field
quantities with frequency/source-specific filenames and labels. The maintained GetDP
files are captured with each run so later edits cannot change an active scan.

## Execution model

One call to `compute` builds one complete two-dimensional Gmsh mesh for each
distinct physical frequency. Each mesh uses that frequency's earth skin depth
for its finite air/earth radius and has its own conformal annular
transformation-to-infinity shell. The final, highest-frequency mesh is the
displayed mesh; every frequency-specific mesh is reused by all of that
frequency's terminal excitations. After preparing all meshes, Julia launches
up to `frequency_workers` standalone GetDP processes concurrently. Each process
handles one frequency: it assembles and factors the system for its first
requested terminal, then updates the right-hand side and reuses those factors
for the remaining terminals. Frequencies with different meshes have separate
systems and factorizations. The equations and mesh sizing are unchanged.

The default is two frequency workers with one BLAS/OpenMP thread per solver.
These are independent OS processes; they do not require multiple Julia threads.
Set `frequency_workers=1` for serial frequency execution. Increase the worker
count only within available memory: every active frequency owns its sparse
factorization. `solver_threads` sets each child's BLAS/OpenMP environment without
changing Julia's environment. GetDP must support `GenerateRHSGroup` and
`SolveAgain`; the package-owned artifact is GetDP 3.5.0 with PETSc complex
arithmetic.

Every attempt has a separate working directory, solver prefix, raw columns,
maps and log. Inputs are explicit command arguments and data files; the solvers
do not connect to ONELAB. One Julia coordinator owns Gmsh, publishes UI progress,
validates columns and writes checksummed checkpoints. Backend calls in the same
Julia process serialize access to the shared Gmsh session. A filesystem lock
prevents two coordinators from writing the same run.

Once every column is validated, Julia assembles `raw/Z.tsv` and `raw/P.tsv` in
frequency/terminal order and publishes the scan completion marker. Validation
checks headers, row counts, identities, frequencies and finite values.

The primitive matrices have dimensions
`(nterminal, nterminal, nfrequency)`. The shared LineCableModels reduction path
applies terminal ordering, bundle merging, Kron reduction, and ideal
transposition to both ``Z`` and the potential-coefficient matrix ``P``. The
backend then obtains ``Y`` by a condition-checked direct solve of
``P Y = I``—there is no additional ``j\omega`` factor.

The returned value is the package-native `LineParameters` in `PhaseDomain`,
with ``Z`` in Ω/m, ``Y`` in S/m, and the exact input frequency vector in Hz.
Pass `options=(trace=true,)` to `compute` to retain primitive ``Z/P`` in result
details. `output_basis=:total` continues to use the shared computation option
and scales by line length.

## Schema authority and reconciliation

The existing LineCableModels typed objects are the sole physical-model
authority. The extension creates only derived tags, surface ownership, mesh
sizes, and GetDP tables; it defines no second cable, material, or project
schema.

Every FEM computation starts with a numeric preflight before model adaptation,
runtime-directory creation, Gmsh initialisation, or meshing. The preflight
rebuilds continuous problem data as `Float64`; when a
`Measurements.Measurement` scalar is present, only its nominal value is
retained. Discrete topology such as terminal assignments, material tags, and
pattern counts remains integral. The caller-owned problem is not mutated.

| FEM datum | Authoritative LineCableModels property | Handling |
|---|---|---|
| Material class and electrical properties | each resolved `PlacedRegion.source.material`: `kind`, `rho`, `eps_r`, `mu_r`, `tan_delta`, `T0`, `alpha` | Reused; resistivity follows the shared temperature-correction option and constant intrinsic loss tangent contributes ``\omega\epsilon\tan\delta`` to conductivity |
| Material geometry and topology | `CableDesign.geometry.regions`, each resolved `PlacedRegion.primitive`, and `CableDesign.geometry.outer` | Requires an area-complete material partition, then adapts it to built-in `gmsh.model.geo` loops and cut-hole surfaces |
| Cable identity and placement | `LineCableSystem.designs`, `CableDesign.cable_id`, `LineCableSystem.positions`, and resolved `LineCableSystem.geometry` | Reused in declared order; stable IDs form physical names |
| Terminal ownership and order | `LineCableSystem.terminal_order`, `terminal_map`, and `connection_order` | Reused exactly; disconnected surfaces of one electrical Group share one terminal physical group |
| Phase, bundle, and grounded-conductor reduction | `LineCableSystem.connection_order` plus shared formulation options | Delegated to the Engine reduction implementation for both ``Z`` and ``P`` |
| Frequencies | `LineParametersProblem.frequencies` | Published as indexed, hidden, read-only ONELAB numbers; each isolated GetDP job also receives its exact physical frequency and matching transformation radii directly |
| Temperature | `LineParametersProblem.temperature` | Reused through material `T0` and `alpha` when temperature correction is enabled |
| Earth material | `LineParametersProblem.earth_props` | One homogeneous horizontal earth half-space is adapted to the FEM domain |
| Optional environment declaration | `LineCableSystem.environment` | `nothing` and `EarthModel` are accepted; other declarations produce a typed unsupported-feature error |
| Line length and output basis | `LineCableSystem.line_length` and shared `compute` options | Per-unit-length is canonical; total basis uses the existing package scaling |
| Propagation constant | backend-owned fixed quasi-TEM constant | A non-`nothing` problem-level `Γ` is rejected rather than silently reinterpreted |
| Mesh resolution | local characteristic lengths derived from each resolved solid, tube, strand, foil, and passive region; per-frequency earth skin depth controls only the exterior domain | Thin internal features remain local and cannot refine unrelated layers or the earth domain |

Disks, ellipses, and cable sectors retain exact Gmsh circle/ellipse arcs;
rectangles and schema polygons retain exact line segments. Annuli, conformal
sector shells, enclosure differences, and
assembly boundaries use shared oriented loops. Circular boundaries are
pre-segmented at sector endpoints and circle contacts, so adjacent materials
reuse the same curve and tangent strands reuse the same point. A shared
material interface takes the smaller of its two local characteristic lengths.
Thin internal foils and strands do not export their size to the cable/earth
boundary. One `Distance`/`Threshold` field per actual cable exterior grows
from that exterior layer's size to `domain_radius/20` using an adjacent-element
growth factor of 1.2; overlapping fields are combined with `Min`. No artificial
refinement rings are introduced. The adapter
rejects an incomplete area partition before starting Gmsh and rejects any
material curve lacking a neighbouring field surface after synchronization,
before meshing or invoking GetDP.

The current FEM domain explicitly rejects vertical earth layers, more than one
earth half-space layer, a problem-supplied propagation constant, unsupported
environment types, incomplete material partitions, and any resolved primitive
without a two-dimensional built-in-`geo` boundary adaptation. These failures use
`LineCableModelsFEMError`, including the owning object ID and offending field,
before Gmsh is touched where possible.

## Mesh lifecycle and diagnostics

`mesh_policy=:reuse` first validates an explicit highest-frequency `mesh_path`,
then checks the fingerprinted repository-local cache for each frequency, and
otherwise generates the missing frequency-specific mesh.
Compatibility checks cover mesh dimension, terminal count, material and
terminal physical groups, and physical names. `mesh_policy=:remesh` always
regenerates and atomically refreshes the matching cache. The fingerprint
includes the serialized problem, stable physical metadata, every local and
exterior mesh size, the physical mesh frequency, transformation radii, growth
law, and Gmsh version.

Runs live under `.linecablemodels/fem/runs/`; cached meshes live under
`.linecablemodels/fem/meshes/`. A successful run directory is deleted after the
result is constructed unless `keep_run_directory=true`. Failed or incomplete
runs are retained, and their typed error reports the path. Retained runs contain
the problem snapshot, immutable GetDP data, mesh snapshot and metadata, raw
tables, maps, logger output, and atomic `run.json` state transitions. Numerical
process logs and attempt metadata live in `attempts/fNNNN-*/`; per-column timing
records separate constraint updates, assembly, solve and output. At
`getdp_verbosity>=4`, each attempt also retains PETSc profiling output.

Field maps are off by default. With `plot_field_maps=true`, nine supplied
quantities are written for every frequency/source pair, with names such as
`bm_f0002_b0003.pos`. Every expected file must exist before the scan succeeds.
Headless execution does not merge them; UI execution merges them only after the
complete numerical scan validates. Map paths are retained in result details
only when the run directory is retained.

The executable resolution order is:

1. `fem_options=(getdp_executable="/absolute/path/to/getdp",)`;
2. the `LINECABLEMODELS_GETDP` environment variable;
3. the package's GetDP 3.5.0 lazy artifact; and
4. `getdp` on `PATH` only when the current platform has no artifact binding.

The artifact currently supports glibc Linux and macOS on x86-64. GetDP 3.5.0
publishes its Windows build only as a ZIP, which Julia's artifact installer
cannot consume directly; Windows therefore uses an installed GetDP selected
explicitly, through the environment variable, or on `PATH`. The same external
selection applies on every other unsupported platform. An explicitly selected
or environment-selected invalid path is an error; it is never silently
replaced by another solver. The backend records the resolved source and path
for provenance, while
resume compatibility uses the executable SHA-256 and reported build identity
instead of its filesystem location. See
[`THIRD_PARTY_NOTICES.md`](https://github.com/Electa-Git/LineCableModels.jl/blob/main/THIRD_PARTY_NOTICES.md)
for GetDP's GPL notice and upstream source location.

A nonzero client failure is reported as a typed
error with its frequency, missing basis indices, retained attempt directory
and GetDP log tail. A failure stops scheduling and terminates/reaps the other
active workers. Completed columns remain available for recovery. A zero exit
code is insufficient without valid completion records and numerical output.
Result details distinguish actual process launches (`getdp_invocations`),
`completed_columns`, and `completed_frequencies`. A fresh complete scan normally
launches one process per frequency; retries add invocations.

Resume an interrupted compatible run with
`options=(resume_run_directory="/path/to/run",)` (or `:latest`). Recovery checks
mesh identities and column checksums, adopts complete attempt outputs, and
requests only missing or invalid terminal columns. The first requested column
always builds fresh factors, even when its terminal index is not one. Worker
count may change during recovery; solver thread settings, physical inputs and
source/executable identities must match. A surviving solver from an interrupted
coordinator prevents retry until it exits. Completed runs are reused read-only
after their aggregate checksums pass. Runs from older solver protocols remain
preserved comparison artifacts and require a fresh computation.

## Optional Gmsh UI

The UI is a visualization and debugging surface, not an input editor:

```julia
interactive_fem = Formulation(
    :LineCableModelsFEM;
    fem_options = (
        ui = true,
        plot_field_maps = true,
    ),
)

parameters = compute(problem, interactive_fem)
```

It publishes read-only problem summaries, separate mesh and solve states, and
status text, completed frequency/column counts, and `Generate mesh` and
`Run model` buttons. `Run model` refuses
to proceed before a valid mesh exists.
Closing the window before solving raises a typed `:not_executed` error that
distinguishes closure before mesh generation from closure after meshing.
The event loop remains active while solver processes run. Closing the window
during solving cancels those processes and retains completed checkpoints. After
a successful scan, validated maps remain visible until the user closes the UI.

The extension finalizes only Gmsh sessions it owns. A caller-owned initialized
session retains its current model, unrelated models and views, Gmsh verbosity
options, and pre-existing `LineCableModels/FEM/` ONELAB parameters.

## Numerical reference validation

The committed `fem_python_quasi_tem.json` fixture freezes development-only
outputs from the supplied Python quasi-TEM prototype. These numerical comparisons
do not execute the prototype, and the backend has no Python dependency.
The two cases use the same copper,
dielectric, earth, geometry, frequency ordering, and reductions as their
Julia runs:

| Case | Frequencies [Hz] | Primitive ``Z`` | Primitive ``P`` | Reduced ``Z`` | Final ``Y`` |
|---|---:|---:|---:|---:|---:|
| One coaxial cable, sheath Kron-reduced | 10, 1,000, 100,000 | 1.2340% | 0.1366% | 0.5746% | 0.2228% |
| Two coaxial cables, bundled cores and Kron-reduced sheaths | 10, 1,000, 100,000 | 1.2017% | 0.1322% | 0.5589% | 0.2155% |

Entries are relative Frobenius norms over the complete frequency scan. The
test limit is 10% at both primitive and reduced levels; the recorded values
are the unmodified comparison results from Gmsh 4.15 and GetDP 3.5.
