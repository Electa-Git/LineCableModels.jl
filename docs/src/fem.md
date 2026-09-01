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
        getdp_executable = "/path/to/getdp",
        gmsh_verbosity = 2,
        getdp_verbosity = 2,
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

GetDP is an external executable. Set `getdp_executable` explicitly or make a
`getdp` executable available on `PATH`. The backend verifies its identity
before starting a solve. No Python runtime or `GetDP.jl` problem generator is
used.

## Execution model

One call to `compute` builds one complete two-dimensional Gmsh mesh for each
distinct physical frequency. Each mesh uses that frequency's earth skin depth
for its finite air/earth radius and has its own conformal annular
transformation-to-infinity shell. The final, highest-frequency mesh is the
displayed mesh; every frequency-specific mesh is reused by all of that
frequency's basis-terminal jobs. Julia launches one isolated GetDP process for
each ordered frequency/basis-terminal pair. Every process appends exactly one
response column to `raw/Z.tsv` and `raw/P.tsv`. Julia publishes the completion
marker only after every job has returned a complete column, then validates the
marker, headers, row count, indices, frequencies, and finite values.

Gmsh otherwise gives external ONELAB clients a shared user socket such as
`.gmshsock2`. The extension assigns every FEM run a unique socket basename in
the user directory (or an ephemeral loopback port on Windows) through
`Solver.SocketName`, then restores the caller's option. Independent Julia
processes can therefore mesh and solve concurrently without racing on the
default socket.

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

Disks, ellipses, circular sectors, and rounded sectors retain exact Gmsh
circle/ellipse arcs; rectangles and schema polygons retain exact line
segments. Annuli, conformal rounded-sector shells, enclosure differences, and
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
tables, maps, logger output, and atomic `run.json` state transitions.

Field maps are off by default. With `plot_field_maps=true`, nine supplied
quantities are written for every frequency/source pair, with names such as
`bm_f0002_b0003.pos`. Every expected file must exist before the scan succeeds.
Headless execution does not merge them; UI execution merges them only after the
complete numerical scan validates. Map paths are retained in result details
only when the run directory is retained.

The executable resolution order is an explicit `getdp_executable`, followed
by the executable named `getdp` on `PATH`; this repository has no additional
sanctioned GetDP preference. A nonzero client failure is reported as a typed
error with its frequency and basis indices, the retained run directory, and
the captured Gmsh log tail. If Gmsh reports an error only after a complete
response column was written, the error is deferred to strict table validation.
A normal client return is still not success until the completion marker and
every raw row have passed strict validation. Result details report one GetDP
invocation per frequency/source pair.

## Optional Gmsh UI

The UI is a visualization and debugging surface, not an input editor:

```julia
interactive_fem = Formulation(
    :LineCableModelsFEM;
    fem_options = (
        ui = true,
        plot_field_maps = true,
        getdp_executable = "/path/to/getdp",
    ),
)

parameters = compute(problem, interactive_fem)
```

It publishes read-only problem summaries, separate mesh and solve states, and
status text plus `Generate mesh` and `Run model` buttons. `Run model` refuses
to proceed before a valid mesh exists.
Closing the window before solving raises a typed `:not_executed` error that
distinguishes closure before mesh generation from closure after meshing. After
a successful scan, validated maps remain visible until the user closes the UI.

The extension finalizes only Gmsh sessions it owns. A caller-owned initialized
session retains its current model, unrelated models and views, Gmsh verbosity
options, and pre-existing `LineCableModels/FEM/` ONELAB parameters.

## Numerical reference validation

The committed `fem_python_quasi_tem.json` fixture freezes development-only
outputs from the supplied Python quasi-TEM prototype; Python is not imported
or executed by the package or its tests. The two cases use the same copper,
dielectric, earth, geometry, frequency ordering, and reductions as their
Julia runs:

| Case | Frequencies [Hz] | Primitive ``Z`` | Primitive ``P`` | Reduced ``Z`` | Final ``Y`` |
|---|---:|---:|---:|---:|---:|
| One coaxial cable, sheath Kron-reduced | 10, 1,000, 100,000 | 1.2340% | 0.1366% | 0.5746% | 0.2228% |
| Two coaxial cables, bundled cores and Kron-reduced sheaths | 10, 1,000, 100,000 | 1.2017% | 0.1322% | 0.5589% | 0.2155% |

Entries are relative Frobenius norms over the complete frequency scan. The
test limit is 10% at both primitive and reduced levels; the recorded values
are the unmodified comparison results from Gmsh 4.15 and GetDP 3.5.
