# Gmsh/GetDP finite-element backend

[`LineCableModelsFEM`](@ref) is the Julia-native finite-element
backend for `LineParametersProblem`. Gmsh is a weak dependency: the public
formulation and option normalizers are always available, while the `compute` method
is activated by loading Gmsh.

```julia
using LineCableModels
using Gmsh

fem = Formulation(
    :LineCableModelsFEM;
    insulation_admittance = formula(:default),
    semicon_admittance = formula(:default),
    earth_properties = formula(:default),
    temperature_dependence = formula(:default),
    options = (
        physics = :quasi_tem,
        reduce_bundle = true,
        kron_reduction = true,
        ideal_transposition = false,
    ),
)

parameters = compute(problem, fem;
    options = (
        mesh_policy = :reuse,
        gmsh_verbosity = 2,
        getdp_verbosity = 2,
        frequency_workers = 2,
        solver_threads = 1,
    ),
)
```

All execution controls pass through `compute(...; options=(...))` and are
validated by `computation_options(LineCableModelsFEM, ...)`. This includes
meshing, workers, native verbosity, Julia milestones, logs, and resume:

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
`plot_field_maps=true`, basis-specific output operations write the nine common field
quantities with frequency/source-specific filenames and labels. The maintained GetDP
files are captured with each run so later edits cannot change an active scan.

The former `fem_options` keyword and `LineCableModelsFEMOptions` struct have been
removed. Move their execution keys into `compute` options and `physics` into
formulation options. Benchmark calls use `reference_options` for FEM execution
controls. Supplemental run metadata is a named tuple under `details(result).fem.run`.

## Physics selection

Set formulation `options=(physics=:quasi_tem,)` for the default independent magnetic and
electric blocks, or `(physics=:quasi_fw,)` for the coupled first-order Maxwell
model. Strings `"quasi-tem"` and `"quasi-fw"` and their `Symbol` values are also
accepted. Julia requires `Symbol("quasi-fw")` for a hyphenated symbol;
`:quasi-fw` is parsed as subtraction.

```julia
fem = Formulation(:LineCableModelsFEM; options=(physics=:quasi_fw,))
parameters = compute(problem, fem)
```

`model.pro` exposes the ONELAB number `LineCableModels/FEM/physics`, with
choices `0 = quasi-tem` and `1 = quasi-fw`. It includes `quasi-tem.pro` or
`quasi-full.pro` accordingly. Headless jobs set the same constant with
`-setnumber Physics 0` or `-setnumber Physics 1`; both use the resolution
`LineCableModelsFEMScan`. The Julia-managed GUI displays this choice read-only,
since the numerical inputs of a run are fixed before meshing. Choose physics
in formulation `options` for a new calculation. Physics is saved with the formulation
inputs and column checkpoints; a run cannot resume under different physics.

The quasi-full option retains ``A_t/\Gamma`` and ``\phi/\Gamma`` in a first-order
``\Gamma\to0`` reduction of the vector-potential and continuity equations.
One axial-current excitation supplies both responses. Its voltage extraction
includes ``j\omega\int(A_t/\Gamma)\cdot d\ell`` along a physical vertical path
from earth infinity to each terminal. The backend constructs these paths from
the first-order triangular mesh, choosing the lowest node of each terminal
contour (lowest x breaks ties). Transverse fields are zero inside equipotential
metal, so paths can cross shields or other terminals. The infinite-shell part
uses the pulled-back edge field and 128 straight segments; integration within
each crossed triangle is exact for its lowest-order edge element. This is a
voltage-path convention, not a claim of path independence in an inductive field.

The full equations, units, boundary conditions, gauge and extraction are
documented in [`LineCableModelsFEM`](@ref), with the potential-equation reference
of [Ciuprina2024](@cite). This backend's 2D longitudinal reduction is distinct
from that paper's 3D ECE implementation. It retains displacement and does not
solve a finite-``\Gamma`` eigenproblem. Both public physics options retain the
specified finite metal conductivity in the axial problem.

Quasi-full additionally retains the scalar-only `Pscalar.tsv` per column,
which is gauge dependent and must not be inverted for Y. With field maps
enabled it also writes `bt_mesh`, `v_local` and `hz_scaled`, for twelve maps
per excitation. Its `e`, `em` and `jm` maps represent ``E_t/\Gamma`` [V], its
magnitude [V], and ``|J_t/\Gamma|`` [A/m], respectively.

## Material laws

FEM selects insulation and semicon admittivity, soil frequency dependence, and
cable-material temperature dependence. Analytical `internal_impedance`,
`insulation_impedance`, `earth_impedance`, `earth_admittance`, and
`pipe_impedance` keywords are rejected, including explicit `:default` values.
Supported enclosing geometry is represented directly in the field domain.

`temperature_dependence=formula(:default)` evaluates each cable material's
resistivity as ``\rho(T)=\rho_0[1+\alpha(T-T_0)]``. `T` comes from
`problem.temperature`; reference resistivity, `T0`, and `alpha` come from the
material. Select `nothing` to retain reference resistivity. This law is shared
with analytical calculations, cable constants, and PSCAD export. The retired
`options.temperature_correction` Boolean is rejected: replace `true` with the
`:default` temperature selection and `false` with `nothing`.

The default law requires a positive finite correction factor and
``|T-T_0|<150`` K. These are limits of this approximation, independent of thermal
rating. Custom laws own their applicability and use the usual contribution-hook
contract; no FEM author registration is required. Passive materials can retain
infinite resistivity. Dielectric constituents are evaluated before radial
aggregation, and polarization loss is not corrected a second time.

`earth_properties` calls the same soil constitutive law as the analytical engine
at each frequency. Evaluated resistivity, permittivity, and permeability feed
both GetDP's soil/infinite-soil coefficients and the skin-depth mesh rule.
`:default` and `nothing` retain static soil. Air uses its explicitly declared
static permittivity and permeability, independently of the soil law. The current
FEM geometry requires one horizontal semi-infinite soil; a non-finite conductive
skin depth is unsupported. Equivalent homogeneous-earth reductions are rejected.
An ordinary `EarthModel` supplied after an external reduction carries no history
from which FEM could detect that prior approximation.

Saved FEM formulation details contain only the four consumed `selections`, their
parameters, numerical options, and hook descriptions. Custom hooks are identified
but marked nonreplayable; saved records do not reconstruct executable closures.
The selected propagation approximation remains recorded separately.

## Field equations and matrix extraction

The default `:quasi_tem` physics evaluates the series and shunt problems at ``\Gamma=0``. It retains
diffusion and displacement in the surrounding media, with phasors proportional
to ``e^{j\omega t}`` and complex admittivity ``\kappa=\sigma+j\omega\epsilon``.
Two independent blocks share one assembled GetDP system and factorization.

The magnetic block solves for the axial vector potential ``A_z`` and one axial
electric unknown ``u_i`` per terminal. In each material it solves

```math
-\nabla\cdot(\mu^{-1}\nabla A_z)+\kappa(j\omega A_z+u_i)=0,
\qquad
I_i=-\int_{\Omega_i}\kappa(j\omega A_z+u_i)\,dS.
```

Here ``u_i`` is supported on its conductor region. Exciting terminal ``s`` with
1 A and imposing zero axial current on the others gives ``Z_{is}=-u_i/I_s``.
Metal conductivity and its internal field remain part of this series problem.

The electric block solves the scalar electrodynamic problem in air, soil and
passive cable materials, excluding conductor interiors:

```math
\nabla\cdot(\kappa\nabla v)+\kappa k^2 v=0,
\qquad k^2=-j\omega\mu\kappa.
```

Each terminal has one equipotential degree of freedom ``V_i``. Prescribing
1 A/m of outward transverse terminal current at terminal ``s``, and zero at the
others, gives the column ``P_{is}=V_i/(1\ \mathrm{A/m})``. GetDP's associated
quantity has the opposite sign, so this drive is imposed as ``Q_s=-1``.
Here ``Q`` is a current per unit length, not an electrostatic charge;
``P`` has units Ω m and ``Y=P^{-1}`` has units S/m.

This is equivalent to prescribing a unit voltage on each source terminal in
turn, grounding the others, and extracting ``Y_{is}=-Q_i/(1\ \mathrm{V})``.
The native regression checks both excitations independently. For any other
set of voltage excitations, terminal currents satisfy ``J=YV``; extracting
``Y`` requires the complete voltage matrix, not just a source-voltage rescaling.

The magnetic and electric blocks are separate unit excitations. Electric
potentials are obtained directly from the scalar operator, without dividing a
magnetically driven potential by a small ``\Gamma``. This also keeps the shunt
calculation independent of metal-interior discretization. Bare conductors have
no passive coating in their electric domain; finite-conductivity metal remains
in their magnetic domain.

The earlier coupled `A_z/u/phi` model retained only the axial vector potential
and used continuity to recover `Phi/Gamma`. At material interfaces that
reduction omits a transverse Ampère balance at the same order in Γ as the
electric response being extracted. With
``\mathbf r=\kappa\nabla_t(\phi/\Gamma)-\mu^{-1}\nabla_t A_z``, continuity
enforces ``\nabla_t\cdot\mathbf r=0``, whereas transverse Ampère requires
``\mathbf r=0``. A material interface makes those conditions inequivalent;
reducing Γ does not remove the error after normalization by Γ.

## Execution model

One call to `compute` builds one complete two-dimensional Gmsh mesh for each
distinct physical frequency. Each mesh uses that frequency's earth skin depth
for its finite air/earth radius and has its own conformal annular
transformation-to-infinity shell. The final, highest-frequency mesh is the
displayed mesh; every frequency-specific mesh is reused by all of that
frequency's terminal excitations. Cable topology is constructed once: successive
meshes retain its vertices and surfaces while updating the exterior circles and
mesh-size fields. Only one native geometry/mesh is retained in memory. After
preparing all meshes, Julia launches
up to `frequency_workers` standalone GetDP processes concurrently. Each process
handles one frequency: it assembles and factors the system for its first
requested terminal, then updates the right-hand side and reuses those factors
for the remaining terminals. Frequencies with different meshes have separate
systems and factorizations. Local mesh sizes also resolve the attenuation and
phase scales of the evaluated air and soil properties.

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
pattern counts remains integral. The caller-owned problem is not mutated. Evaluated material-law outputs pass the
same checked nominal `Float64` boundary before transport; finite overflow is
rejected. Analytical scalar and uncertainty propagation remain unchanged.

| FEM datum | Authoritative LineCableModels property | Handling |
|---|---|---|
| Material class and electrical properties | each resolved `PlacedRegion.source.material`: `kind`, `rho`, `eps_r`, `mu_r`, `tan_delta`, `T0`, `alpha` | Reused; resistivity follows the selected shared temperature law and constant intrinsic loss tangent contributes ``\omega\epsilon\tan\delta`` to conductivity |
| Material geometry and topology | `CableDesign.geometry.regions`, each resolved `PlacedRegion.primitive`, and `CableDesign.geometry.outer` | Requires an area-complete material partition, then adapts it to built-in `gmsh.model.geo` loops and cut-hole surfaces |
| Cable identity and placement | `LineCableSystem.designs`, `CableDesign.cable_id`, `LineCableSystem.positions`, and resolved `LineCableSystem.geometry` | Reused in declared order; stable IDs form physical names |
| Terminal ownership and order | `LineCableSystem.terminal_order`, `terminal_map`, and `connection_order` | Reused exactly; disconnected surfaces of one electrical Group share one terminal physical group |
| Phase, bundle, and grounded-conductor reduction | `LineCableSystem.connection_order` plus shared formulation options | Delegated to the Engine reduction implementation for both ``Z`` and ``P`` |
| Frequencies | `LineParametersProblem.frequencies` | Published as indexed, hidden, read-only ONELAB numbers; each isolated GetDP job also receives its exact physical frequency and matching transformation radii directly |
| Temperature | `LineParametersProblem.temperature` | Prescribed input to the selected temperature law; the default uses material `T0` and `alpha` |
| Earth material | `LineParametersProblem.earth_props` | Declared air plus one horizontal soil half-space; the soil law is evaluated per frequency |
| Optional environment declaration | `LineCableSystem.environment` | `nothing` and `EarthModel` are accepted; other declarations produce a typed unsupported-feature error |
| Line length and output basis | `LineCableSystem.line_length` and shared `compute` options | Per-unit-length is canonical; total basis uses the existing package scaling |
| Propagation constant | backend-owned ``\Gamma\to0`` limit, with independent or first-order coupled fields selected by formulation `options.physics` | A non-`nothing` problem-level `Γ` is rejected rather than silently reinterpreted |
| Mesh resolution | local characteristic lengths derived from each resolved solid, tube, strand, foil, and passive region; per-frequency earth skin depth controls the exterior domain, and air/soil propagation scales constrain surrounding-medium resolution | Thin internal features remain local and cannot refine unrelated layers or the earth domain |

Disks, ellipses, and cable sectors retain exact Gmsh circle/ellipse arcs;
rectangles and schema polygons retain exact line segments. Annuli, conformal
sector shells, enclosure differences, and
assembly boundaries use shared oriented loops. Circular boundaries are
pre-segmented at sector endpoints and circle contacts, so adjacent materials
reuse the same curve and tangent strands reuse the same point. Compacted strand
polygons are used unchanged. Touching hole boundaries are partitioned into
connected filler faces; metal-metal seams are excluded from filler boundaries.
All filler faces retain their declared material. Equal evaluated material laws
share a physical material group, independently of geometric strand identity and
electrical terminal groups. A shared
material interface takes the smaller of its two local characteristic lengths.
Thin internal foils and strands do not export their size to the cable/earth
boundary. One `Distance`/`Threshold` field per actual cable exterior grows
from that exterior layer's size to `domain_radius/20` using an adjacent-element
growth factor of 1.2. Additional fields restrict the surrounding-medium size to
``h\leq 1/(8|q|)``, where ``q=\sqrt{j\omega\mu\kappa}``, within six attenuation
lengths of cable exteriors and the air/soil interface. The bound resolves both
decay and phase; it transitions back to the existing domain size beyond that
distance. A lossless medium keeps its phase-resolution bound across the domain.
Gmsh `Restrict` fields apply each bound to its own air or soil surfaces, and
`Min` combines overlapping fields. No artificial refinement rings are
introduced. The adapter
rejects an incomplete area partition before starting Gmsh and rejects any
internal material curve without exactly two adjacent surfaces after
synchronization. After meshing, boundary-edge incidence and material coverage
are checked before a mesh can be cached or passed to GetDP. Nonempty but partial
material meshes are rejected.

Rectangular stranded cores supply their occupied disk boundary directly from
physical resolution. FEM uses the same boundary as preview, analytical
flattening and subsequent layers. Complete bounded formations are recognised
using the same floating-point area tolerance as enclosure resolution, not an
independent engineering fill-fraction cutoff; retained filler is not replaced
by an expanded conductor.

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
terminal physical groups, physical names, complete boundary incidence, material
areas (with curved-boundary discretization allowances), and conductor ownership.
Owned MSH 4.1 files retain all boundary elements, including same-material seams;
explicit mesh files must retain these elements too. `mesh_policy=:remesh` always
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

Maps `az`, `b`, `bm`, `ez`, `jz`, and `rhoj2` describe the axial 1 A drive.
Maps `e`, `em`, and `jm` describe the transverse 1 A/m drive: respectively
``-\nabla v``, its magnitude, and ``|\kappa\nabla v|`` in the surrounding media.
Their view labels identify the drive. These axial and transverse fields belong
to different excitations and do not form a single full-wave field vector.

The executable resolution order is:

1. `compute(...; options=(getdp_executable="/absolute/path/to/getdp",))`;
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
for calculation records, while
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
preserved comparison artifacts and require a fresh computation. Indexed soil and
declared-air coefficients use run-input schema 7 and solver protocol 3; older
schemas cannot resume. Evaluated cable, soil, and air coefficients participate
in solve reuse identity. Numerically identical laws can share a solve while
retaining separate selection calculation records and independent result arrays.

## Optional Gmsh UI

The UI is a visualization and debugging surface, not an input editor:

```julia
interactive_fem = Formulation(:LineCableModelsFEM)
parameters = compute(problem, interactive_fem;
    options = (
        ui = true,
        plot_field_maps = true,
    ),
)
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

The two-bare-wire Gauntlet benchmark runs fresh native FEM against the default
and Xue analytical formulations at eight frequencies from 0.1 Hz to 1 MHz.
Its geometry contains only two metal disks: radius 4.25 cm, depth 1 m, separation
1 m, in soil with resistivity 0.1 Ω m. There is no insulation or fitted FEM
table. Run it through `dev/run_two_bare_wires.jl` for the shared report and
PlotBuilder comparison plots.

Native regressions compare every complex self and mutual entry at 100 kHz and
1 MHz, where a matrix norm alone can hide a mutual-admittance error. They also
check reciprocity, independent unit-current and unit-voltage extraction,
invariance of shunt admittance to metal conductivity, factorization reuse,
and overhead and mixed conductor layouts. Existing coaxial-capacitance,
constitutive-law, enclosure, reduction and recovery tests cover the surrounding
backend behavior.

The committed `fem_python_quasi_tem.json` retains historical outputs from the
earlier coupled Python prototype for two insulated coaxial cases. Its optional
comparisons are legacy checks with a 10% matrix-norm tolerance, not the reference
for the bare-wire electric formulation. The backend has no Python dependency.
