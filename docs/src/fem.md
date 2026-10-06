# Gmsh/GetDP finite-element backend

[`FEM.LineCableModelsFEM`](@ref) is the finite-element backend for
`LineParametersProblem`. Loading the weak dependency Gmsh activates its
formulation, computation options, solver, native file readers and detached
export. `LineCableModelsGmshExt` imports and extends the main package's interfaces
to provide these definitions, without a dependency from the main package
on the extension.

Use `Formulation(:LineCableModelsFEM)` after `using Gmsh`. The concrete
formulation, error, mesh and field types are no longer main-package exports.
For explicit type access, use
`FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)` and
`FEM.LineCableModelsFEM`, `FEM.LineCableModelsFEMError`, `FEM.FEMMesh` or
`FEM.FEMFieldMap`. Ordinary callers use `compute`, `import_data`, `export_data`
and `plot` without direct access to extension types. Spatial plotting activates when
Gmsh and a Makie backend are loaded, in either order.

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
        physics = :helmholtz,
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
validated by `computation_options(FEM.LineCableModelsFEM, ...)`. This includes
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
region tags, terminal names, physical geometry, evaluated material coefficients
and normalized controls. Shared native `parameters.pro` computes exterior
dimensions, transverse roots, PML strengths and physical mesh targets;
`geometry.geo` and `mesh.geo` construct and mesh that case.
Maintained GetDP files own material-domain bindings, the equations and the
parameterized field-map operation. Material conductivity is the effective
real part of the selected complex admittivity: dielectric losses are already
included, so GetDP does not add another loss-tangent contribution.

Each solver job selects its material coefficients by frequency index. When
`plot_field_maps=true`, basis-specific output operations write scalar, vector and material-region field
quantities with frequency/source-specific filenames and labels. The maintained GetDP
files are captured with each run so later edits cannot change an active scan.

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
signature; no FEM author registration is required. Passive materials can retain
infinite resistivity. Dielectric constituents are evaluated before radial
aggregation, and polarization loss is not corrected a second time.

`earth_properties` calls the same soil constitutive law as the analytical engine
at each frequency. Evaluated resistivity, permittivity, and permeability feed
both GetDP's physical-soil/PML coefficients and the base-root sizing rule.
`:default` and `nothing` retain static soil. Air uses its explicitly declared
static permittivity and permeability, independently of the soil law. The current
FEM geometry requires one horizontal semi-infinite soil; lossless earth uses its wavelength and the native sizing ceiling. Equivalent homogeneous-earth reductions are rejected.
An ordinary `EarthModel` supplied after an external reduction carries no history
from which FEM could detect that prior approximation.

Saved FEM formulation details contain the four consumed `selections`, their
parameters and numerical options, alongside requested and resolved records.
Custom selections retain their identities and data; saved records never
reconstruct executable methods.
The selected propagation approximation remains recorded separately.

## Field equations and voltage measurement

The supported model is Helmholtz with phasors ``\exp(j\omega t-\Gamma z)``.
The prescribed longitudinal propagation constant ``\Gamma`` has units m⁻¹;
zero is the default. The normalized potentials ``A_z=a``, ``A_t=\Gamma b``
and ``\phi=\Gamma v`` retain every ``\Gamma^2`` term and have a regular
zero-Γ limit without numerical division by Γ. Finite metal retains the axial
field and total-current constraint; its metal contour is equipotential.
The complete equations and terminal conditions are specified in
[`FEM.LineCableModelsFEM`](@ref), with the potential-equation reference
[Ciuprina2024](@cite).

One imposed axial current excites each column. The native measurement is

```math
P_{ij}=\frac{1}{I_j}\int_{\ell_i}(\nabla v+j\omega b)\cdot d\ell,
\qquad Z_{ij}=K_{ij}+\Gamma^2P_{ij},\qquad Y=P^{-1},
```

where ``K_{ij}=-U_i/I_j`` is the longitudinal drive coefficient.
``P`` is in Ω·m, ``Z`` in Ω/m and ``Y`` in S/m. GetDP assembles the complex
primitive matrices before applying connection ordering, bundle reduction,
Kron reduction and ideal transposition. The scalar-only `Pscalar.tsv` is a
separate gauge-dependent diagnostic and cannot be inverted for Y.

Coordinates are x horizontal, y vertical and z axial, with the interface at
y=0. Each independent floating line runs upwards from the local interface for
an overhead receiver, or from the outer bottom boundary for a buried receiver.
It ends inside its own terminal metal, including the appropriate shell or
sector in a cable with multiple terminals. The horizontal displacement is
`min(1e-5,0.01*d)` metres, where `d` is the terminal's smallest metal dimension.
Its points and nodes are independent of the two-dimensional mesh and gauge tree.
No receiver columns divide the bottom PML.

For each excitation GetDP stores the full complex vector
`grad v + j*omega*bt` over media and PML in memory, then integrates its vertical
component with four-point line quadrature. There are no endpoint scalar terms
or extra stretch factors. Metal interiors have no stored field and contribute
zero. Physical line sizing uses `MeasurementLineSizeFactor=0.25` times the
minimum exterior target from every intersected cable neighbourhood, exterior
contour, wave and decay fields and earth interface layer. Conductor-interior skin
sizes do not constrain the measurement line. Buried PML segments use
`ceil(N_bottom/MeasurementLineSizeFactor)` intervals with the bottom progression.
Managed mesh metadata records the achieved outside-metal line/2D edge ratios;
a prescribed size factor alone does not prove measurement convergence.

## Physical earth sizing length

The earth base root is
``q_{0,e}=\sqrt{j\omega\mu_e(\sigma_e+j\omega\epsilon_e)}``, with
nonnegative real part. Its physical sizing length is

```math
L_0=\min\left(\frac{1}{\Re q_{0,e}},\frac{2\pi}{|q_{0,e}|}\right),
\qquad L_{\rm cap}=\sqrt{\frac{2\rho_{\rm cap}}{\omega\mu_e}},
\qquad \rho_{\rm cap}=10^5\ {\rm \Omega\,m}.
```

A zero real part means infinite decay length without evaluating division by
zero. The wavelength then supplies the physical length. The reference-earth
ceiling limits the near-lossless quasi-static extreme; it does not alter
material values, Γ, field equations or roots. Native sizing uses

```math
L=\min(L_0,L_{\rm cap}),\qquad
D=\max\left(R_{\rm layout},\ mL\frac{|q_{0,e}|}
{\max(|q_{0,e}|,|q_e|)}\right).
```

`domain_size_factor` is the positive dimensionless multiplier ``m``; its default
is 2. A larger prescribed transverse root can shorten the domain, while a
small or cutoff root cannot enlarge it. The layout bound includes the existing
5 m minimum and cable clearances. Each PML thickness is D times its directional
the native directional thickness factors. When the ceiling binds, the native flag is published
and managed execution warns
“earth too resistive for FEM domain sizing; results not qualified”.

## Staged coordinate stretching

Each medium has
``q_m^2=j\omega\mu_m(\sigma_m+j\omega\epsilon_m)-\Gamma^2`` and selected
root ``q_m=a_m+jb_m``. The sheet rule is per medium: choose the root with
``\Im q_m\geq0`` when ``\Im\Gamma<b_{0,m}``, otherwise the root with
``\Re q_m\geq0``. Here ``b_{0,m}=\Im q_{0,m}``. Squared-root real parts alone
do not distinguish radiation from a slow lossy wave. Air conductivity is a
native input, and the two permeabilities are independent.

For outward normalized layer coordinate u, the staged profile is

```math
s_d(u)=1+A_d\left(u^3-j\frac94\eta_d u^8\right).
```

Real stretching starts cubically and wave absorption is delayed to degree
eight. A medium's need for imaginary stretching is
``\nu_m=\operatorname{clamp}((b_m-a_m)/b_m,0,1)`` for ``b_m>0`` and zero
otherwise. Side layers use air and earth, top uses air and bottom uses earth.
Their ``\eta_d`` is the maximum participating need, capped at 1 and, for every
participating ``b_m<0``, at ``2a_m/(9|b_m|)``. Media with positive ``\nu_m`` require imaginary stretching, and the shared
side coefficient matches the interface.

Let ``C_m=(a_{0,m}+b_{0,m})/b_{0,m}``,
``d_{m,d}=(a_m+\eta_d b_m)/C_m`` and
``\widehat d_{m,d}=\max(d_{m,d},0.1b_{0,m})``. Strengths are calibrated by

```math
T=-\log(\mathrm{pml\_reflection})/2,\qquad
A_d=\frac{4T}{L_d\min_{m\in M_d}\widehat d_{m,d}}.
```

The default reflection target is ``10^{-3}``. The `0.1*b0` sizing-rate floor
caps the numerical strength to prevent factorization failure. The physical
root and terminal-conductance accuracy criterion remain unchanged. Floor activation is informational.
The flag uses ``d<(1-10^{-9})0.1b_0`` to ignore last-bit rounding.

Near DC, the added real stretch is limited to 100 domain half-widths.
Where the air field is quasi-static over that whole extent, the stretch is
purely real and ends at the homogeneous Dirichlet boundary. Read-only ONELAB
values identify directions where this extent cap is active. Near-DC
conductance signs are reported without an accuracy gate under this extent limit.
The original layer thickness adds to the capped stretched increment.

## Loss-aware PML intervals

`PmlSideLayers`, `PmlTopLayers` and `PmlBottomLayers` (default 16) supplies minimum normal interval counts, rather than fixed
counts. With ``X_d=L_d(1+A_d/4-j\eta_d A_d/4)`` for the layer alone and
``E_{m,d}=\Re(q_mX_d)``, the native bound is

```math
\Phi_{m,d}=|q_m||X_d|
\begin{cases}\min(1,T/E_{m,d}),&E_{m,d}>0,\\1,&E_{m,d}\leq0,\end{cases}
\qquad
N_d=\max\left(N_{\min,d},
\left\lceil\frac{\max_{m\in M_d}(\mathrm{PPW}_m\Phi_{m,d})}{2\pi}\right\rceil\right).
```

The fixed base coefficient is `PmlPointsPerWavelength=10`.
For positive ``a_m``,
`PPW_m=10/MeshSizeFactor*clamp(sqrt(0.1*abs(q_m)/a_m),1,3)`; for nonpositive ``a_m``, it is `30/MeshSizeFactor`.
The clamp bounds are read-only ONELAB values. The decay factor restricts the interval bound to the part of the stretched
layer before attenuation T. Exact cutoff has
``|q_m|=0`` and contributes zero to the bound.

Effective side, top and bottom counts control both transfinite divisions and
progression ``\exp(g_d/N_d)``. The default grading is `log(20)`; zero is uniform.
Inward curves reverse the ratio. PML elements default to quadrangles with
four-point quadrature; physical materials use triangles. Gmsh evaluates the interval counts from these native coefficients.

## Physical mesh and earth interface

Julia serializes the cable CAD, material ownership and terminal anchors.
Gmsh computes local disk, annulus, strand, sector and passive-media targets from
those facts and current native coefficients. Shared contour points use the
smallest adjacent target. Conductor controls set angular accuracy, elements
per skin depth, radial growth, graded depth and wall divisions. Cable boundaries
extend their local targets into passive media. Physical regions remain
triangular; the PML grid is transfinite.

The bulk target is `mesh_size_factor*max(layout_radius,L)/20`. A medium's wave
target is `min(MeshBulk,mesh_size_factor/(8*abs(q_m)))`. Its decay footprint is
at most `min(2*resolution_radius,6/a_m)` for positive ``a_m``. Excited air retains its
fine wave target inside that footprint; lossy-earth exterior growth is
described below. The wave cap constrains remote sizes and PML
tangential counts only when the physical box boundary lies inside the footprint;
remote sizes otherwise follow the bulk size and growth law. This does not remove the
fine target near a conductor or the interface.

Air is excited when any cable exterior extends into air, when the earth has
zero transverse decay rate, or when a buried cable's exterior lies within
`6/Re(q_earth)` of the interface. Otherwise, air sizes follow the bulk and remote
growth law without air wave sources or a wave cap on the physical box. Air also
leaves the top and side PML wave-resolution bounds: the top retains its interval
minimum, and the side uses earth alone. Root selection and stretch calibration
remain unchanged. The read-only ONELAB value and native observations report
whether air is excited.

When the earth interface layer is active, the earth wave-distance field
uses cable contours without the projected footprint sources. Excited-air distance
sources and all conductor targets are unchanged.
For lossy earth with a wave target below the remote target, the exterior
size grows geometrically with the native distance ``d``:

```math
h(d)=\min(\mathrm{MeshRemoteEarth},\mathrm{MeshWaveEarth}\exp(a_{\rm earth}d/2)).
```

The distance retains footprint sources when the interface layer is inactive.

An earth-only triangular BoundaryLayer field on the physical interface has
first size `MeshWaveEarth` and growth ratio 1.4. Its effective thickness is

```math
t_{\rm layer}=\min\left(\mathrm{MeshDecayEarth},
\frac12\min_{p\ {\rm buried}}(|y_p|-r_p-h_p)\right),
```

where ``r_p`` is the cable radius and ``h_p`` its cable-neighbourhood size.
Without buried cables the candidate is `MeshDecayEarth`. The layer is active
only when `MeshWaveEarth < MeshRemoteEarth` and its candidate thickness is at
least `2*MeshWaveEarth`; otherwise its effective thickness is zero.
Conductor boundary layers remain active and air and conductor surfaces are excluded
from the earth layer. Read-only native values publish the effective thickness
and the flag for a clipped or omitted layer; clipping alone does not warn.

Vertical box and side-PML edges grade from
`min(existing_last_size,MeshBulk,medium_wave_size)` with the existing native
edge law. `InterfaceRefinementFactor=1` retains the cable footprint law.
Interface mesh seeds come from one cable-centre source plus the central seeds.
They are merged within `max(1e-9,1e-12*domain_span)` metres, and also within
0.1 times the smaller local target. A native hard check rejects shorter
interface segments. Measurement endpoints are not interface mesh seeds.

## Qualification and comparison limits

The per-medium net exponent is
``\Re(q_m[D+L_d(1+A_d/4-j\eta_d A_d/4)])``. A nonpositive exponent outside
exact cutoff rejects the unsupported prescription with an explicit native
error. Deliberately capped real maps permit zero net attenuation. Exact
cutoff solves with its zero exponent published and warns that the transverse
problem is singular. Below-target attenuation and sizing-floor activation
are observations, not qualification warnings. The earth-resistive sizing
ceiling still warns that results are not qualified.

Managed execution reads frequency-aligned native observations into
`details(result).data.fem.pml_observations` and warns once per run, including
completed-run replay. The observations retain effective intervals, directional η, cap, cutoff and
floor flags, target and net exponents, the sizing-ceiling flag, and earth-layer
thickness and clipping status. Physical mesh and voltage measurement still
require convergence checks.

Exact air cutoff is singular. A prescribed Γ can be non-passive; supplying it does not guarantee a
passive forward mode. Prescriptions that produce a nonpositive net exponent
outside cutoff are unsupported. The analytical `:unified` root convention is
unchanged; an incoming-sheet reference is excluded from accuracy comparisons.

For comparison with `:unified`, its mean-field receiver represents each
conductor by a single line source. Its declared receiver range is
``|\kappa_m r_p|\lesssim0.1``. Reception at thick conductors near the interface
is outside that scope; a FEM/reference discrepancy there is not by itself a
FEM convergence estimate. Testing prescribed-Γ 2D boundary-value problems at
300 MHz does not extend the analytical engine's quasi-TEM applicability limit.
The analytical EarthModel keeps its requirement of infinite air resistivity.

## Options

Scientific selections (`physics=:helmholtz`, prescribed Γ and reductions)
belong to the formulation. Computation and detached export share one validator.
Unknown options fail. Removed expert keywords report the native override to use.

| Group | Options and defaults |
|---|---|
| Physics | `mesh_size_factor=1.0`, `domain_size_factor=2.0`, `pml_reflection=1e-3` |
| Workflow | `mesh_policy=:reuse`, `mesh_path=nothing`, `keep_run_directory=false`, `resume_run_directory=nothing`, `getdp_executable=nothing`, `output_basis=:pul`, `on_result=nothing` |
| Resources | `frequency_workers=clamp(Sys.CPU_THREADS ÷ 4, 1, 8)`, `solver_threads=1` |
| Solver | `linear_solver=:mumps` (also `:gmres`), `gmres_iterations_max=20`, `gmres_relative_tolerance=1e-12`, `gmres_absolute_tolerance=0.0` |
| Diagnostics | `plot_field_maps=false`, `verbosity=(default=0,)`, `gmsh_verbosity=2`, `getdp_verbosity=2`, `log_file=nothing`, `trace=false`, `timing=false` |
| Expert | `overrides=(;)`, a named tuple of numeric native parameters |

`mesh_size_factor` is one resolution scale: larger values coarsen conductor
angular spacing, first normal step and bulk sizes, wave and earth-layer targets,
measurement-line sizes and the PML points-per-wavelength bound. Minimum circle
segments and the minimum PML interval counts stay fixed. Factor one preserves
the default mesh. The circle segment count set by geometric area tolerance is
independent of `mesh_size_factor`. It does not change physical decay lengths or domain extent.
`domain_size_factor` multiplies the physical earth sizing length after its
resistive ceiling and prescribed-Γ adjustment. `pml_reflection` sets nominal
normal-wave attenuation. Near DC the added real stretch is limited to 100
domain half-widths and becomes purely real where the air field is quasi-static
throughout that extent; native observations report activation.

`mesh_policy=:reuse` reuses compatible meshes within a run; `:remesh`
regenerates them. Supplying `mesh_path` with `:remesh` raises `ArgumentError`.
Very large meshes may require fewer frequency workers to bound memory.

Expert overrides are validated against one declaration of types and ranges.
The native guards also reject invalid CLI and ONELAB edits. The exported
`Expert` group uses the same native parameter names and defaults.

| Native override | Default | Allowed values or meaning |
|---|---|---|
| `PmlSideThicknessFactor`, `PmlTopThicknessFactor`, `PmlBottomThicknessFactor` | `1.0` | Positive relative thickness |
| `PmlSideLayers`, `PmlTopLayers`, `PmlBottomLayers` | `16` | Positive minimum interval counts; native wave bound may raise them |
| `PmlSideGrading`, `PmlTopGrading`, `PmlBottomGrading` | `log(20)` | Nonnegative grading exponents |
| `ExteriorMeshSizeFactor` | `1.0` | Remote-buffer factor, at least one |
| `InterfaceRefinementFactor` | `1.0` | Cable-projected interface footprint factor, at least one |
| `ConductorGeometryTolerance` | `1e-3` | Positive relative circular area error |
| `ConductorSkinDepthElements` | `6.0` | Positive normal divisions per conductor skin depth |
| `ConductorMeshGrowth` | `sqrt(1.25)` | Normal growth, at least one |
| `ConductorSkinDepths` | `5.0` | Positive graded depth, limited by conductor width |
| `ConductorThicknessElements` | `4` | Positive integer wall divisions |
| `VolumeQuadrature` | `12` | Triangle points: 4, 7, 12 or 13 |
| `PhysicalVolumeQuadrature` | `3` | Triangle points: 3, 4, 7, 12 or 13 |
| `PmlQuadrangles` | `1` | 1 for quadrangles, 0 for triangles |
| `PmlQuadrature` | `4` | Quadrangle points: 4, 9 or 16 |
| `MumpsOrdering` | `0` | AMD; -1 retains solver default, other codes: 2, 3, 4, 5, 6, 7 |
| `PetscPrealloc` | `0` | 0 retains solver default; otherwise positive row allocation |
| `MumpsErrorAnalysis` | `2` | 0 off, 1 full, 2 backward errors |
| `MumpsRefinementMax` | `2` | Nonnegative refinement limit |
| `MumpsBackwardErrorTolerance` | `1e-12` | Positive backward-error target |

Three physical quadrature points integrate first-order field products exactly
for piecewise-constant materials. PML quadrature remains separate.
The MUMPS forward-error warning budget is a fixed 0.01. Diagnostics do not
trigger automatic solver switches or refinement studies. Failed and interrupted
runs retain their logs. Mesh identity includes physical inputs, Γ, mesh controls
and native sources; quadrature alone does not change the mesh prescription.

## Export a detached ONELAB project

```julia
export_data(:onelab, problem, fem;
    file_name="onelab_model/model.pro",
    options=(domain_size_factor=2.0, overrides=(PmlSideLayers=48, PmlTopLayers=48, PmlBottomLayers=48, MumpsOrdering=0),))
```

Detached `options` accepts physics, overrides, solver controls, `solver_threads`,
`plot_field_maps`, `output_basis` and native verbosity. When verbosity is omitted,
Gmsh inherits its caller setting and GetDP retains its existing stage verbosity
(3 for mesh-only, 4 for solving). Managed workers,
retention, resume, callbacks and logging options are rejected.

The exported entry, physical CAD, evaluated material values and native sources
are self-contained. Managed Julia runs use these same sources.
Use the generated README for caller-specific filenames. A native frequency
scan remeshes and solves each exported case in order. Indexed frequency, Γ
and material coefficients can be edited in ONELAB; this does not evaluate a
new Julia material law. Derived counts, layer thicknesses, sizing flags and
attenuation are read-only observations. Reset ONELAB's database to recover
file defaults after edits. Native execution does not rewrite its source files.

Primitive/reduced matrices and indexed raw columns are written below each
result directory; `completed.txt` is written last. A failed or stopped run has
no valid completion marker. The optional field maps can be viewed by Gmsh or
read with `import_data`. Their PML values are analytic continuations rather
than physical energy or loss densities.

The consolidated manual tool in `test/manual/calculations/fem_validation/`
runs FEM/reference comparisons and homogeneous closed-form checks using current
options. It is excluded from automated discovery. Reports preserve signed
components, scale-normalized errors, warnings, native timing and DOFs; numerical
agreement is separate from applicability and convergence.

## Headless computation and saved-file inspection

Computation uses the Gmsh API and standalone GetDP processes. It opens no window
and uses no ONELAB parameters. The former `ui` computation option has been removed.
Interrupting computation stops and reaps its solver processes and preserves
completed checkpoints.

Write and retain field maps when a later inspection is needed:

```julia
parameters = compute(problem, Formulation(:LineCableModelsFEM);
    options=(plot_field_maps=true, keep_run_directory=true))
run = details(parameters).data.fem.run
```

The runtime root is `.linecablemodels/fem` under the working directory captured
when computation starts. Each system has its own directory, and each run has an
exclusively created `YYYYmmdd-HHMMSS-xxxx` directory. System identifiers are
sanitized for filenames; the four-character suffix distinguishes runs created
in the same second. A retained run is self-contained:

```text
.linecablemodels/fem/<system_id>/<timestamp>-<suffix>/
    input/
    mesh/       frequency_XXXX.msh and frequency_XXXX.json
    results/    adopted matrices, column checkpoints and timing
    logs/       mesh and solver diagnostics
    maps/       only when field maps are requested
    work/       only retained after interruption or failure
    run.json
    coordinator.lock
```

Native mesh observations are folded into each mesh's JSON sidecar. There is no
cross-run mesh cache: `mesh_policy=:reuse` shares equal fingerprints within the
run and reuses valid meshes when resuming it. Repeated computations generate new
meshes unless `mesh_path` supplies a `.msh` file or a previous managed run directory.
For a supplied run, each frequency matches the fingerprints in its mesh sidecars;
matching meshes are linked or copied into the new run and validated, while unmatched
frequencies mesh normally. The new run remains intact if the source run is deleted.
A missing run directory or one without mesh sidecars is rejected. `resume_run_directory=:latest` selects the newest compatible run
for the same system. Explicit resume paths must be inside that system’s runtime
directory. Successful solver scratch is removed after adoption;
failed scratch remains available for resume. Because adopted solver scratch is
deleted, a corrupted result column is solved again on resume. A successful run is deleted unless
`keep_run_directory=true`. A failed run is retained and its path is included in
the error. Deleting its directory removes all files produced by that run.
Plotting cannot recover deleted files.

Load Gmsh and a Makie backend to inspect saved files independently:

```julia
using Gmsh, CairoMakie
using LineCableModels: plot

mesh = import_data(:msh, run.run_directory; frequency_index=1)
preview(problem.system; mesh)
plot(mesh)
plot("e_f0001_b0001.pos"; component=2, part=:real)
```

`frequency_index` follows the result's frequency order, starting at one.
Pass a saved run directory, its `mesh` directory, or a mesh-file path within
that directory. Every frequency has a `frequency_XXXX.msh` file and matching
JSON sidecar; existing retained filenames are resolved from their metadata.
Omitting the index reads the specified file; for a run directory, it reads the
last frequency in the saved order. `mesh.provenance` retains the run directory,
frequency index, frequency [Hz], terminal identifiers, evaluated earth inputs
and prescribed Γ. Older sidecars without earth inputs or Γ report those values
as `nothing`. The mesh plot and preview show the run name, index and Hz in their
overview. Direct standalone
`.msh` import remains available without run metadata:

```julia
mesh = import_data(:msh, "external_mesh.msh")
```

These operations do not resolve or launch GetDP, generate a mesh, or start a
computation. Imported `FEMMesh` and `FEMFieldMap` values are detached from Gmsh;
they remain usable after the reader releases its session. The reader preserves
caller-owned models, views, and options and finalizes only sessions it owns.
No Python runtime is required.

Native mesh import preserves node/element tags, element order and physical
names. It accepts both ASCII and binary formats supported by Gmsh. Display does
not require the physical groups that make a mesh admissible to the FEM solver.
The current renderer supports first-order planar triangles, quadrangles and
line boundaries. Higher-order and volume meshes remain importable but require
an explicit rendering capability instead of silent linearization or flattening.

For physical-group colors and click diagnostics, load an interactive backend:

```julia
using Gmsh, GLMakie
using LineCableModels: plot

mesh = import_data(:msh, "detached-study/study.msh")
p = plot(mesh; color_by=:physical, inspect=:element, backend=:gl)
# The same controls are available on the geometry preview:
g = preview(problem.system; mesh, mesh_color_by=:physical,
    mesh_inspect=:node, backend=:gl)
```

The **Color** menu selects uniform edges, physical-group colors, or entity
colors. Each distinct complete surface membership has its own category;
overlapping groups are retained together. Names are qualified by dimension and
tag, with explicit unnamed and ungrouped entries. Boundary groups remain in
the inspection record. Legends page through categories when needed.

The **Inspect** menu selects elements, nodes, or off. A click highlights the
selection and retains its diagnostics in the sidebar; dragging keeps the native
Makie navigation behavior. **Copy details** copies the full record when the
onscreen text is abbreviated. Node records show original tags, native Cartesian
coordinates, incident element tags, entities, and all physical groups. Element
records show original tags, type/order, entity, connectivity, memberships,
centroid, and geometric area/perimeter/edge lengths, or length for a line.
Both display triangles of a quadrangle select the original quadrangle. On a
shared edge, an explicitly stored boundary line takes precedence. Coincident
nodes retain separate IDs; picking chooses the lowest native tag at that point.
Zoom in when several small elements share the same screen pixels.

The defaults remain `color_by=:uniform, inspect=:none`. Preview geometry uses
outlines while mesh coloring or inspection is active. Cairo renders categories
and legends for static output; click interaction has been validated with GLMakie.
These diagnostics use only imported mesh metadata and geometry; coefficients,
constraint equations and field values are outside the mesh inspector.

The existing preview labels horizontal position as `y [m]` and height as `z [m]`.
The solver stores these as mesh x and y respectively. Field component 1 is
horizontal, 2 vertical, and 3 axial. Meshes produced here use meters; when an
external file uses millimeters, import with `coordinate_scale=1e-3`.

Native POS field views retain element-local values and output-step times.
`import_data(:pos, path)` returns one map, or a vector if the file has several
views. Select a single view with `view=1`. Owned harmonic fields include
`phasor=real,imag` in their label. For older known harmonic files use
`representation=:complex`; two unlabeled steps are otherwise independent real
samples. List-based scalar, vector and tensor values remain available in the
imported blocks; field rendering currently accepts triangle and quadrangle
samples. Model-based views require export to native list-based POS first.

Field plotting accepts `part=:real`, `:imag`, `:magnitude`, or `:phase`.
For vectors, choose `component=1`, `2`, or `3`, or omit the component with
`:magnitude` to color by the Euclidean norm of the phasor. Magnitudes are not
converted to RMS. Phase uses radians and is undefined at zero. Independent
real views select their recorded sample with `step=1`.

Field plots use a square viewport with equal spatial scaling. The field label,
units and selected component or part appear horizontally on the right. A
horizontal colorbar below the viewport follows its width when the window is
resized. Native axis attributes and `colorbar_position` / `colorbar_attributes`
remain available for presentation changes.

```julia
field = import_data(:pos, "e_f0001_b0001.pos")
plot(field; component=2, part=:real, geometry=problem.system, mesh)
plot(field; part=:magnitude, colorscale=log10)
plot(field; component=2, part=:real, arrows=true, arrow_stride=20,
    arrow_attributes=(normalize=true, lengthscale=0.005))
```

Arrows display the transverse real or imaginary vector at selected element
centers. `arrow_stride` controls element sampling. `arrow_attributes` accepts
native Makie options. The example normalizes directions and draws 5 mm arrows;
without normalization, `lengthscale` converts field-value magnitudes to plotted
lengths. These presentation choices do not change retained values.
Geometry is drawn as outlines
over field colors. Coincident field vertices retain their separate element-side
values; interpolation does not average across material interfaces. Nonfinite
samples and nonpositive logarithmic samples are not colored, and the plot status
reports their count. Original samples are preserved.

Native labels retain source normalization, physical units, scaled Helmholtz
quantities and the distinction between physical fields and PML analytic
continuations. The existing Makie shell supplies controls, legends, colorbars,
axis manipulation and SVG export through the returned `UIPlot`.

## Extension-owned definitions

```@docs
FEM.LineCableModelsFEM
FEM.LineCableModelsFEMError
FEM.FEMMesh
FEM.FEMElementBlock
FEM.FEMFieldMap
FEM.FEMFieldBlock
```

The detached ONELAB export enables automatic checks. Changing the selected
frequency or a geometry/mesh input rebuilds the geometry and generates and saves
the selected mesh without a field solve. `Mesh/Current mesh` publishes the case
index and frequency [Hz] only after mesh generation and saving succeed; its
status is `No mesh` during input checking or after a failed rebuild. Run uses the
same native geometry and mesh prescription before solving.
