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
        physics = :quasi_fw,
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
`plot_field_maps=true`, basis-specific output operations write scalar, vector and material-region field
quantities with frequency/source-specific filenames and labels. The maintained GetDP
files are captured with each run so later edits cannot change an active scan.

The former `fem_options` keyword and `LineCableModelsFEMOptions` struct have been
removed. Move their execution keys into `compute` options and `physics` into
formulation options. Benchmark calls use `reference_options` for FEM execution
controls. Supplemental run metadata is a named tuple under `details(result).data.fem.run`.

## Physics selection

The physics selector defaults to `options=(physics=:quasi_fw,)`, the coupled
prescribed-Γ Maxwell model. This is the supported field model. The string
`"quasi-fw"` and `Symbol("quasi-fw")` are also accepted. Julia requires
`Symbol("quasi-fw")` for a hyphenated symbol;
`:quasi-fw` is parsed as subtraction.

```julia
fem = Formulation(:LineCableModelsFEM; options=(physics=:quasi_fw,))
parameters = compute(problem, fem)
```

`model.pro` defaults to `Physics=1` and includes `quasi-full.pro`.
The native selector accepts `-setnumber Physics 1` and uses the resolution
`LineCableModelsFEMScan`. Select physics in formulation
`options` before computation. Physics is saved with the formulation inputs and
column checkpoints; a run cannot resume under different physics.

The quasi-full option accepts a prescribed complex ``Γ`` [1/m]:

```julia
fem = Formulation(:LineCableModelsFEM; options=(physics=:quasi_fw, Γ=0.01+0.02im))
# A vector follows problem.frequencies, including repeated frequencies.
```

The default is zero. The variables ``A_t/\Gamma`` and ``\phi/\Gamma``
retain the complete Maxwell equations, including the complex square ``\Gamma^2``.
At exactly zero the equations reduce to their regular normalized limit without
numerical division by Γ.
Finite Γ also changes series extraction: ``Z=-U/I+\Gamma^2P``.
The transverse PML wavenumbers are ``\sqrt{j\omega\mu\kappa-\Gamma^2}``;
mesh and PML coefficients are resolved for the prescribed frequency/Γ pair.
The common complex stretch ray is chosen to damp both exterior media.
At a transverse cutoff the base stretch remains finite; the zero transverse
wavenumber has algebraic truncation error, which must be checked by refinement.
One axial-current excitation supplies both responses. Its voltage extraction
includes ``j\omega\int(A_t/\Gamma)\cdot d\ell`` along physical vertical paths.
GetDP coordinates are x horizontal, y vertical and z axial, with interface y=0.
For every source column, an overhead receiving row subtracts its local scalar
surface trace and integrates from that surface projection to the conductor.
Buried rows retain the deep-earth reference and path. No surface constraint is added.
The model imposes zero scalar potential on the entire outer boundary,
at the outside of the PML, including its air side.

Every terminal has an explicit vertical measurement curve from its reference
to its lowest CAD contour vertex (lowest x breaks ties). These curves and their
reference points are physical groups in the native geometry. Gmsh embeds the
paths in their material surfaces and makes buried PML paths shared boundaries
of transfinite blocks. Every measurement element is an edge of the electric
field mesh. GetDP evaluates the `BF_Edge` trace directly and integrates it with
the four-point line rule `I2` in `integration.pro`. The positive vertical
direction fixes circulation independently of CAD edge numbering. Conductor
interiors are excluded from the path groups. No stored-field interpolation,
external triangle clipping, or separate measurement mesh is used.

The full equations, units, boundary conditions, gauge and extraction are
documented in [`LineCableModelsFEM`](@ref), with the potential-equation reference
of [Ciuprina2024](@cite). This backend's 2D longitudinal reduction is distinct
from that paper's 3D ECE implementation. The model retains the
specified finite metal conductivity in the axial problem and treat the metal
transversely as equipotential terminals.

Quasi-full additionally retains the scalar-only `Pscalar.tsv` per column,
which is gauge dependent and must not be inverted for Y. With field maps
enabled it also writes `bt_mesh`, `v_local` and `hz_scaled` per excitation.
The maps also include `material_region`, `pml_mask`,
`jt` and `jt_mesh`. `jt` contains the physical
components and `jt_mesh` the transformed flux density used in the weak form.
Values inside the PML are analytic continuations, identified by `pml_mask=1`;
they are not physical loss or energy densities. For quasi-fw, `e`, `em` and `jm` maps represent ``E_t/\Gamma`` [V], its
magnitude [V], and ``|J_t/\Gamma|`` [A/m], respectively. Multiply the
complex transverse fields by the prescribed Γ to recover physical fields.
The `b` map includes the axial component ``B_z=ΓCb`` for finite Γ.

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
both GetDP's physical-soil/PML coefficients and the skin-depth mesh rule.
`:default` and `nothing` retain static soil. Air uses its explicitly declared
static permittivity and permeability, independently of the soil law. The current
FEM geometry requires one horizontal semi-infinite soil; a non-finite conductive
skin depth is unsupported. Equivalent homogeneous-earth reductions are rejected.
An ordinary `EarthModel` supplied after an external reduction carries no history
from which FEM could detect that prior approximation.

Saved FEM formulation details contain the four consumed `selections`, their
parameters and numerical options, alongside requested and resolved records.
Custom selections retain their identities and data; saved records never
reconstruct executable methods.
The selected propagation approximation remains recorded separately.

## Field equations and matrix extraction

The coupled model uses phasors ``e^{jωt-Γz}`` and complex admittivity
``κ=σ+jωε``. Its exact normalized potentials are ``A_z=a``, ``A_t=Γb`` and
``φ=Γv``. The exterior Maxwell equations retain the longitudinal ``Γ²``
terms and reduce regularly at zero, as specified in [`LineCableModelsFEM`](@ref).

Each terminal supplies one axial current excitation. Finite metal retains
its axial ``a/u`` field and total-current constraint, with
``E_z=-jωa-u+Γ²v``. Metal surfaces are transverse equipotential terminals.

Native voltage paths run from each receiver's reference to its conductor:

```math
P_{ij}=\frac{v_i-v_{R_i}+j\omega\int_{\ell_i}b\cdot d\ell}{I_j},
\qquad Z_{ij}=-U_i/I_j+\Gamma^2P_{ij},\qquad Y=P^{-1}.
```

``P`` has units Ω m, ``Z`` Ω/m and ``Y`` S/m. The complex voltage matrix is
assembled before reductions and full inversion. Its scalar and inductive
contributions use the same axial excitation and voltage reference. The
scalar diagnostic alone is gauge dependent and cannot replace this voltage.

## Execution model

One call to `compute` resolves a two-dimensional Gmsh mesh for each
frequency/Γ case. The physical air/earth rectangle has half-width
``L=\max(R_{\mathrm{layout}},m\delta)`` and extends from ``y=-L`` to ``y=L``,
with conductive soil skin depth ``\delta=\sqrt{\rho/(\pi f\mu)}``.
`domain_skin_depths=2.0` sets ``m``. A finite Cartesian PML surrounds this box;
ordinary volume Jacobians apply throughout. Its thickness defaults to ``L``
on each side, above and below. Set `pml_thickness` to one positive length in
metres or `(side, top, bottom)`; `pml_thickness_factor` multiplies the resolved
lengths. `pml_layers=128` sets the number of intervals normal to each PML;
`pml_layers=(side, top, bottom)` sets independent directional counts. Cartesian
patches align their boundaries with the stretch onsets and the air/earth
interface; corners inherit both adjacent counts. With two triangles per cell,
the four corners contain ``4N_s(N_t+N_b)`` triangles. The default
`pml_element_family=:triangle` retains this construction. Select
`:quadrangle` to use Gmsh's native transfinite recombination on the PML
surfaces alone. It removes the cell diagonals, preserving the graded boundary
nodes, physical-region triangles and voltage-path edges. The same fixed grid
then has half as many PML elements. This is an element-family choice, not
automatic coarsening or an accuracy claim.

`pml_grading=(192/191)*log(1536)` sets a fixed dimensionless exponent, also
accepting a `(side, top, bottom)` tuple. Nodes follow
``d_i=t\operatorname{expm1}(gi/N)/\operatorname{expm1}(g)`` from the physical
interface towards the wall. Gmsh receives ``N+1`` nodes with progression
``\exp(g/N)``. Zero grading gives uniform spacing; one interval uses only the
endpoints. Exponents so small that the native ratio rounds to one also give
uniform spacing. Unrepresentable node spacing is a construction error.
The fixed default preserves the former normalized distribution at 192 intervals.
Other counts, including the unchanged API default of 128, now sample that fixed
shape instead of changing the exponent with the count. Grading does not change
PML thickness or absorption strength.

Alternatively, prescribe a mesh from physical scales:

```julia
mesh = (
    domain_skin_depths = 24.0,
    pml_resolution = (interpolation_cells=72, coefficient_change=0.12),
    mesh_size_factor = 3.0,
    exterior_mesh_size_factor = 8.0,
)
result = compute(problem, fem; options=mesh)
export_data(:onelab, problem, fem; file_name="detached/study.pro", mesh_options=mesh)
```

`pml_resolution` is optional and mutually exclusive with explicitly supplied
`pml_layers` or `pml_grading`. Its interpolation density depends on the material
Γ-dependent transverse propagation constants, PML thickness, distance from conductors, and the cubic
stretch. Propagating, grazing and evanescent scalar modes contribute to this
prescribed density. A second density, ``|\partial_u\log s|/\eta``, resolves
transformed-coefficient variation, where ``u=d/t`` and
`coefficient_change` sets ``\eta``. Larger `interpolation_cells` or smaller
`coefficient_change` generally prescribe more cells. These are numerical
resolution controls, not accuracy guarantees for ``G`` or other observables.
Geometric strip fitting approximates the density; `coefficient_change` is not
a strict bound on every realized cell's coefficient change.

The constructor integrates that density once and creates native geometric
strips at six equal density increments and the physical stretch transition
``b u^3=1``. It does not run a trial FEM problem, compare solutions, reject a
scientific result or refine after solving. A positive conductor-to-PML
clearance is required to define the modal length scale. Corners and buried
voltage paths share the same strip endpoints and native nodes. Managed
meshing and detached export consume the same resolved strips; interval counts,
endpoints and ratios enter mesh/resume identity and mesh metadata.

For the structured physical-density prescription at the retained
24-skin-depth domain, the earlier qualification preserved all 363
conductance signs in the 79-frequency two-wire, three-wire, screen, tube and
sector selection. The seven two-wire scans took 561.6 s versus 889.2 s with
144/144/96 intervals (one batch each, four frequency workers, one solver thread).
First-call compilation was 28.7 s versus 26.4 s. This is a measured 36.8%
elapsed reduction, not a general timing guarantee. Magnitude errors remain:
the worst two-wire G relative error was 76.9%, and one screen mutual G changed
by 24.8% relative to its retained reference (9.8e-15 S/m absolute). Fixed
controls remain available. See the
[qualification record](../../test/manual/calculations/fem_pml_physical_mesh.md)
for signed component comparisons and scope. Two-skin-depth domains were tested
with both 24-skin-depth and 2-skin-depth PML thicknesses at 0.1 and 21.54 Hz.
Both choices gave wrong signs for all four G entries at both frequencies.
The qualified manual preset therefore retains the 24-skin-depth domain.

`pml_reflection=1e-10` specifies the continuous normal-wave round-trip amplitude
used to calculate the cubic stretch strength, not a discretization error bound.
The profile is ``s=1+(1-j\eta)b(d/t)^3``, with ``\eta=1`` at Γ=0.
For finite Γ, the common positive ray slope ``\eta`` is reduced if needed so
``\operatorname{Re}[q_m(1-j\eta)]>0`` in both media; the strengths use these
Γ-dependent decay rates. Its real part stretches distance and
attenuates evanescent modes; its imaginary sign corresponds to ``\exp(+j\omega t)``. Side stretching is shared across
the air/earth interface; bottom stretching uses the soil propagation constant.
Conductivity, permittivity and permeability are all transformed consistently.

For a domain study use, for example,
`compute(problem, fem; options=(domain_skin_depths=1.5, pml_thickness=10.0))`.
A fixed PML thickness isolates the physical-domain change. Changing only
`domain_skin_depths` also changes the default PML thickness. Physical mesh-size
targets are independent of that option. `mesh_size_factor` scales those targets,
while `pml_layers` controls normal PML resolution. `volume_quadrature` selects
4, 7, 12 (default), or 13 triangle integration points.
`physical_volume_quadrature=nothing` inherits this rule; an explicit value
of 3, 4, 7, 12 or 13 changes only physical-material triangles.
Triangular PML continues to use `volume_quadrature`. Quadrilateral PML uses
`pml_quadrature=9` Gauss–Legendre points, with 4 and 16 also available.
The line-integral rule is unchanged. Mesh-affecting controls
enter cache identity; all controls and captured sources enter resume identity.
These defaults require convergence checks for the requested material and spectrum.

For explicit directional controls, use for example
`options=(pml_layers=(144,144,96), pml_grading=(192/191)*log(1536))` in `compute`, or
the same tuple as `mesh_options` in `export_data(:onelab, ...)`. The generated
detached `.geo` contains the same oriented transfinite constraints. Counts and
grading enter mesh and resume identity. An explicit `mesh_path` retains its
existing caller-owned mesh; these controls do not remesh that override.
An ordinary computation uses the requested controls once and returns its result;
it does not compare discretizations or choose a scientific operating point.

For example, these explicit settings are shared by both execution paths:

```julia
discretization = (pml_element_family=:quadrangle,
    physical_volume_quadrature=3, pml_quadrature=9)
result = compute(problem, fem; options=discretization)
entry = export_data(:onelab, problem, fem;
    file_name="detached/study.pro", mesh_options=discretization)
```

The detached UI exposes the PML element family and integration rules under
`Mesh` and `Numerics`. Remesh after changing the element family. Supplied
meshes must already have the requested PML element type. Integration rules
change assembly without changing mesh identity. The native implementation uses
[Gmsh transfinite recombination](https://gmsh.info/doc/texinfo/gmsh.html#t6) and
[GetDP integration criteria and rules](https://getdp.info/doc/texinfo/getdp.html#Integration).

Harmonic assembly combines conductivity and permittivity before integration:
``j\omega\sigma-\omega^2\epsilon=j\omega(\sigma+j\omega\epsilon)``.
The ``(a,a)``, ``(u_r,a)``, ``(a,u_r)`` and ``(u_r,u_r)`` trial/test pairs
retain their coefficients and supports: conductivity vanishes outside the
lossy regions, and ``u_r`` has support only in finite metal. The prescribed-Γ
couplings, terminal constraints and voltage extraction are unchanged. This is
one consolidated harmonic implementation; no legacy solver branch is retained.

The manual two-wire runner exposes these options directly. Its explicit
`--compare-native --save` mode uses the ten-frequency scan with the retained
reference mesh factors 3/8, compares the frozen pre-consolidation assembly
against the implemented choices, and saves signed R/X/G/B, absolute G
differences and timings. `--review-native` reads those results into interactive
GLMakie plots without invoking FEM. Comparisons report discrepancies without
deciding scientific acceptance or promoting defaults.

The aerial two-wire conductance study exposed low-frequency sign errors with
96 intervals in every direction. Its fixed `(144,144,96)` prescription, together
with `domain_skin_depths=24` and the retained physical/conductor targets, passed
the complete 0.1 Hz–1 MHz radius/resistivity study and the selected buried/mixed
and cable-shape checks. This establishes the observed signs for those fixtures,
not a general accuracy bound: the smallest conductance magnitude still differed
by up to 85% from its analytical reference. Exact controls, component errors,
costs and detached preservation are recorded in
`test/manual/calculations/fem_pml_conductance_cost.md`.

The final, highest-frequency mesh is displayed; every frequency-specific mesh
is reused by all of that frequency's terminal excitations. Cable topology is
constructed once: successive meshes retain its vertices and surfaces while
updating the exterior rectangle, PML patches and mesh-size fields. Only one native geometry/mesh is retained in memory. After
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
arithmetic. The mixed system is diagonally equilibrated by PETSc before
factorization and restored afterwards, preserving physical residuals and
subsequent source right-hand sides. This is needed for the large coefficient
range of the low-frequency PML; it does not change the field equations.

`mumps_ordering=nothing` and `petsc_prealloc=nothing` retain the native solver
defaults. To request AMD ordering, set `mumps_ordering=0`. Other native MUMPS
codes are 2 (AMF), 3 (Scotch), 4 (PORD), 5 (METIS), 6 (QAMD), and 7 (automatic);
their availability depends on the solver build. `petsc_prealloc=256`, for example,
requests space for 256 sparse entries per ordinary row. A larger estimate may
reduce allocation work but consumes more memory. These are prescribed execution
controls, with no automatic ordering search or retry. Compare raw results and
timings for the intended workload before choosing them.

Every attempt has a separate working directory, solver prefix, raw columns,
maps and log. Inputs are explicit command arguments and data files; the solvers
do not connect to ONELAB. One Julia coordinator owns Gmsh, logs progress,
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
runtime-directory creation, Gmsh initialization, or meshing. The preflight
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
| Frequencies | `LineParametersProblem.frequencies` | Retained in immutable input arrays; each isolated GetDP job receives its exact physical frequency and matching domain/PML parameters directly |
| Temperature | `LineParametersProblem.temperature` | Prescribed input to the selected temperature law; the default uses material `T0` and `alpha` |
| Earth material | `LineParametersProblem.earth_props` | Declared air plus one horizontal soil half-space; the soil law is evaluated per frequency |
| Optional environment declaration | `LineCableSystem.environment` | `nothing` and `EarthModel` are accepted; other declarations produce a typed unsupported-feature error |
| Line length and output basis | `LineCableSystem.line_length` and shared `compute` options | Per-unit-length is the default; total basis scales Z and Y by line length |
| Propagation constant | prescribed `options.Γ` for `physics=:quasi_fw`; exact normalized zero limit | Scalar or frequency-aligned vector; problem-level `Gamma` remains unsupported |
| Mesh resolution | local characteristic lengths derived from each resolved solid, tube, strand, foil, and passive region; per-frequency earth skin depth controls the exterior domain, and air/soil propagation scales constrain surrounding-medium resolution | Thin internal features remain local and cannot refine unrelated layers or the earth domain |

Disks, ellipses, and cable sectors retain exact Gmsh circle/ellipse arcs;
rectangles and schema polygons retain exact line segments. Annuli, conformal
sector shells, enclosure differences, and
assembly boundaries use shared oriented loops. Circular boundaries are
pre-segmented at sector endpoints and circle contacts, so adjacent materials
reuse the same curve and tangent strands reuse the same point. Compacted strand
polygons are used unchanged. Touching hole boundaries are partitioned into
connected filler faces; metal-metal seams are excluded from filler boundaries.
Annular wire-ring compartments use that same face tracing: a wire-wire point
contact separates the inner and outer filler lobes into distinct CAD faces,
without adding clearance or changing any material boundary. Separated wires
retain their connecting filler gap. All filler faces retain their declared
material. Equal evaluated material laws
share a physical material group, independently of geometric strand identity and
electrical terminal groups. A shared
material interface takes the smaller of its two local characteristic lengths.
Thin internal foils and strands do not export their size to the cable/earth
boundary. One `Distance`/`Threshold` field per actual cable exterior grows
the element-size target from that exterior layer's size with distance, using
a prescribed slope of 0.2. The central bulk cap is
``R_{\mathrm{resolution}}/20``, where
``R_{\mathrm{resolution}}=\max(R_{\mathrm{layout}},\delta)``.
This resolution scale is independent of `domain_skin_depths`;
`mesh_size_factor` scales the bulk cap and the local cable sizes.
`exterior_mesh_size_factor=1.0` retains that bulk cap throughout the domain.
Larger values permit gradual coarsening beyond ``2R_{\mathrm{resolution}}``
without multiplying conductor or central bulk targets. The remote air cap also
retains the phase-resolution bound. PML tangential divisions use these remote
sizes; the side edges are graded towards the air/soil interface. This control
does not change the physical-domain extent, PML stretch, or `pml_layers`.
Additional fields restrict the surrounding-medium size to
``h\leq 1/(8|q|)``, where ``q=\sqrt{j\omega\mu\kappa-\Gamma^2}``, within six attenuation
lengths of cable exteriors and the air/soil interface, capped at
``2R_{\mathrm{resolution}}``. The bound resolves both decay and phase;
it transitions back to the surrounding-medium size beyond that distance.
Gmsh `Restrict` fields apply each bound to its own air or soil surfaces, and
`Min` combines overlapping fields. In both air and soil, the interface
contribution is localized to each cable's projection. A cable centred at
``(x_c,y_c)`` with outer radius ``r`` contributes
an interface segment of halfwidth ``a(|y_c|+r)``, where
`interface_refinement_factor=1.0` supplies ``a``. Increase this factor (at least
one) to retain more interface refinement. Distance to these finite segments and
to the actual cable contours controls the existing size transition. This uses
native Gmsh `MathEval`, `Distance`, `Min` and `Threshold` fields.

Constant size fields remain unchanged. The local sizes and decay distances use
``q=\sqrt{j\omega\mu\kappa-\Gamma^2}`` at the prescribed frequency and complex Γ.
The footprint prescription does not inspect solutions or trigger validation,
fallback, or refinement loops. It does not assert that a distant interface
field is zero; resolution remains a prescribed numerical choice.
Conductor targets, domain/PML construction and voltage-path divisions are
independent of `interface_refinement_factor`.

Qualification covered three two-wire placements at 0.1 Hz, 1 kHz and 1 MHz for
Γ=0, and both endpoints for Γ=0.99γearth. Maximum individual R/X/G/B changes were
1.624% at Γ=0 and 0.235% at finite Γ, retaining signs. These compare independent
meshes with the same equations; they are not analytical-error bounds. The
finite-Γ comparisons use the exact metal-drive substitution described below
on both meshes. The original reference files and their unstable coefficients
remain separately recorded. At mixed 1 MHz a pre-existing 2.0655% difference
between independently constructed managed and detached meshes remains in weak
X21; the managed before/after localization difference is 0.000155%. Agreement
on a shared mesh does not establish interchangeability of independent meshes.
Full evidence and measured costs are in
`test/manual/calculations/fem_local_interface_mesh/` and
`.linecablemodels/fem/local-interface-mesh/` in the development checkout.

In finite metal the internal unknown is ``w_i=u_i-\Gamma^2v_i``, so the axial
field is ``E_z=-j\omega a-w_i``. The physical drive is recovered as
``u_i=w_i+\Gamma^2v_i`` for series-coefficient extraction. This exact substitution
cancels identical metal basis contributions before assembly. It prevents the
large drive/potential cancellation that caused finite-Γ mutual coefficients to
depend strongly on remeshing and LU ordering. Exterior equations, gauges, PML,
and native voltage integration are unchanged. At Γ=0 the substitution is the
identity. Analytical discrepancies of some weak mutual coefficients remain;
the correction does not certify physical accuracy for arbitrary cases.

The PML alone uses structured graded
quadrilateral patches split into triangles, with matching air/soil interfaces.
The physical domain retains its geometry-driven unstructured mesh. The adapter
rejects an incomplete area partition before starting Gmsh and rejects any
internal material curve without exactly two adjacent surfaces after
synchronization. After meshing, boundary-edge incidence and material coverage
are checked before a mesh can be cached or passed to GetDP. Nonempty but partial
material meshes are rejected.

Round conductors, explicit round screen wires, annular metal walls and cable sectors use
local conductor constraints. These are independent of exterior coarsening:

| Computation option | Default | Prescribed construction |
|---|---:|---|
| `conductor_geometry_tolerance` | `1e-3` | Relative circle-area target; the default requests 96 segments per full circle |
| `conductor_skin_depth_elements` | `6.0` | Maximum first normal spacing of one sixth of the conductor skin depth |
| `conductor_mesh_growth` | `sqrt(1.25)` | Ratio between successive normal layer widths; at least one |
| `conductor_skin_depths` | `5.0` | Graded depth in skin depths, limited to 80% of a disk radius or 45% of a wall thickness on each face |
| `conductor_thickness_elements` | `4` | Minimum prescribed bulk divisions through an annular wall |

The minimum circle count uses dyadic multiples of 12 and the polygon area-deficit bound
``2\pi^2/(3N^2)``. A 3 mm wire and an 85 mm wire therefore receive the same
default angular resolution. The existing local characteristic length also caps
each round conductor's edge length; a thin foil can therefore have more edges
than a solid wire at the same geometry tolerance. `mesh_size_factor` scales
this local target, while `conductor_geometry_tolerance` retains the angular
fidelity bound. Normal spacing uses each region's selected complex
admittivity and permeability at the current frequency:
``\delta=1/\operatorname{Re}\sqrt{j\omega\mu(\sigma+j\omega\epsilon)}``.
Gmsh `BoundaryLayer` fields apply only inside the owning conductor. They grade
both faces of a tube and generate triangles, with a bounded unstructured bulk
fill. When skin depth exceeds the section width, bulk sizing resolves the
section without normal layers. Weakly attenuating media retain a bulk phase
bound. Shared conductor/material interfaces remain conforming.

Inside passive cable media, a native Gmsh `Extend` field carries the actual
boundary-edge resolution into the dielectric, growing toward the cable's outer
mesh-size target. Its `Restrict` field confines it to the owning cable's passive
surfaces. This resolves the transition beside small screen wires and thin foil
without enabling global boundary-size extension in air, soil or the PML. Detached
exports retain the same fields; `MeshSizeFactor` scales the local edge cap and the
extension distance together. These prescribed sizes do not certify shunt accuracy.

Convex `Sector` conductors use native triangular transfinite strips joining the
unchanged physical contour to a copy scaled by 0.1 about the exact centroid.
The inner patch is unstructured. The longest spoke determines the common radial
count from the prescribed first step and growth; the strips extend to the core
instead of ending after `conductor_skin_depths`. Tangential spacing uses exact
arc radii and a straight-side cap of `r_back/10` at the default tolerance.
Sectors use twice the round angular count (192 at the default); tightening the
geometry tolerance also tightens the straight-side cap. The core size is at
most `r_back/30`. Internal partition edges belong to the same metal and are
excluded from the terminal contour. No copied or averaged voltage measurement
is introduced. This construction does not extend to arbitrary concave shapes,
sector-shaped metal shells or touching metal assemblies.

Pass these options to `compute` or to `export_data(:onelab, ...; mesh_options)`.
They prescribe geometry and mesh sizes, not scientific error guarantees. The
engine does not compare reference results, reject results for scientific error,
retry, or tighten these controls automatically.

Rectangular stranded cores supply their occupied disk boundary directly from
physical resolution. FEM uses the same boundary as preview, analytical
flattening and subsequent layers. Complete bounded formations are recognized
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

New meshes use binary MSH 4.1 in both managed and detached execution. This
preserves node coordinates, connectivity and physical groups while reducing
native read/write time. Existing ASCII meshes remain readable and eligible for
reuse; the storage format does not change the geometric cache identity.

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
exterior mesh size, conductor controls, the physical mesh frequency,
transformation radii, growth law, and Gmsh version.

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
Headless execution does not merge them. Detached ONELAB exports merge selected
maps after their numerical columns validate. Map paths are retained in result details
only when the run directory is retained.

Maps `az`, `b`, `bm`, `ez`, `jz`, and `rhoj2` describe the axial 1 A drive.
The `jz` map contains total axial current, including displacement current in
air and lossless dielectrics; only `rhoj2` is restricted to dissipative regions.
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

The retained run contains frequency-specific meshes under `mesh/` and native
field maps under `maps/`. Successful runs are otherwise removed according to
`keep_run_directory`; plotting cannot recover deleted files.

Load Gmsh and a Makie backend to inspect saved files independently:

```julia
using Gmsh, CairoMakie
using LineCableModels: plot

mesh = import_data(:msh, "frequency_0001.msh")
preview(problem.system; mesh)
plot(mesh)
plot("e_f0001_b0001.pos"; component=2, part=:real)
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

Native labels retain source normalization, physical units, scaled quasi-fw
quantities and the distinction between physical fields and PML analytic
continuations. The existing Makie shell supplies controls, legends, colorbars,
axis manipulation and SVG export through the returned `UIPlot`.

Native voltage extraction is checked with
`julia --project=test --startup-file=no test/runtests.jl fem_native_measurements fem_quasi_full`.
The fixtures check complex line orientation, conforming field traces, piecewise
field integration, source normalization and excitation-dependent extraction without
field-map files. The historical contour-averaged study no longer has a solve
entrypoint. Its saved `matrices.csv` can still be plotted with
`julia --project=test test/manual/calculations/run_fem_voltage_reference.jl --plot STUDY_DIRECTORY`;
those plots retain the original reference and averaging conventions.

## Export a detached ONELAB project

Load `Gmsh` and export the existing problem and FEM formulation:

```julia
entry = export_data(:onelab, problem, fem;
    file_name="detached-study/study.pro", mesh_options=(pml_layers=128,))
```

The system convenience method accepts `earth_props`, `frequencies` and
`temperature` and constructs the same problem before export:

```julia
entry = export_data(:onelab, system, fem;
    earth_props=earth, frequencies=[50.0, 10000.0], temperature=20.0,
    file_name="detached-study/study.pro")
```

Export builds native geometry without meshing, executing GetDP or opening a
window. The returned absolute `.pro` path is the ONELAB entrypoint. Its directory
contains native Gmsh/GetDP sources and usage instructions. After export it uses
only Gmsh/GetDP/ONELAB, with no Julia or Python runtime. Geometry and evaluated
material frequency cases are fixed; choose an exported case in ONELAB.
Each case also selects its prescribed Γ and matching PML coefficients.
Re-export after changing Γ so the field operator, mesh and PML remain consistent.
Export accepts `solver_options=(mumps_ordering=0, petsc_prealloc=256)`
separately from `mesh_options`. The detached project's Numerics controls expose
the same settings. Its ordering value −1 and allocation value 0 retain native
defaults; Julia uses `nothing` for these choices.

| File | Model information or responsibility |
| --- | --- |
| `study_data.pro` | Named material coefficients and units; frequency cases; terminal identities, tags and connections; source amplitudes and normalization; interactive controls |
| `study.geo`, `geometry/*.geo` | Native geometry, physical groups, voltage paths, mesh fields and transfinite PML constraints |
| `formulations/quasi-full.pro` | Regions/domains, constraints, spaces, field equations, voltage measurements and maps |
| `formulations/materials.pro`, `pml.pro` | Constitutive assignments and coordinate stretching |
| `formulations/line-parameters.pro` | Native bundle/Kron reduction, transposition and arbitrary-size `P Y = I` solves |
| `formulations/onelab.pro` | Selected-case setup, complete native resolution, tables and ONELAB publication |
| `views.geo` | Native cleanup of this project's derived views before checks and repeated solves |

Native model edits are consumed on the next run and source files are never
rewritten during execution. Only `DefineConstant` controls use ONELAB's
persisted interactive values; reset the database to restore file defaults.
When relocating an already opened project, remove its derived `.db` if it
contains stale session paths. The named material and connection assignments
remain file-authoritative. Source amplitudes also normalize Z and P: scaling
them changes fields while preserving normalized coefficients.

Open `study.pro` in Gmsh. Its normal GetDP solver registration is the only
executable setting. Select the frequency, formulation, basis and mesh settings,
then Run. Basis zero computes full matrices; a positive terminal gives diagnostic
columns. Run action also supports mesh-only execution. Check parses without
solving; native Stop controls GetDP. The GetDP thread control sets its `-nt`
option, separately from any linked numerical library's threading policy.

Enable **Run frequency scan** for one Run to mesh and solve every exported
frequency in order. ONELAB runs the native Gmsh/GetDP sequence for each case;
the selected formulation, basis and other settings apply throughout. The manual
frequency dropdown is hidden while scanning and keeps its previous selection
when scanning is disabled. Each scan starts at the first case, including after
Stop. Tables and maps remain under each frequency's result directory; the Results
panel, displayed fields and working mesh correspond to the last processed case.
Mesh-only scans visit all cases and retain the last mesh.

```sh
gmsh study.geo -setnumber BuildMesh 1 -0
getdp study.pro -msh study.msh -solve LineCableModelsFEM
```

Use the actual entry filename. Select another frequency with
`-setnumber FrequencyIndex 2` on both commands. GetDP accepts
`-setnumber Physics 1` (quasi-full-wave, the default),
`-setnumber BasisTerminal 1`, and `-setnumber PlotFieldMaps 0`.
The `-msh` option also accepts an existing compatible mesh. The Mesh panel exposes
`MeshSizeFactor` and `ExteriorMeshSizeFactor` with the same values and meanings as
Julia's `mesh_size_factor` and `exterior_mesh_size_factor`. They recompute the
size fields, graded measurement paths and PML tangential divisions when meshing.
`InterfaceRefinementFactor` matches Julia's `interface_refinement_factor` and
widens the local interface footprints in both air and soil. The native fields
use the edited physical/exterior size factors for local targets and remote caps.
`MeshRefinements` applies uniform refinement. Measurement curves conform to the
field mesh and refine with it. PML normal strips retain their exported prescription.
The native conductor controls appear under **Mesh**. They recompute conductor
sizes from the selected frequency and editable material arrays. `MeshSizeFactor`
changes background targets without coarsening the prescribed conductor layers
or circle resolution. Editing a frequency/material does not reconstruct the
exported exterior/PML geometry; re-export when those domain dimensions must
change.

The two direct commands above run one selected case. ONELAB also supports a
native batch scan using its configured GetDP executable:

```sh
gmsh study.pro -setnumber RunFrequencyScan 1 -run
```

Derived files are under `results/fNNNN-physics-bNNNN/`: native indexed columns
in `raw/jobs/`, `.pos` fields in `maps/`, and named tables in `matrices/`.
Primitive and reduced P are inverse-admittance coefficients in ohm m;
GetDP solves `P Y = I` with no extra `j omega` factor. Total output multiplies Z
and Y by line length. The completion marker is written last and invalidated on
Check/Run; partial files from failed or stopped calculations are not completed
results. Gmsh views and the separate `import_data(:msh, ...)` /
`import_data(:pos, ...)` Makie inspection workflow consume these native artifacts.

Voltage curves are conforming field edges; refining them changes the field
mesh too. Agreement between two execution routes is not a physical accuracy
certificate. See the exported README for the reference and PML conventions.

`overwrite=true` replaces recorded export files and removes obsolete owned
files while preserving unrelated files. Export itself does not execute a solver
or open a GUI. `test/manual/fem/export_onelab_toy.jl` emits a small mixed-wire toy;
its eight PML layers demonstrate execution and preservation, rather than
engineering convergence.
