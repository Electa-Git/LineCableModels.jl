function description(::Type{LineCableModelsFEM}, slot::Union{Val{:reduce_bundle},Val{:kron_reduction},Val{:ideal_transposition}}, value::Bool; compact::Bool=false)
    description(LineParametersFormulation, slot, value; compact)
end

"""Describe the active FEM field equations without executing a field solve."""
description(::Type{LineCableModelsFEM},::Val{:physics},value::Symbol;compact::Bool=false) =
    description(LineCableModelsFEM,Val(:physics),Val(value);compact)
description(::Type{LineCableModelsFEM},::Val{:physics},::Val{:helmholtz};compact::Bool=false) =
    "Helmholtz"
description(::Type{LineCableModelsFEM},::Val{:Γ},value::Number;compact::Bool=false) =
    "Γ="*string(value)*" m⁻¹"
description(::Type{LineCableModelsFEM},::Val{:Γ},value::AbstractVector;compact::Bool=false) =
    "Γ=["*join(value,", ")*"] m⁻¹ (frequency order)"

function formulation_options(::Type{LineCableModelsFEM}, record::FormulationOptions)::FormulationOptions
    options = record.data
    physics = get(options, :physics, :helmholtz)
    physics isa Union{Symbol, AbstractString} || throw(ArgumentError(
        "physics must be :helmholtz"))
    physics = Symbol(physics)
    physics === :helmholtz || throw(ArgumentError(
        "physics must be :helmholtz (also accepted as the string \"helmholtz\")"))
    Γ = get(options, :Γ, 0)
    Γ isa Union{Number,AbstractVector} || throw(ArgumentError(
        "FEM Γ must be a finite scalar or a frequency-aligned vector [1/m]"))
    values = Γ isa Number ? (Γ,) : Γ
    !isempty(values) && all(x -> x isa Number && !(x isa Bool) && isfinite(x), values) ||
        throw(ArgumentError("FEM Γ must be a finite scalar or nonempty finite vector [1/m]"))
    reductions = formulation_options(LineParametersFormulation,
        FormulationOptions(; (key => value for (key, value) in pairs(options) if key ∉ (:physics,:Γ))...))
    return FormulationOptions(; reductions.data..., physics, Γ=Γ isa AbstractVector ? copy(Γ) : Γ)
end

"""
$(TYPEDSIGNATURES)

Validate Gmsh/GetDP execution controls supplied to `compute(...; options)`.
The field model is selected separately by `formulation_options(LineCableModelsFEM, ...)`.

# Arguments

- `owner`: The `LineCableModelsFEM` computation owner.
- `options`: Caller-supplied named tuple. Supported keys and defaults are:
  - `plot_field_maps=false`: Emit spatial field maps for every solve.
  - `mesh_policy=:reuse`: Reuse compatible meshes. `:remesh` regenerates them.
  - `mesh_path=nothing`: Optional existing `.msh` path.
  - `domain_size_factor=2.0`: Multiplier of the physical earth sizing length \\[dimensionless\\].
    `min(1/real(q0_e),2pi/abs(q0_e))` \\[m\\], subject to the native ceiling
    and prescribed-Γ adjustment below. The half-width is `max(layout_radius, m*L*abs(q0_e)/max(abs(q0_e),abs(q_e)))`,
    where `m` is this option and `q0_e`, `q_e` are the earth roots \\[1/m\\]
    at zero and prescribed Γ. The sizing length \\[m\\] is
    `L=min(1/real(q0_e),2pi/abs(q0_e),sqrt(2e5/(omega*mu_earth)))`.
    The decay length is infinite when `real(q0_e)=0`. The last term is the
    reference-earth ceiling at resistivity `1e5` \\[ohm m\\]. When it binds,
    native observations mark results unqualified and the run warns.
    Must be finite and positive. Changing it retains local mesh-size targets.
  - `pml_thickness_factor=1.0`: Positive relative thickness
    \\[dimensionless\\], or a `(side, top, bottom)` tuple. Each thickness is
    the native frequency-dependent physical half-width times this factor.
  - `pml_layers=16`: Minimum positive normal interval count, or a
    `(side, top, bottom)` tuple. The native wave bound may raise each count to
    `max(Nmin, ceil(max(PPW_m*Phi_m)/(2*pi)))`. Gmsh uses one more node per curve.
    The native `PPW_m = 10*clamp(sqrt(0.1*abs(q_m)/real(q_m)),1,3)`
    uses the upper clamp when `real(q_m)=0`. Clamp bounds are read-only ONELAB
    values. The default interval floor remains 48.
    For layer length `X=L*(1+A/4-im*eta*A/4)` \\[m\\], root `q` \\[1/m\\],
    exponent `E=real(q*X)` and target `T=-log(pml_reflection)/2`,
    `Phi=abs(q)*abs(X)*min(1,T/E)` when `E>0`, otherwise `abs(q)*abs(X)`.
    Side layers use both media. Top uses air. Bottom uses earth. The constant
    10 points per wavelength is native and fixed. Exact cutoff gives `Phi=0`.
  - `pml_grading=log(20)`: Nonnegative exponent \\[dimensionless\\], or a
    `(side, top, bottom)` tuple. Native ratio is `exp(g/N)` with the effective
    interval count.  Zero is uniform.
  - `pml_reflection=1e-3`: Nominal normal-wave round-trip amplitude target
    \\[dimensionless\\], strictly between zero and one.  Not an error bound on Y.
    Near DC, the added real stretch is limited to 100 domain half-widths.
    Where the air field is quasi-static over that whole extent, the stretch
    is purely real. Read-only ONELAB values identify active directions.
  - `mesh_size_factor=1.0`: Positive physical-region mesh-size multiplier
    \\[dimensionless\\], independent of the normal PML grading.
  - `exterior_mesh_size_factor=1.0`: Maximum remote-to-central bulk element-size
    ratio \\[dimensionless\\], finite and at least one. Grading starts beyond
    twice the resolution radius.  conductor targets remain fixed. Air sizes
    remain bounded by the wavelength. Also controls tangential PML divisions.
  - `interface_refinement_factor=1.0`: Multiplier for the cable's projected
    interface footprint \\[dimensionless\\], finite and at least one. Larger
    values retain more interface refinement in both air and soil. Local wave
    sizes and transition distances retain their frequency, material and Γ dependence.
  - `volume_quadrature=12`: Points for volume integration in triangles (4, 7, 12 or 13).
  - `physical_volume_quadrature=3`: Override triangle integration in
    physical materials with 3, 4, 7, 12 or 13 points. Three points integrate
    first-order field products exactly for piecewise-constant materials. `nothing` inherits
    `volume_quadrature`. Triangular PML keeps `volume_quadrature` independently.
  - `pml_element_family=:quadrangle`: PML cells, either `:triangle` or
    `:quadrangle`. Quadrangles recombine the same native transfinite grid.
    Physical materials and voltage-path edges retain their construction.
  - `pml_quadrature=4`: Gauss-Legendre points per PML quadrangle (4, 9 or 16).
    Unused for triangular PML. These controls select prescribed quadrature rules.
  - `conductor_geometry_tolerance=1e-3`: Requested relative circle-area
    construction tolerance \\[dimensionless\\].  selects at least 96 segments
    per complete circle at the default. Sector strips use twice this angular
    count and a corresponding straight-side cap. This is not a result acceptance test.
  - `conductor_skin_depth_elements=6.0`: Reciprocal first normal step in a
    conducting region's skin depth \\[dimensionless\\].
  - `conductor_mesh_growth=sqrt(1.25)`: Normal layer-size ratio, at least one
    \\[dimensionless\\]. A value of one requests uniform layers.
  - `conductor_skin_depths=5.0`: Requested graded depth in skin depths
    \\[dimensionless\\], bounded by the disk or annular section. Sector strips
    instead extend to a fixed internal core.
  - `conductor_thickness_elements=4`: Minimum prescribed divisions through a
    tubular wall \\[dimensionless\\]. These controls grade circular conductors
    (including screen wires), annuli and convex cable sectors.
  - `keep_run_directory=false`: Retain successful run artifacts.
  - `getdp_executable=nothing`: Executable override. Otherwise resolve the
    environment override, package artifact, or unsupported-platform `PATH` fallback.
  - `gmsh_verbosity=2`, `getdp_verbosity=2`: Native message levels from 0 through 5.
  - `frequency_workers=clamp(Sys.CPU_THREADS ÷ 4, 1, 8)`: Maximum concurrent
    frequency solver processes. Very large meshes may require fewer workers
    to bound memory use.
  - `solver_threads=1`: BLAS and OpenMP threads per solver process.
  - `linear_solver=:mumps`: Native direct MUMPS or `:gmres`. GMRES starts
    from zero with fixed LU/MUMPS preconditioning and no MUMPS refinement.
  - `gmres_iterations_max=20`: Positive maximum total GMRES iteration count.
  - `gmres_relative_tolerance=1e-12`, `gmres_absolute_tolerance=0.0`: Relative
    and absolute residual targets for the PETSc-scaled system in GMRES mode.
    The relative target is positive. Zero disables the absolute threshold.
    Native convergence and the explicitly recomputed residual are reported.
  - `mumps_ordering=0`: Native MUMPS `ICNTL(7)` ordering code, defaulting to AMD:
    0 (AMD), 2 (AMF), 3 (Scotch), 4 (PORD), 5 (METIS), 6 (QAMD), or 7
    (automatic). Availability depends on the GetDP/MUMPS build. `nothing`
    leaves the solver default unchanged. No alternative ordering is retried.
  - `petsc_prealloc=nothing`: Optional positive sparse-row allocation estimate
    passed to GetDP's `-petsc_prealloc`. Larger estimates can reduce reallocations
    at the expense of memory. `nothing` retains GetDP's default.
  - `mumps_error_analysis=2`: In direct mode, native `ICNTL(11)`: 0 disables diagnostics,
    2 reports backward errors, and 1 also estimates solution sensitivity.
    Estimates concern the PETSc-scaled field system, not terminal G accuracy.
  - `mumps_refinement_max=2`: Nonnegative maximum native iterative-refinement
    corrections (`ICNTL(10)`) in direct mode. Zero disables refinement.
    Reuses the LU factors. MUMPS refinement and error-analysis controls apply to
    direct MUMPS, not its use as a GMRES preconditioner.
  - `mumps_backward_error_tolerance=1e-12`: Positive finite stopping target
    (`CNTL(2)`) for the sum of MUMPS backward errors. Managed execution warns
    when enabled diagnostics report an unmet target, without rejecting results.
  - `mumps_forward_error_tolerance=0.01`: Positive finite warning budget for
    the estimated scaled-solution forward error, used with error analysis 1.
    Detached GetDP reports native estimates and targets. Its `.pro` language
    does not expose MUMPS estimates for target-based warning comparisons.
  - `verbosity=(default=0,)`: Julia logging levels from 0 through 2.
  - `output_basis=:pul`: Per-unit-length matrices. `:total` scales by line length.
  - `trace=false`: Retain primitive matrices.
  - `timing=false`: Retain fresh complete-scan measurements in result details.
  - `on_result=nothing`: Optional callback `(problem, index, result)`.
  - `log_file=nothing`: Optional Julia log path.
  - `resume_run_directory=nothing`: Resume a compatible path or `:latest`.

# Returns

- A fixed-key [`ComputationOptions`](@ref) named tuple. Paths and integer controls
  are normalized. `output_basis` and `trace` are lowered to `Val` values.

# Errors

- `ArgumentError`: Unknown keys and invalid control values, including empty paths.
"""
function computation_options(::Type{LineCableModelsFEM}, record::ComputationOptions)::ComputationOptions
    options = record.data
    defaults = (plot_field_maps=false, mesh_policy=:reuse,
        mesh_path=nothing, domain_size_factor=2.0,
        pml_thickness_factor=1.0, pml_layers=16,
        pml_grading=log(20), pml_reflection=1e-3,
        mesh_size_factor=1.0, exterior_mesh_size_factor=1.0,
        interface_refinement_factor=1.0, volume_quadrature=12,
        physical_volume_quadrature=3, pml_element_family=:quadrangle,
        pml_quadrature=4,
        conductor_geometry_tolerance=1e-3, conductor_skin_depth_elements=6.0,
        conductor_mesh_growth=sqrt(1.25), conductor_skin_depths=5.0,
        conductor_thickness_elements=4,
        keep_run_directory=false, getdp_executable=nothing,
        gmsh_verbosity=2, getdp_verbosity=2, frequency_workers=clamp(Sys.CPU_THREADS ÷ 4, 1, 8), solver_threads=1,
        linear_solver=:mumps, gmres_iterations_max=20,
        gmres_relative_tolerance=1e-12, gmres_absolute_tolerance=0.0,
        mumps_ordering=0, petsc_prealloc=nothing,
        mumps_error_analysis=2, mumps_refinement_max=2,
        mumps_backward_error_tolerance=1e-12, mumps_forward_error_tolerance=0.01,
        log_file=nothing, resume_run_directory=nothing)
    standard_keys = (:verbosity, :output_basis, :trace, :on_result, :timing)
    unknown = setdiff(keys(options), (keys(defaults)..., standard_keys...))
    isempty(unknown) || throw(ArgumentError(
        "unknown LineCableModelsFEM computation options: $(Tuple(unknown))"))
    standard = computation_options(LineCableModelsCoaxial,
        ComputationOptions(; (key => value for (key, value) in pairs(options) if key in standard_keys)...))
    normalized = merge(defaults,
        (; (key => value for (key, value) in pairs(options) if key in keys(defaults))...))
    for name in (:plot_field_maps, :keep_run_directory)
        getproperty(normalized, name) isa Bool || throw(ArgumentError("$name must be Bool"))
    end
    for name in (:conductor_geometry_tolerance, :conductor_skin_depth_elements,
                 :conductor_skin_depths, :conductor_mesh_growth)
        value = getproperty(normalized, name)
        value isa Real && !(value isa Bool) && isfinite(value) &&
            0 < value <= floatmax(Float64) && Float64(value) > 0 ||
            throw(ArgumentError("$name must be finite, positive and Float64-representable"))
    end
    normalized.conductor_mesh_growth >= 1 ||
        throw(ArgumentError("conductor_mesh_growth must be at least one"))
    normalized.mesh_policy in (:reuse, :remesh) || throw(ArgumentError(
        "mesh_policy must be :reuse or :remesh"))
    radius_factor = normalized.domain_size_factor
    radius_factor isa Real && !(radius_factor isa Bool) &&
        isfinite(radius_factor) && 0 < radius_factor <= floatmax(Float64) &&
        Float64(radius_factor) > 0 || throw(ArgumentError(
            "domain_size_factor must be a finite positive Float64-representable number"))
    mesh_factor = normalized.mesh_size_factor
    mesh_factor isa Real && !(mesh_factor isa Bool) &&
        isfinite(mesh_factor) && 0 < mesh_factor <= floatmax(Float64) &&
        Float64(mesh_factor) > 0 || throw(ArgumentError(
            "mesh_size_factor must be a finite positive Float64-representable number"))
    exterior_factor = normalized.exterior_mesh_size_factor
    exterior_factor isa Real && !(exterior_factor isa Bool) &&
        isfinite(exterior_factor) && 1 <= exterior_factor <= floatmax(Float64) ||
        throw(ArgumentError("exterior_mesh_size_factor must be finite and at least one"))
    interface_factor = normalized.interface_refinement_factor
    interface_factor isa Real && !(interface_factor isa Bool) &&
        isfinite(interface_factor) && 1 <= interface_factor <= floatmax(Float64) ||
        throw(ArgumentError("interface_refinement_factor must be finite and at least one"))
    thickness_factor = normalized.pml_thickness_factor
    thickness_factor = thickness_factor isa Real ? (thickness_factor,thickness_factor,thickness_factor) : thickness_factor
    thickness_factor isa Tuple && length(thickness_factor) == 3 &&
        all(t -> t isa Real && !(t isa Bool) && isfinite(t) && 0 < t <= floatmax(Float64) && Float64(t) > 0,thickness_factor) ||
        throw(ArgumentError("pml_thickness_factor must be positive or a positive (side, top, bottom) tuple"))
    thickness_factor = Float64.(thickness_factor)
    reflection = normalized.pml_reflection
    reflection isa Real && !(reflection isa Bool) && 0 < reflection < 1 &&
        0 < Float64(reflection) < 1 || throw(ArgumentError(
            "pml_reflection must be strictly between zero and one"))
    layers = normalized.pml_layers
    layers = layers isa Integer ? (layers,layers,layers) : layers
    layers isa Tuple && length(layers) == 3 &&
        all(n -> n isa Integer && !(n isa Bool) && 1 <= n < typemax(Cint),layers) ||
        throw(ArgumentError("pml_layers must be a positive integer or a (side, top, bottom) tuple; each count plus one must fit Gmsh's Cint node count"))
    layers = Int.(layers)
    grading = normalized.pml_grading
    grading = grading isa Real ? (grading,grading,grading) : grading
    grading isa Tuple && length(grading) == 3 &&
        all(g -> g isa Real && !(g isa Bool) && isfinite(g) && 0 <= g <= floatmax(Float64) && (iszero(g) || Float64(g) > 0),grading) ||
        throw(ArgumentError("pml_grading must be a finite nonnegative exponent or a (side, top, bottom) tuple"))
    grading = Float64.(grading)
    normalized.volume_quadrature isa Integer &&
        normalized.volume_quadrature in (4, 7, 12, 13) || throw(ArgumentError(
            "volume_quadrature must be 4, 7, 12 or 13 triangle points"))
    physical_quadrature = normalized.physical_volume_quadrature
    physical_quadrature === nothing || physical_quadrature isa Integer &&
        physical_quadrature in (3, 4, 7, 12, 13) || throw(ArgumentError(
            "physical_volume_quadrature must be nothing or 3, 4, 7, 12 or 13 triangle points"))
    normalized.pml_element_family in (:triangle, :quadrangle) || throw(ArgumentError(
        "pml_element_family must be :triangle or :quadrangle"))
    normalized.pml_quadrature isa Integer && normalized.pml_quadrature in (4, 9, 16) ||
        throw(ArgumentError("pml_quadrature must be 4, 9 or 16 Gauss-Legendre points"))
    for name in (:mesh_path, :getdp_executable, :log_file)
        value = getproperty(normalized, name)
        value isa Union{Nothing, AbstractString} || throw(ArgumentError(
            "$name must be a path string or nothing"))
        value === nothing || !isempty(value) || throw(ArgumentError("$name cannot be empty"))
    end
    for name in (:gmsh_verbosity, :getdp_verbosity)
        value = getproperty(normalized, name)
        value isa Integer && !(value isa Bool) && value in 0:5 || throw(ArgumentError(
            "$name must be an integer from 0 through 5"))
    end
    for name in (:frequency_workers, :solver_threads, :conductor_thickness_elements)
        value = getproperty(normalized, name)
        value isa Integer && !(value isa Bool) && 1 <= value <= typemax(Int) ||
            throw(ArgumentError("$name must be a positive integer"))
    end
    resume = normalized.resume_run_directory
    ordering = normalized.mumps_ordering
    ordering === nothing || ordering isa Integer && !(ordering isa Bool) &&
        ordering in (0, 2, 3, 4, 5, 6, 7) || throw(ArgumentError(
            "mumps_ordering must be nothing or a native ordering code: 0, 2, 3, 4, 5, 6, 7"))
    prealloc = normalized.petsc_prealloc
    prealloc === nothing || prealloc isa Integer && !(prealloc isa Bool) &&
        1 <= prealloc <= typemax(Cint) || throw(ArgumentError(
            "petsc_prealloc must be nothing or a positive Cint-representable integer"))
    normalized.linear_solver in (:mumps, :gmres) ||
        throw(ArgumentError("linear_solver must be :mumps or :gmres"))
    iterations = normalized.gmres_iterations_max
    iterations isa Integer && !(iterations isa Bool) && 1 <= iterations <= typemax(Cint) ||
        throw(ArgumentError("gmres_iterations_max must be a positive Cint-representable integer"))
    analysis = normalized.mumps_error_analysis
    analysis isa Integer && !(analysis isa Bool) && analysis in (0, 1, 2) ||
        throw(ArgumentError("mumps_error_analysis must be 0 (off), 1 (full), or 2 (backward errors)"))
    refinement = normalized.mumps_refinement_max
    refinement isa Integer && !(refinement isa Bool) && 0 <= refinement <= typemax(Cint) ||
        throw(ArgumentError("mumps_refinement_max must be a nonnegative Cint-representable integer"))
    for name in (:mumps_backward_error_tolerance, :mumps_forward_error_tolerance,
            :gmres_relative_tolerance, :gmres_absolute_tolerance)
        value = getproperty(normalized, name)
        allow_zero = name === :gmres_absolute_tolerance
        value isa Real && !(value isa Bool) && isfinite(value) &&
            0 <= value <= floatmax(Float64) &&
            (Float64(value) > 0 || allow_zero && iszero(value)) ||
            throw(ArgumentError("$name must be finite, $(allow_zero ? "nonnegative" : "positive") and Float64-representable"))
    end
    resume === nothing || resume === :latest || resume isa AbstractString && !isempty(resume) ||
        throw(ArgumentError("resume_run_directory must be nothing, :latest, or a nonempty path string"))
    return ComputationOptions(; standard.data...,
        plot_field_maps=normalized.plot_field_maps,
        mesh_policy=normalized.mesh_policy,
        mesh_path=normalized.mesh_path === nothing ? nothing : String(normalized.mesh_path),
        domain_size_factor=Float64(radius_factor),
        pml_thickness_factor=thickness_factor,
        pml_layers=layers, pml_grading=grading,
        pml_reflection=Float64(reflection), mesh_size_factor=Float64(mesh_factor),
        exterior_mesh_size_factor=Float64(exterior_factor),
        interface_refinement_factor=Float64(interface_factor),
        volume_quadrature=Int(normalized.volume_quadrature),
        physical_volume_quadrature=physical_quadrature === nothing ? nothing : Int(physical_quadrature),
        pml_element_family=normalized.pml_element_family,
        pml_quadrature=Int(normalized.pml_quadrature),
        conductor_geometry_tolerance=Float64(normalized.conductor_geometry_tolerance),
        conductor_skin_depth_elements=Float64(normalized.conductor_skin_depth_elements),
        conductor_mesh_growth=Float64(normalized.conductor_mesh_growth),
        conductor_skin_depths=Float64(normalized.conductor_skin_depths),
        conductor_thickness_elements=Int(normalized.conductor_thickness_elements),
        keep_run_directory=normalized.keep_run_directory,
        getdp_executable=normalized.getdp_executable === nothing ? nothing : String(normalized.getdp_executable),
        gmsh_verbosity=Int(normalized.gmsh_verbosity), getdp_verbosity=Int(normalized.getdp_verbosity),
        frequency_workers=Int(normalized.frequency_workers), solver_threads=Int(normalized.solver_threads),
        linear_solver=normalized.linear_solver, gmres_iterations_max=Int(iterations),
        gmres_relative_tolerance=Float64(normalized.gmres_relative_tolerance),
        gmres_absolute_tolerance=Float64(normalized.gmres_absolute_tolerance),
        mumps_ordering=ordering === nothing ? nothing : Int(ordering),
        petsc_prealloc=prealloc === nothing ? nothing : Int(prealloc),
        mumps_error_analysis=Int(analysis), mumps_refinement_max=Int(refinement),
        mumps_backward_error_tolerance=Float64(normalized.mumps_backward_error_tolerance),
        mumps_forward_error_tolerance=Float64(normalized.mumps_forward_error_tolerance),
        log_file=normalized.log_file === nothing ? nothing : String(normalized.log_file),
        resume_run_directory=resume isa AbstractString ? String(resume) : resume)
end
