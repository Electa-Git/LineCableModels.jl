function formulation_options(
        ::Type{LineParametersFormulation},
        record::FormulationOptions
)::FormulationOptions
    options = record.data
    allowed = (
        :reduce_bundle,
        :kron_reduction,
        :ideal_transposition
    )
    unknown = filter(key -> key ∉ allowed, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown line-parameter formulation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge(
        (
            reduce_bundle = true,
            kron_reduction = true,
            ideal_transposition = true
        ),
        options
    )
    all(name -> getproperty(normalized, name) isa Bool,
        (:reduce_bundle, :kron_reduction, :ideal_transposition)) || throw(ArgumentError(
        "reduction and transposition options must be Bool",
    ))
    return FormulationOptions(;
        reduce_bundle = normalized.reduce_bundle,
        kron_reduction = normalized.kron_reduction,
        ideal_transposition = normalized.ideal_transposition
    )
end

function computation_options(
        ::Type{LineCableModelsCoaxial},
        record::ComputationOptions
)::ComputationOptions
    options = record.data
    allowed = (:verbosity, :output_basis, :trace, :on_result, :timing)
    unknown = filter(key -> key ∉ allowed, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown LineCableModelsCoaxial computation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge(
        (verbosity = (default = 0,), output_basis = :pul,
            trace = false, on_result = nothing, timing = false),
        options
    )
    levels = verbosity(normalized.verbosity)
    basis_value = normalized.output_basis
    basis_value in (:pul, :total) || throw(ArgumentError(
        "output_basis must be :pul or :total; got $(repr(basis_value))",
    ))
    normalized.trace isa Bool || throw(ArgumentError("trace must be Bool"))
    normalized.timing isa Bool || throw(ArgumentError("timing must be Bool"))
    return ComputationOptions(;
        verbosity = levels,
        output_basis = Val(basis_value),
        trace = Val(normalized.trace),
        on_result = normalized.on_result,
        timing = normalized.timing
    )
end


description(::Type{<:Union{LineParametersFormulation,LineCableModelsFEM}},
    ::Val{:reduce_bundle},value::Bool;compact::Bool=false) = "bundle reduction="*string(value)
description(::Type{<:Union{LineParametersFormulation,LineCableModelsFEM}},
    ::Val{:kron_reduction},value::Bool;compact::Bool=false) = "Kron reduction="*string(value)
description(::Type{<:Union{LineParametersFormulation,LineCableModelsFEM}},
    ::Val{:ideal_transposition},value::Bool;compact::Bool=false) = "ideal transposition="*string(value)

"""Describe the active FEM field equations without executing a field solve."""
description(::Type{LineCableModelsFEM},::Val{:physics},value::Symbol;compact::Bool=false) =
    description(LineCableModelsFEM,Val(:physics),Val(value);compact)
description(::Type{LineCableModelsFEM},::Val{:physics},::Val{Symbol("quasi-fw")};compact::Bool=false) =
    "quasi-full-wave"
description(::Type{LineCableModelsFEM},::Val{:Γ},value::Number;compact::Bool=false) =
    "Γ="*string(value)*" m⁻¹"
description(::Type{LineCableModelsFEM},::Val{:Γ},value::AbstractVector;compact::Bool=false) =
    "Γ=["*join(value,", ")*"] m⁻¹ (frequency order)"

function formulation_options(::Type{LineCableModelsFEM}, record::FormulationOptions)::FormulationOptions
    options = record.data
    physics = get(options, :physics, :quasi_fw)
    physics isa Union{Symbol, AbstractString} || throw(ArgumentError(
        "physics must be :quasi_fw"))
    physics = Symbol(replace(String(physics), '_' => '-'))
    physics === Symbol("quasi-fw") || throw(ArgumentError(
        "physics must be :quasi_fw (also accepted as hyphenated strings or Symbols)"))
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
  - `mesh_policy=:reuse`: Reuse compatible meshes; `:remesh` regenerates them.
  - `mesh_path=nothing`: Optional existing `.msh` path.
  - `domain_skin_depths=2.0`: Minimum physical-domain half-width in conductive
    earth skin depths \\[dimensionless\\]; the layout can require a larger box.
    Must be finite and positive. Changing it retains local mesh-size targets.
  - `pml_thickness=nothing`: Cartesian PML thickness \\[m\\]. A positive number
    applies to every side; a positive `(side, top, bottom)` tuple sets each
    direction. `nothing` uses the physical-domain half-width at each frequency.
  - `pml_thickness_factor=1.0`: Positive multiplier \\[dimensionless\\] applied
    to the resolved thicknesses; permits frequency-dependent thickness studies.
  - `pml_layers=128`: Positive interval count normal to each PML, or a
    `(side, top, bottom)` tuple of counts. Gmsh receives one more node per curve.
  - `pml_grading=(192/191)*log(1536)`: Nonnegative grading exponent
    \\[dimensionless\\], or a `(side, top, bottom)` tuple. For `N` intervals,
    the native progression ratio is `exp(g/N)`; zero gives uniform spacing.
    The fixed default reproduces the former distribution at 192 intervals.
    Other counts sample this fixed shape instead of changing its exponent.
  - `pml_resolution=nothing`: Optional physical mesh prescription, e.g.
    `(interpolation_cells=72, coefficient_change=0.12)`. The interpolation
    density resolves propagation, attenuation and source clearance; the
    coefficient density resolves logarithmic stretch variation. Both controls
    are dimensionless resolution parameters, not field-error tolerances.
    Counts and native geometric strips are calculated once per frequency.
    This option is mutually exclusive with explicit `pml_layers`/`pml_grading`.
    It performs no trial solves, scientific validation or adaptive refinement.
  - `pml_reflection=1e-10`: Nominal normal-wave round-trip amplitude target
    \\[dimensionless\\], strictly between zero and one; not an error bound on Y.
  - `mesh_size_factor=1.0`: Positive physical-region mesh-size multiplier
    \\[dimensionless\\], independent of the normal PML grading.
  - `exterior_mesh_size_factor=1.0`: Maximum remote-to-central bulk element-size
    ratio \\[dimensionless\\], finite and at least one. Grading starts beyond
    twice the resolution radius; conductor targets remain fixed. Air sizes
    remain bounded by the wavelength. Also controls tangential PML divisions.
  - `interface_refinement_factor=1.0`: Multiplier for the cable's projected
    interface footprint \\[dimensionless\\], finite and at least one. Larger
    values retain more interface refinement in both air and soil. Local wave
    sizes and transition distances retain their frequency/material/Γ dependence.
  - `volume_quadrature=12`: Triangle volume integration points (4, 7, 12 or 13).
  - `physical_volume_quadrature=nothing`: Override triangle integration in
    physical materials with 3, 4, 7, 12 or 13 points. `nothing` inherits
    `volume_quadrature`; triangular PML keeps `volume_quadrature` independently.
  - `pml_element_family=:triangle`: PML cells, either `:triangle` or
    `:quadrangle`. Quadrangles recombine the same native transfinite grid;
    physical materials and voltage-path edges retain their construction.
  - `pml_quadrature=9`: Gauss-Legendre points per PML quadrangle (4, 9 or 16).
    Unused for triangular PML. These are prescribed rules, not error tolerances.
  - `conductor_geometry_tolerance=1e-3`: Requested relative circle-area
    construction tolerance \\[dimensionless\\]; selects at least 96 segments
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
  - `getdp_executable=nothing`: Executable override; otherwise resolve the
    environment override, package artifact, or unsupported-platform `PATH` fallback.
  - `gmsh_verbosity=2`, `getdp_verbosity=2`: Native message levels from 0 through 5.
  - `frequency_workers=2`: Maximum concurrent frequency solver processes.
  - `solver_threads=1`: BLAS and OpenMP threads per solver process.
  - `mumps_ordering=nothing`: Optional native MUMPS `ICNTL(7)` ordering code:
    0 (AMD), 2 (AMF), 3 (Scotch), 4 (PORD), 5 (METIS), 6 (QAMD), or 7
    (automatic). Availability depends on the GetDP/MUMPS build. `nothing`
    leaves the solver default unchanged; no alternative ordering is retried.
  - `petsc_prealloc=nothing`: Optional positive sparse-row allocation estimate
    passed to GetDP's `-petsc_prealloc`. Larger estimates can reduce reallocations
    at the expense of memory. `nothing` retains GetDP's default.
  - `verbosity=(default=0,)`: Julia logging levels from 0 through 2.
  - `output_basis=:pul`: Per-unit-length matrices; `:total` scales by line length.
  - `trace=false`: Retain primitive matrices.
  - `timing=false`: Retain fresh complete-scan measurements in result details.
  - `on_result=nothing`: Optional callback `(problem, index, result)`.
  - `log_file=nothing`: Optional Julia log path.
  - `resume_run_directory=nothing`: Resume a compatible path or `:latest`.

# Returns

- A fixed-key [`ComputationOptions`](@ref) named tuple. Paths and integer controls
  are normalized; `output_basis` and `trace` are lowered to `Val` values.

# Errors

- `ArgumentError`: Unknown keys, invalid controls, or empty paths.
"""
function computation_options(::Type{LineCableModelsFEM}, record::ComputationOptions)::ComputationOptions
    options = record.data
    defaults = (plot_field_maps=false, mesh_policy=:reuse,
        mesh_path=nothing, domain_skin_depths=2.0,
        pml_thickness=nothing, pml_thickness_factor=1.0, pml_layers=128,
        pml_grading=(192/191)*log(1536), pml_resolution=nothing, pml_reflection=1e-10,
        mesh_size_factor=1.0, exterior_mesh_size_factor=1.0,
        interface_refinement_factor=1.0, volume_quadrature=12,
        physical_volume_quadrature=nothing, pml_element_family=:triangle,
        pml_quadrature=9,
        conductor_geometry_tolerance=1e-3, conductor_skin_depth_elements=6.0,
        conductor_mesh_growth=sqrt(1.25), conductor_skin_depths=5.0,
        conductor_thickness_elements=4,
        keep_run_directory=false, getdp_executable=nothing,
        gmsh_verbosity=2, getdp_verbosity=2, frequency_workers=2, solver_threads=1,
        mumps_ordering=nothing, petsc_prealloc=nothing,
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
    radius_factor = normalized.domain_skin_depths
    radius_factor isa Real && !(radius_factor isa Bool) &&
        isfinite(radius_factor) && 0 < radius_factor <= floatmax(Float64) &&
        Float64(radius_factor) > 0 || throw(ArgumentError(
            "domain_skin_depths must be a finite positive Float64-representable number"))
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
    thickness_factor isa Real && !(thickness_factor isa Bool) &&
        isfinite(thickness_factor) && 0 < thickness_factor <= floatmax(Float64) &&
        Float64(thickness_factor) > 0 || throw(ArgumentError(
            "pml_thickness_factor must be a finite positive Float64-representable number"))
    reflection = normalized.pml_reflection
    reflection isa Real && !(reflection isa Bool) && 0 < reflection < 1 &&
        0 < Float64(reflection) < 1 || throw(ArgumentError(
            "pml_reflection must be strictly between zero and one"))
    thickness = normalized.pml_thickness
    thickness = thickness isa Real ? (thickness, thickness, thickness) : thickness
    if thickness !== nothing
        thickness isa Tuple && length(thickness) == 3 &&
            all(t -> t isa Real && !(t isa Bool) && isfinite(t) &&
                0 < t <= floatmax(Float64) && Float64(t) > 0, thickness) ||
            throw(ArgumentError("pml_thickness must be nothing, a positive thickness, or a (side, top, bottom) tuple in metres"))
        thickness = Float64.(thickness)
    end
    resolution = normalized.pml_resolution
    if resolution !== nothing
        resolution isa NamedTuple &&
            isempty(setdiff(keys(resolution),(:interpolation_cells,:coefficient_change))) ||
            throw(ArgumentError("pml_resolution must be a named tuple with interpolation_cells and coefficient_change"))
        resolution = merge((interpolation_cells=72,coefficient_change=0.12),resolution)
        cells, change = resolution.interpolation_cells, resolution.coefficient_change
        cells isa Integer && !(cells isa Bool) && 1 <= cells < typemax(Cint) ||
            throw(ArgumentError("pml_resolution.interpolation_cells must be a positive Cint-representable count"))
        change isa Real && !(change isa Bool) && isfinite(change) &&
            0 < change <= floatmax(Float64) && Float64(change) > 0 ||
            throw(ArgumentError("pml_resolution.coefficient_change must be finite, positive and Float64-representable"))
        all(name -> get(options,name,nothing) === nothing,(:pml_layers,:pml_grading)) ||
            throw(ArgumentError("pml_resolution cannot be combined with explicit pml_layers or pml_grading"))
        resolution = (interpolation_cells=Int(cells),coefficient_change=Float64(change))
    end
    layers = resolution === nothing ? normalized.pml_layers : nothing
    grading = resolution === nothing ? normalized.pml_grading : nothing
    if resolution === nothing
        layers = layers isa Integer ? (layers, layers, layers) : layers
        layers isa Tuple && length(layers) == 3 &&
            all(n -> n isa Integer && !(n isa Bool) && 1 <= n < typemax(Cint), layers) ||
            throw(ArgumentError("pml_layers must be a positive integer or a (side, top, bottom) tuple; each count plus one must fit Gmsh's Cint node count"))
        layers = Int.(layers)
        grading = grading isa Real ? (grading, grading, grading) : grading
        grading isa Tuple && length(grading) == 3 &&
            all(g -> g isa Real && !(g isa Bool) && isfinite(g) &&
                0 <= g <= floatmax(Float64) && (iszero(g) || Float64(g) > 0), grading) ||
            throw(ArgumentError("pml_grading must be a finite nonnegative Float64-representable exponent or a (side, top, bottom) tuple"))
        grading = Float64.(grading)
    end
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
    resume === nothing || resume === :latest || resume isa AbstractString && !isempty(resume) ||
        throw(ArgumentError("resume_run_directory must be nothing, :latest, or a nonempty path string"))
    return ComputationOptions(; standard.data...,
        plot_field_maps=normalized.plot_field_maps,
        mesh_policy=normalized.mesh_policy,
        mesh_path=normalized.mesh_path === nothing ? nothing : String(normalized.mesh_path),
        domain_skin_depths=Float64(radius_factor),
        pml_thickness=thickness, pml_thickness_factor=Float64(thickness_factor),
        pml_layers=layers, pml_grading=grading, pml_resolution=resolution,
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
        mumps_ordering=ordering === nothing ? nothing : Int(ordering),
        petsc_prealloc=prealloc === nothing ? nothing : Int(prealloc),
        log_file=normalized.log_file === nothing ? nothing : String(normalized.log_file),
        resume_run_directory=resume isa AbstractString ? String(resume) : resume)
end
