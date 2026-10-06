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

# One declaration owns the expert names, defaults, types and allowed ranges.
# Each override names exactly one native parameter.
const FEM_NATIVE_OVERRIDES = (
    (names=(:PmlSideThicknessFactor,:PmlTopThicknessFactor,:PmlBottomThicknessFactor), field=:pml_thickness_factor, mesh=true, default=1., kind=:real, range=(0.,Inf), closed=false),
    (names=(:PmlSideLayers,:PmlTopLayers,:PmlBottomLayers), field=:pml_layers, mesh=true, default=16, kind=:integer, range=(1,typemax(Cint)-1), closed=true),
    (names=(:PmlSideGrading,:PmlTopGrading,:PmlBottomGrading), field=:pml_grading, mesh=true, default=log(20), kind=:real, range=(0.,Inf), closed=true),
    (names=(:ExteriorMeshSizeFactor,), field=:exterior_mesh_size_factor, mesh=true, default=1., kind=:real, range=(1.,Inf), closed=true),
    (names=(:InterfaceRefinementFactor,), field=:interface_refinement_factor, mesh=true, default=1., kind=:real, range=(1.,Inf), closed=true),
    (names=(:ConductorGeometryTolerance,), field=:conductor_geometry_tolerance, mesh=true, default=1e-3, kind=:real, range=(0.,Inf), closed=false),
    (names=(:ConductorSkinDepthElements,), field=:conductor_skin_depth_elements, mesh=true, default=6., kind=:real, range=(0.,Inf), closed=false),
    (names=(:ConductorMeshGrowth,), field=:conductor_mesh_growth, mesh=true, default=sqrt(1.25), kind=:real, range=(1.,Inf), closed=true),
    (names=(:ConductorSkinDepths,), field=:conductor_skin_depths, mesh=true, default=5., kind=:real, range=(0.,Inf), closed=false),
    (names=(:ConductorThicknessElements,), field=:conductor_thickness_elements, mesh=true, default=4, kind=:integer, range=(1,typemax(Cint)), closed=true),
    (names=(:VolumeQuadrature,), field=:volume_quadrature, mesh=false, default=12, kind=:integer, range=(4,7,12,13), closed=:choices),
    (names=(:PhysicalVolumeQuadrature,), field=:physical_volume_quadrature, mesh=false, default=3, kind=:integer, range=(3,4,7,12,13), closed=:choices),
    (names=(:PmlQuadrangles,), field=:pml_element_family, mesh=true, default=1, kind=:integer, range=(0,1), closed=:choices),
    (names=(:PmlQuadrature,), field=:pml_quadrature, mesh=false, default=4, kind=:integer, range=(4,9,16), closed=:choices),
    (names=(:MumpsOrdering,), field=:mumps_ordering, mesh=false, default=0, kind=:integer, range=(-1,0,2,3,4,5,6,7), closed=:choices),
    (names=(:PetscPrealloc,), field=:petsc_prealloc, mesh=false, default=0, kind=:integer, range=(0,typemax(Cint)), closed=true),
    (names=(:MumpsErrorAnalysis,), field=:mumps_error_analysis, mesh=false, default=2, kind=:integer, range=(0,1,2), closed=:choices),
    (names=(:MumpsRefinementMax,), field=:mumps_refinement_max, mesh=false, default=2, kind=:integer, range=(0,typemax(Cint)), closed=true),
    (names=(:MumpsBackwardErrorTolerance,), field=:mumps_backward_error_tolerance, mesh=false, default=1e-12, kind=:real, range=(0.,Inf), closed=false),
)
const FEM_MUMPS_FORWARD_ERROR_BUDGET = 0.01

function _fem_override_value(spec, name, value)
    valid_type = spec.kind === :integer ? value isa Integer : value isa Real
    valid_type && !(value isa Bool) && isfinite(value) ||
        throw(ArgumentError("native override $name must be a finite $(spec.kind) value"))
    allowed = spec.closed === :choices ? value in spec.range :
        (spec.closed ? first(spec.range) <= value : first(spec.range) < value) && value <= last(spec.range)
    allowed || throw(ArgumentError("native override $name is outside its allowed range $(spec.range)"))
    spec.kind === :integer && return Int(value)
    converted = Float64(value)
    isfinite(converted) && (iszero(value) || converted > 0) ||
        throw(ArgumentError("native override $name must be Float64-representable"))
    return converted
end

function _fem_expert_controls(overrides)
    overrides isa NamedTuple || throw(ArgumentError("overrides must be a NamedTuple of native parameter names and numeric values"))
    allowed = Tuple(Iterators.flatten(spec.names for spec in FEM_NATIVE_OVERRIDES))
    unknown = setdiff(keys(overrides),allowed)
    isempty(unknown) || throw(ArgumentError("unknown native FEM overrides: $(Tuple(unknown))"))
    values = map(FEM_NATIVE_OVERRIDES) do spec
        map(name -> _fem_override_value(spec,name,get(overrides,name,spec.default)),spec.names)
    end
    return (;
        pml_thickness_factor=map(Float64,values[1]),
        pml_layers=map(Int,values[2]),
        pml_grading=map(Float64,values[3]),
        exterior_mesh_size_factor=Float64(only(values[4])),
        interface_refinement_factor=Float64(only(values[5])),
        conductor_geometry_tolerance=Float64(only(values[6])),
        conductor_skin_depth_elements=Float64(only(values[7])),
        conductor_mesh_growth=Float64(only(values[8])),
        conductor_skin_depths=Float64(only(values[9])),
        conductor_thickness_elements=Int(only(values[10])),
        volume_quadrature=Int(only(values[11])),
        physical_volume_quadrature=Int(only(values[12])),
        pml_element_family=Int(only(values[13])),
        pml_quadrature=Int(only(values[14])),
        mumps_ordering=Int(only(values[15])),
        petsc_prealloc=Int(only(values[16])),
        mumps_error_analysis=Int(only(values[17])),
        mumps_refinement_max=Int(only(values[18])),
        mumps_backward_error_tolerance=Float64(only(values[19])),
    )
end

# The flat record is private lowering to the existing native writers and workers.
# Public keyword names are checked before lowering, so these are not aliases.
"""
$(TYPEDSIGNATURES)

Validate FEM computation options and lower expert overrides to native parameters.
The field formulation and prescribed `Γ` belong to `Formulation`.

# Keywords

- Physics: `mesh_size_factor=1.0` scales every native resolution target, including
  conductor angular and normal sizes, interface and wave sizes, measurement lines,
  and the PML points-per-wavelength bound. Larger values coarsen the mesh; fixed
  minimum circle segments and PML intervals remain fixed. The circle segment count
  set by geometric area tolerance is independent of `mesh_size_factor`. `domain_size_factor=2.0`
  multiplies the earth sizing length `min(1/Re(q₀), 2π/abs(q₀))`, subject to the
  resistive ceiling. `pml_reflection=1e-3` sets the nominal attenuation target.
  Near DC, the added real stretch is limited to 100 domain half-widths and is purely
  real where the air field is quasi-static over that extent; read-only ONELAB values
  report activation.
- Workflow: `mesh_policy=:reuse`, `mesh_path=nothing`,
  `keep_run_directory=false`, `resume_run_directory=nothing`,
  `getdp_executable=nothing`, `output_basis=:pul`, `on_result=nothing`.
  `mesh_path` accepts one `.msh` file or a managed run directory whose meshes
  are matched by fingerprint; unmatched frequencies mesh normally. It cannot
  be combined with `mesh_policy=:remesh`.
  Resume accepts a run directory or `:latest` for the same system.
- Resources: `frequency_workers=clamp(Sys.CPU_THREADS ÷ 4, 1, 8)` and
  `solver_threads=1`. Very large meshes may require fewer workers to bound memory.
- Solver: `linear_solver=:mumps` (also `:gmres`), `gmres_iterations_max=20`,
  `gmres_relative_tolerance=1e-12`, `gmres_absolute_tolerance=0.0`.
- Diagnostics: `plot_field_maps=false`, `verbosity=(default=0,)`,
  `gmsh_verbosity=2`, `getdp_verbosity=2`, `log_file=nothing`, `trace=false`,
  `timing=false`.
- Expert: `overrides=(;)` contains numeric native parameter names and values,
  validated against `FEM_NATIVE_OVERRIDES`. For example,
  `overrides=(PmlSideLayers=48, PmlTopLayers=48, PmlBottomLayers=48,
  PhysicalVolumeQuadrature=3, MumpsOrdering=0)`. Each name controls one native parameter.
  The default physical rule has 3 points, which integrate first-order field
  products exactly for piecewise-constant materials; PML quadrature is separate.
  The MUMPS forward-error warning budget is fixed at 0.01.

# Returns

A validated `ComputationOptions` record for native execution. Removed public
expert keywords and unknown names raise `ArgumentError`; use native overrides.
"""
function computation_options(::Type{LineCableModelsFEM}, record::ComputationOptions)::ComputationOptions
    options=record.data
    defaults=(mesh_size_factor=1.,domain_size_factor=2.,pml_reflection=1e-3,
        mesh_policy=:reuse,mesh_path=nothing,keep_run_directory=false,resume_run_directory=nothing,
        getdp_executable=nothing,frequency_workers=clamp(Sys.CPU_THREADS ÷ 4,1,8),solver_threads=1,
        linear_solver=:mumps,gmres_iterations_max=20,gmres_relative_tolerance=1e-12,gmres_absolute_tolerance=0.,
        plot_field_maps=false,gmsh_verbosity=2,getdp_verbosity=2,log_file=nothing,overrides=(;))
    standard_keys=(:verbosity,:output_basis,:trace,:on_result,:timing)
    for spec in FEM_NATIVE_OVERRIDES
        haskey(options,spec.field) && throw(ArgumentError(
            "$(spec.field) was removed; use overrides=($(first(spec.names))=...,)"))
    end
    haskey(options,:mumps_forward_error_tolerance) && throw(ArgumentError(
        "mumps_forward_error_tolerance was removed; the forward-error warning budget is fixed at $FEM_MUMPS_FORWARD_ERROR_BUDGET"))
    unknown=setdiff(keys(options),(keys(defaults)...,standard_keys...))
    isempty(unknown) || throw(ArgumentError("unknown LineCableModelsFEM computation options: $(Tuple(unknown))"))
    standard=computation_options(LineCableModelsCoaxial,
        ComputationOptions(; verbosity=get(options,:verbosity,(default=0,)),
            output_basis=get(options,:output_basis,:pul),trace=get(options,:trace,false),
            on_result=get(options,:on_result,nothing),timing=get(options,:timing,false)))
    controls=merge(defaults,options)
    expert=_fem_expert_controls(controls.overrides)
    for name in (:mesh_size_factor,:domain_size_factor,:gmres_relative_tolerance,:gmres_absolute_tolerance)
        value=getproperty(controls,name);allow_zero=name===:gmres_absolute_tolerance
        value isa Real && !(value isa Bool) && isfinite(value) &&
            0 <= value <= floatmax(Float64) && (Float64(value)>0 || allow_zero && iszero(value)) ||
            throw(ArgumentError("$name must be finite, $(allow_zero ? "nonnegative" : "positive") and Float64-representable"))
    end
    reflection=controls.pml_reflection
    reflection isa Real && !(reflection isa Bool) && 0<Float64(reflection)<1 ||
        throw(ArgumentError("pml_reflection must be strictly between zero and one"))
    for name in (:plot_field_maps,:keep_run_directory)
        getproperty(controls,name) isa Bool || throw(ArgumentError("$name must be Bool"))
    end
    controls.mesh_policy in (:reuse,:remesh) || throw(ArgumentError("mesh_policy must be :reuse or :remesh"))
    controls.mesh_path !== nothing && controls.mesh_policy===:remesh &&
        throw(ArgumentError("mesh_path cannot be combined with mesh_policy=:remesh"))
    for name in (:mesh_path,:getdp_executable,:log_file)
        value=getproperty(controls,name)
        value===nothing || value isa AbstractString && !isempty(value) ||
            throw(ArgumentError("$name must be a nonempty path string or nothing"))
    end
    if controls.mesh_path !== nothing && (isdir(controls.mesh_path) ||
            !isfile(controls.mesh_path) && !endswith(lowercase(controls.mesh_path),".msh"))
        directory=joinpath(controls.mesh_path,"mesh")
        isdir(directory) && any(name->occursin(r"^frequency_\d{4,}\.json$",name),readdir(directory)) ||
            throw(ArgumentError("mesh_path run directory has no mesh sidecars: $(controls.mesh_path)"))
    end
    resume=controls.resume_run_directory
    resume===nothing || resume===:latest || resume isa AbstractString && !isempty(resume) ||
        throw(ArgumentError("resume_run_directory must be nothing, :latest, or a nonempty path string"))
    for name in (:gmsh_verbosity,:getdp_verbosity)
        value=getproperty(controls,name)
        value isa Integer && !(value isa Bool) && value in 0:5 || throw(ArgumentError("$name must be an integer from 0 through 5"))
    end
    for name in (:frequency_workers,:solver_threads,:gmres_iterations_max)
        value=getproperty(controls,name);ceiling=name===:gmres_iterations_max ? typemax(Cint) : typemax(Int)
        value isa Integer && !(value isa Bool) && 1<=value<=ceiling || throw(ArgumentError("$name must be a positive representable integer"))
    end
    controls.linear_solver in (:mumps,:gmres) || throw(ArgumentError("linear_solver must be :mumps or :gmres"))
    return ComputationOptions(;
        verbosity=standard.data.verbosity,output_basis=standard.data.output_basis,
        trace=standard.data.trace,on_result=standard.data.on_result,timing=standard.data.timing,
        mesh_size_factor=Float64(controls.mesh_size_factor),
        domain_size_factor=Float64(controls.domain_size_factor),
        pml_reflection=Float64(controls.pml_reflection),
        mesh_policy=controls.mesh_policy,
        mesh_path=controls.mesh_path isa AbstractString ? String(controls.mesh_path) : controls.mesh_path,
        keep_run_directory=controls.keep_run_directory,
        resume_run_directory=controls.resume_run_directory isa AbstractString ? String(controls.resume_run_directory) : controls.resume_run_directory,
        getdp_executable=controls.getdp_executable isa AbstractString ? String(controls.getdp_executable) : controls.getdp_executable,
        frequency_workers=Int(controls.frequency_workers),
        solver_threads=Int(controls.solver_threads),
        linear_solver=controls.linear_solver,
        gmres_iterations_max=Int(controls.gmres_iterations_max),
        gmres_relative_tolerance=Float64(controls.gmres_relative_tolerance),
        gmres_absolute_tolerance=Float64(controls.gmres_absolute_tolerance),
        plot_field_maps=controls.plot_field_maps,
        gmsh_verbosity=Int(controls.gmsh_verbosity),
        getdp_verbosity=Int(controls.getdp_verbosity),
        log_file=controls.log_file isa AbstractString ? String(controls.log_file) : controls.log_file,
        pml_thickness_factor=expert.pml_thickness_factor,
        pml_layers=expert.pml_layers,
        pml_grading=expert.pml_grading,
        exterior_mesh_size_factor=expert.exterior_mesh_size_factor,
        interface_refinement_factor=expert.interface_refinement_factor,
        conductor_geometry_tolerance=expert.conductor_geometry_tolerance,
        conductor_skin_depth_elements=expert.conductor_skin_depth_elements,
        conductor_mesh_growth=expert.conductor_mesh_growth,
        conductor_skin_depths=expert.conductor_skin_depths,
        conductor_thickness_elements=expert.conductor_thickness_elements,
        volume_quadrature=expert.volume_quadrature,
        physical_volume_quadrature=expert.physical_volume_quadrature,
        pml_element_family=expert.pml_element_family,
        pml_quadrature=expert.pml_quadrature,
        mumps_ordering=expert.mumps_ordering,
        petsc_prealloc=expert.petsc_prealloc,
        mumps_error_analysis=expert.mumps_error_analysis,
        mumps_refinement_max=expert.mumps_refinement_max,
        mumps_backward_error_tolerance=expert.mumps_backward_error_tolerance,
    )
end
