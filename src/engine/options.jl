function formulation_options(
        ::Type{LineParametersFormulation},
        options::NamedTuple
)::FormulationOptions
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
    return (
        reduce_bundle = normalized.reduce_bundle,
        kron_reduction = normalized.kron_reduction,
        ideal_transposition = normalized.ideal_transposition
    )
end

function computation_options(
        ::Type{LineCableModelsCoaxial},
        options::NamedTuple
)::ComputationOptions
    allowed = (:verbosity, :output_basis, :trace, :on_result)
    unknown = filter(key -> key ∉ allowed, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown LineCableModelsCoaxial computation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge(
        (verbosity = (default = 0,), output_basis = :pul,
            trace = false, on_result = nothing),
        options
    )
    verbosity_values = normalized.verbosity
    verbosity_values isa NamedTuple || throw(ArgumentError(
        "verbosity must be a named tuple",
    ))
    haskey(verbosity_values, :default) || throw(ArgumentError(
        "verbosity must define a default level",
    ))
    all(value -> value isa Integer && value in 0:2, values(verbosity_values)) ||
        throw(ArgumentError("verbosity levels must be integers from 0 to 2"))
    basis_value = normalized.output_basis
    basis_value in (:pul, :total) || throw(ArgumentError(
        "output_basis must be :pul or :total; got $(repr(basis_value))",
    ))
    normalized.trace isa Bool || throw(ArgumentError("trace must be Bool"))
    levels = NamedTuple{keys(verbosity_values)}(Int.(values(verbosity_values)))
    return (
        verbosity = levels,
        output_basis = Val(basis_value),
        trace = Val(normalized.trace),
        on_result = normalized.on_result
    )
end

function verbosity(options::NamedTuple, key::Symbol)
    haskey(options, :verbosity) || throw(ArgumentError(
        "computation options do not define verbosity",
    ))
    return get(options.verbosity, key, options.verbosity.default)
end

function formulation_options(::Type{LineCableModelsFEM}, options::NamedTuple)::FormulationOptions
    physics = get(options, :physics, :quasi_tem)
    physics isa Union{Symbol, AbstractString} || throw(ArgumentError(
        "physics must be :quasi_tem or :quasi_fw"))
    physics = Symbol(replace(String(physics), '_' => '-'))
    physics in (Symbol("quasi-tem"), Symbol("quasi-fw")) || throw(ArgumentError(
        "physics must be :quasi_tem or :quasi_fw (also accepted as hyphenated strings or Symbols)"))
    reductions = formulation_options(LineParametersFormulation,
        Base.structdiff(options, (; physics)))
    return (; reductions..., physics)
end

"""
$(TYPEDSIGNATURES)

Validate Gmsh/GetDP execution controls supplied to `compute(...; options)`.
The field model is selected separately by `formulation_options(LineCableModelsFEM, ...)`.

# Arguments

- `owner`: The `LineCableModelsFEM` computation owner.
- `options`: Caller-supplied named tuple. Supported keys and defaults are:
  - `ui=false`: Open the Gmsh graphical interface.
  - `plot_field_maps=false`: Emit spatial field maps for every solve.
  - `mesh_policy=:reuse`: Reuse compatible meshes; `:remesh` regenerates them.
  - `mesh_path=nothing`: Optional existing `.msh` path.
  - `keep_run_directory=false`: Retain successful run artifacts.
  - `getdp_executable=nothing`: Executable override; otherwise resolve the
    environment override, package artifact, or unsupported-platform `PATH` fallback.
  - `gmsh_verbosity=2`, `getdp_verbosity=2`: Native message levels from 0 through 5.
  - `frequency_workers=2`: Maximum concurrent frequency solver processes.
  - `solver_threads=1`: BLAS and OpenMP threads per solver process.
  - `verbosity=(default=0,)`: Julia logging levels from 0 through 2.
  - `output_basis=:pul`: Per-unit-length matrices; `:total` scales by line length.
  - `trace=false`: Retain primitive matrices.
  - `on_result=nothing`: Optional callback `(problem, index, result)`.
  - `log_file=nothing`: Optional Julia log path.
  - `resume_run_directory=nothing`: Resume a compatible path or `:latest`.

# Returns

- A fixed-key [`ComputationOptions`](@ref) named tuple. Paths and integer controls
  are normalized; `output_basis` and `trace` are lowered to `Val` values.

# Errors

- `ArgumentError`: Unknown keys, invalid controls, or empty paths.
"""
function computation_options(::Type{LineCableModelsFEM}, options::NamedTuple)::ComputationOptions
    defaults = (ui=false, plot_field_maps=false, mesh_policy=:reuse,
        mesh_path=nothing, keep_run_directory=false, getdp_executable=nothing,
        gmsh_verbosity=2, getdp_verbosity=2, frequency_workers=2, solver_threads=1,
        log_file=nothing, resume_run_directory=nothing)
    standard_keys = (:verbosity, :output_basis, :trace, :on_result)
    unknown = setdiff(keys(options), (keys(defaults)..., standard_keys...))
    isempty(unknown) || throw(ArgumentError(
        "unknown LineCableModelsFEM computation options: $(Tuple(unknown))"))
    standard = computation_options(LineCableModelsCoaxial,
        Base.structdiff(options, defaults))
    normalized = merge(defaults, Base.structdiff(options, standard))
    for name in (:ui, :plot_field_maps, :keep_run_directory)
        getproperty(normalized, name) isa Bool || throw(ArgumentError("$name must be Bool"))
    end
    normalized.mesh_policy in (:reuse, :remesh) || throw(ArgumentError(
        "mesh_policy must be :reuse or :remesh"))
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
    for name in (:frequency_workers, :solver_threads)
        value = getproperty(normalized, name)
        value isa Integer && !(value isa Bool) && 1 <= value <= typemax(Int) ||
            throw(ArgumentError("$name must be a positive integer"))
    end
    resume = normalized.resume_run_directory
    resume === nothing || resume === :latest || resume isa AbstractString && !isempty(resume) ||
        throw(ArgumentError("resume_run_directory must be nothing, :latest, or a nonempty path string"))
    return (; standard...,
        ui=normalized.ui, plot_field_maps=normalized.plot_field_maps,
        mesh_policy=normalized.mesh_policy,
        mesh_path=normalized.mesh_path === nothing ? nothing : String(normalized.mesh_path),
        keep_run_directory=normalized.keep_run_directory,
        getdp_executable=normalized.getdp_executable === nothing ? nothing : String(normalized.getdp_executable),
        gmsh_verbosity=Int(normalized.gmsh_verbosity), getdp_verbosity=Int(normalized.getdp_verbosity),
        frequency_workers=Int(normalized.frequency_workers), solver_threads=Int(normalized.solver_threads),
        log_file=normalized.log_file === nothing ? nothing : String(normalized.log_file),
        resume_run_directory=resume isa AbstractString ? String(resume) : resume)
end
