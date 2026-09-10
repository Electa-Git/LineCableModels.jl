# Engine-owned formulation hierarchy.
"""
$(TYPEDEF)

Select the LineCableModels backend for concentric coaxial cable assemblies.

Nonconcentric cable parts must reach this backend through an equivalent
concentric representation supplied by DataModel.
"""
struct LineCableModelsCoaxial end

"""
$(TYPEDEF)

Supertype for formulation values that select a complete numerical backend.
"""
abstract type AbstractFormulationBackend <: AbstractFormulation end

"""
$(TYPEDEF)

Supertype for backend-owned formulation option records.
"""
abstract type AbstractFormulationOptions end

"""
$(TYPEDEF)

Configure execution owned by the Gmsh/GetDP finite-element backend.

The problem supplies geometry, reference materials, frequencies, and prescribed
temperature. The formulation selects constitutive laws and matrix reductions;
this record contains FEM execution controls.

$(TYPEDFIELDS)
"""
struct LineCableModelsFEMOptions <: AbstractFormulationOptions
    "Open the optional Gmsh graphical interface."
    ui::Bool
    "Generate field-map files for every frequency and basis terminal."
    plot_field_maps::Bool
    "Mesh selection policy, either `:reuse` or `:remesh`."
    mesh_policy::Symbol
    "Optional existing Gmsh mesh path."
    mesh_path::Union{Nothing, String}
    "Retain a successful run directory."
    keep_run_directory::Bool
    "Optional GetDP executable override."
    getdp_executable::Union{Nothing, String}
    "Gmsh message verbosity from 0 through 5."
    gmsh_verbosity::Int
    "GetDP message verbosity from 0 through 5."
    getdp_verbosity::Int
    "Maximum number of independent frequency solver processes."
    frequency_workers::Int
    "BLAS and OpenMP thread limit for each GetDP process."
    solver_threads::Int
end

"""
$(TYPEDSIGNATURES)

Construct validated finite-element execution options.

# Keywords

- `ui=false`: Open the Gmsh graphical interface.
- `plot_field_maps=false`: Emit spatial field maps for every solve.
- `mesh_policy=:reuse`: Reuse a compatible mesh or generate one. Use
  `:remesh` to regenerate it unconditionally.
- `mesh_path=nothing`: Optional existing `.msh` file.
- `keep_run_directory=false`: Retain successful run artifacts.
- `getdp_executable=nothing`: Let the Gmsh extension resolve GetDP from its
  environment override, package-owned artifact, or unsupported-platform
  `PATH` fallback. A supplied path overrides all three.
- `gmsh_verbosity=2`: Gmsh message verbosity from 0 through 5.
- `getdp_verbosity=2`: GetDP message verbosity from 0 through 5.
- `frequency_workers=2`: Maximum concurrent GetDP frequency processes. Use
  `1` for sequential frequency execution; terminal factors are still reused.
- `solver_threads=1`: BLAS and OpenMP threads per GetDP process. Set the
  worker and thread counts together to fit available memory and CPU resources.

# Returns

- A validated [`LineCableModelsFEMOptions`](@ref) value.

# Errors

- `ArgumentError`: An option has an unsupported value or an empty path.
"""
function LineCableModelsFEMOptions(;
        ui::Bool = false,
        plot_field_maps::Bool = false,
        mesh_policy::Symbol = :reuse,
        mesh_path::Union{Nothing, AbstractString} = nothing,
        keep_run_directory::Bool = false,
        getdp_executable::Union{Nothing, AbstractString} = nothing,
        gmsh_verbosity::Integer = 2,
        getdp_verbosity::Integer = 2,
        frequency_workers::Integer = 2,
        solver_threads::Integer = 1
)
    mesh_policy in (:reuse, :remesh) || throw(ArgumentError(
        "mesh_policy must be :reuse or :remesh; got $(repr(mesh_policy))",
    ))
    normalized_mesh_path = mesh_path === nothing ? nothing : String(mesh_path)
    normalized_getdp = getdp_executable === nothing ? nothing :
                       String(getdp_executable)
    normalized_mesh_path === "" && throw(ArgumentError("mesh_path cannot be empty"))
    normalized_getdp === "" && throw(ArgumentError(
        "getdp_executable cannot be empty",
    ))
    gmsh_verbosity in 0:5 || throw(ArgumentError(
        "gmsh_verbosity must be an integer from 0 through 5",
    ))
    getdp_verbosity in 0:5 || throw(ArgumentError(
        "getdp_verbosity must be an integer from 0 through 5",
    ))
    for (name, value) in ((:frequency_workers, frequency_workers),
        (:solver_threads, solver_threads))
        !(value isa Bool) && 1 <= value <= typemax(Int) ||
            throw(ArgumentError("$name must be a positive integer"))
    end
    return LineCableModelsFEMOptions(
        ui,
        plot_field_maps,
        mesh_policy,
        normalized_mesh_path,
        keep_run_directory,
        normalized_getdp,
        Int(gmsh_verbosity),
        Int(getdp_verbosity),
        Int(frequency_workers),
        Int(solver_threads)
    )
end

"""
$(TYPEDEF)

Select the Julia-native Gmsh/GetDP quasi-TEM finite-element backend.

`options` stores the shared LineCableModels reduction and temperature policy.
`execution` stores only Gmsh/GetDP execution controls.

$(TYPEDFIELDS)
"""
struct LineCableModelsFEM{M <: NamedTuple, O <: NamedTuple, D <: NamedTuple} <:
       AbstractFormulationBackend
    "Shared scientific formula selections, independent of FEM execution controls."
    methods::M
    "Shared line-parameter formulation options."
    options::O
    "Requested formula definitions retained for provenance."
    definitions::D
    "Finite-element execution options."
    execution::LineCableModelsFEMOptions
end

"""
$(TYPEDEF)

Report a finite-element adaptation, mesh, solve, or result-contract failure.

$(TYPEDFIELDS)
"""
struct LineCableModelsFEMError <: Exception
    "Failure category."
    category::Symbol
    "Stable identifier of the object that owns the failure."
    object_id::String
    "Field or derived datum that failed validation."
    field::Symbol
    "Human-readable failure description."
    message::String
    "Retained run directory, or `nothing` before run creation."
    run_directory::Union{Nothing, String}
end

function LineCableModelsFEMError(
        category::Symbol,
        object_id,
        field::Symbol,
        message::AbstractString;
        run_directory::Union{Nothing, AbstractString} = nothing
)
    path = run_directory === nothing ? nothing : String(run_directory)
    return LineCableModelsFEMError(
        category, String(object_id), field, String(message), path
    )
end

function Base.showerror(io::IO, error::LineCableModelsFEMError)
    print(
        io,
        "LineCableModelsFEMError(",
        error.category,
        ", object=",
        repr(error.object_id),
        ", field=:",
        error.field,
        "): ",
        error.message
    )
    error.run_directory === nothing || print(
        io, "; retained run directory: ", error.run_directory
    )
end

"""
$(TYPEDEF)

Supertype for Engine impedance formulations.
"""
abstract type AbstractImpedanceFormulation <: AbstractFormulation end
abstract type InternalImpedanceFormulation <: AbstractImpedanceFormulation end
abstract type PipeImpedanceFormulation <: AbstractImpedanceFormulation end
abstract type InsulationImpedanceFormulation <: AbstractImpedanceFormulation end
abstract type EarthImpedanceFormulation <: AbstractImpedanceFormulation end

abstract type AbstractAdmittanceFormulation <: AbstractFormulation end
abstract type InsulationAdmittanceFormulation <: AbstractAdmittanceFormulation end
abstract type SemiconAdmittanceFormulation <: AbstractAdmittanceFormulation end
abstract type EarthAdmittanceFormulation <: AbstractAdmittanceFormulation end

"Return whether an earth formulation consumes homogeneous or stratified media."
function media end

"Declare source-owned physical hook defaults and admitted overrides for an equation binding."
function hooks end

"Resolve scalar or homogeneous three-field selections through their formula owner."
function Formulation(::Type{F}, selected) where {F <: Union{
        EarthImpedanceFormulation, EarthAdmittanceFormulation}}
    return F(selected)
end

function Formulation(::Type{F}, selected::NamedTuple) where {F <: Union{
        EarthImpedanceFormulation, EarthAdmittanceFormulation}}
    names = (:air, :earth, :mixed)
    length(selected) == 3 && all(in(keys(selected)), names) || throw(ArgumentError(
        "homogeneous earth selections require exactly air, earth and mixed"))
    return map(F, NamedTuple{names}(selected))
end

"Resolve the selected formula for exact source and target layer indices."
Formulation(selected::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
    ::Val{S}, ::Val{T}) where {S, T} = selected

Formulation(selected::NamedTuple{(:air, :earth, :mixed)}, ::Val{1}, ::Val{1}) = selected.air
Formulation(selected::NamedTuple{(:air, :earth, :mixed)}, ::Val{2}, ::Val{2}) = selected.earth
Formulation(selected::NamedTuple{(:air, :earth, :mixed)}, ::Val{1}, ::Val{2}) = selected.mixed
Formulation(selected::NamedTuple{(:air, :earth, :mixed)}, ::Val{2}, ::Val{1}) = selected.mixed

function Formulation(::NamedTuple{(:air, :earth, :mixed)}, ::Val{S}, ::Val{T}) where {S, T}
    throw(ArgumentError(
        "homogeneous selection is not defined for source in layer $S and target in layer $T"))
end

function validate(selected::NamedTuple{(:air, :earth, :mixed)}, earth::EarthModel)
    validate(earth)
    !earth.vertical_layers && length(earth.layers) == 2 &&
        all(layer -> isinf(layer.thickness), earth.layers) || throw(ArgumentError(
        "air/earth/mixed selections require physical air and one homogeneous soil half-space; use a scalar formulation for a layered model"))
    return selected
end

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        pair::EarthPair)
    return only(validate(formula, (pair,)))
end

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        pairs::Union{Tuple, AbstractVector{<:EarthPair}})
    equations = map(pairs) do pair
        validate(pair)
        equation = validate(FormulaMethod(formula, pair))
        validate(pair, equation)
        equation
    end
    identities = unique(equations)
    bindings = map(identities) do equation
        declared = hooks(equation)
        all(in(declared.configurable), keys(formula.hooks)) || throw(ArgumentError(
            "an explicit physical hook is unused by $equation"))
        selected_hooks = merge(declared.defaults, formula.hooks)
        defaults = selected_hooks.contribution === nothing ? computation_options(equation) :
                   computation_options(equation, selected_hooks.contribution)
        (equation = equation, kind = typeof(first(equation.arguments)).parameters[1],
            hooks = selected_hooks, defaults = defaults)
    end
    admitted = union((keys(binding.defaults) for binding in bindings)...)
    unknown = setdiff(keys(formula.options), admitted)
    isempty(unknown) || throw(ArgumentError(
        "unused numerical sections $(Tuple(unknown)) for required cases of :$(formula_id(formula))"))
    resolved = map(bindings) do binding
        names = Tuple(intersect(keys(formula.options), keys(binding.defaults)))
        options = computation_options(binding.equation, binding.defaults, formula.options[names])
        (equation = binding.equation, kind = binding.kind,
            hooks = binding.hooks, options = options)
    end
    return map(equation -> resolved[findfirst(==(equation), identities)], equations)
end

# Equation-specific geometric restrictions extend the existing validation protocol.
validate(pair::EarthPair, ::FormulaMethod) = pair

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation}, count::Integer)
    count in formula.assumptions.layers || throw(DimensionMismatch(
        "formula :$(formula_id(formula)) requires $(formula.assumptions.layers) media including air; received $count"))
    return formula
end

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation}, earth::EarthModel)
    validate(earth)
    earth.vertical_layers &&
        throw(ArgumentError("earth-return equations require horizontal interfaces or an explicit EquivalentHomogeneous reduction"))
    validate(formula, length(earth.layers))
    return formula
end

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        rho::AbstractVector, epsilon::AbstractVector, mu::AbstractVector, thickness)
    length(rho) == length(epsilon) == length(mu) ||
        throw(DimensionMismatch("material vectors must align"))
    validate(formula, length(rho))
    all(x -> x > 0 && !isnan(x), rho) ||
        throw(DomainError(rho, "resistivities must be positive, including infinite air resistivity"))
    all(x -> isfinite(x) && !iszero(x), epsilon) && all(x -> isfinite(x) && x > 0, mu) ||
        throw(DomainError((epsilon, mu),
            "permittivities must be nonzero and finite; permeabilities positive and finite"))
    epsilon[1] > 0 || throw(DomainError(epsilon[1], "air permittivity must be positive"))
    restriction = formula.assumptions.permittivity
    restriction in (:positive, :nonzero) ||
        throw(ArgumentError("unknown source permittivity restriction"))
    restriction === :positive && !all(>(0), epsilon) &&
        throw(DomainError(epsilon,
            "formula :$(formula_id(formula)) requires positive permittivity; the earth-material data type permits artificial negative values"))
    if thickness === nothing
        media(formula) === Val(:stratified) && length(rho) > 2 &&
            throw(DimensionMismatch("stratified equations require aligned physical layer thicknesses"))
    else
        length(thickness) == length(rho) ||
            throw(DimensionMismatch("layer thicknesses must align with materials"))
        isinf(first(thickness)) && isinf(last(thickness)) &&
        all(x -> isfinite(x) && x > 0, @view(thickness[2:(end - 1)])) ||
            throw(DomainError(thickness,
                "air and bottom half-spaces must be infinite; internal layers positive and finite"))
        media(formula) === Val(:homogeneous) && length(thickness) != 2 &&
            throw(DimensionMismatch("homogeneous equations have no internal soil interfaces"))
    end
    return formula
end

"Route an explicit external formulation tag to its `Val` dispatch method."
Formulation(backend::Symbol; kwargs...) = Formulation(Val(backend); kwargs...)

function _fem_execution_options(options::LineCableModelsFEMOptions)
    return options
end

function _fem_execution_options(options::NamedTuple)
    return LineCableModelsFEMOptions(; options...)
end

function _fem_formulation(
        insulation_admittance, semicon_admittance, earth_properties, temperature_dependence,
        options::NamedTuple,
        fem_options::Union{NamedTuple, LineCableModelsFEMOptions}
)
    methods = (
        insulation_admittance = InsulationAdmittance.Formula(insulation_admittance),
        semicon_admittance = SemiconAdmittance.Formula(semicon_admittance),
        earth_properties = earth_properties === nothing ? nothing :
                           Earth.FrequencyDependent.Formula(earth_properties),
        temperature_dependence = temperature_dependence === nothing ? nothing :
                                 TemperatureDependent.Formula(temperature_dependence),
    )
    definitions = (; insulation_admittance, semicon_admittance, earth_properties,
        temperature_dependence)
    return LineCableModelsFEM(methods, formulation_options(LineCableModelsFEM, options),
        definitions, _fem_execution_options(fem_options))
end

"""
$(TYPEDSIGNATURES)

Construct the Gmsh/GetDP finite-element formulation. FEM owns its field equations
and selects four material laws. Each law, `options`, and `fem_options` accepts a
scalar or an explicit `Grid`/`Gridspace`; varying inputs return a
`Gridspace{LineCableModelsFEM}`.

# Keywords

- `insulation_admittance`: Insulation admittivity law; `:default` is lossless.
- `semicon_admittance`: Semicon admittivity law; `:default` is lossless.
- `earth_properties`: Soil frequency-dependent constitutive law; `:default` and
  `nothing` preserve the declared static soil. Equivalent-earth reductions are
  unsupported. Air uses its declared static properties.
- `temperature_dependence`: Cable-material resistivity law; `:default` selects
  the linear law and `nothing` retains reference resistivity. Operating
  temperature belongs to `LineParametersProblem`.
- `options=(;)`: Shared bundle, Kron, and ideal-transposition reductions.
- `fem_options=(;)`: A `LineCableModelsFEMOptions` value or equivalent named tuple.
- `combine=:product`: Product or zip composition among varying inputs.

Analytical impedance/admittance kernel keywords are rejected. Supported enclosure
geometry is represented directly in the FEM domain.
"""
function Formulation(
        ::Val{:LineCableModelsFEM};
        insulation_admittance = formula(:default),
        semicon_admittance = formula(:default),
        earth_properties = formula(:default),
        temperature_dependence = formula(:default),
        options = (;),
        fem_options = (;),
        combine::Symbol = :product
)
    return parameterize(
        LineCableModelsFEM,
        _fem_formulation,
        (insulation_admittance, semicon_admittance, earth_properties,
            temperature_dependence, options, fem_options);
        combine
    )
end

function LineCableModelsFEM(; kwargs...)
    return Formulation(Val(:LineCableModelsFEM); kwargs...)
end

function validate(binding::FormulaMethod, reduction::EquivalentHomogeneous.AbstractRule)
    throw(ArgumentError("$binding does not admit equivalent-earth reduction :$(formula_id(reduction))"))
end

"""Expose FEM constitutive/admittance selections, reductions and numerical execution controls."""
function Base.NamedTuple(value::LineCableModelsFEM)
    record = function (selected)
        selected === nothing && return nothing
        selected isa Symbol && return NamedTuple(formula(selected))
        selected isa NamedTuple && return map(record,selected)
        return NamedTuple(selected)
    end
    execution=NamedTuple{fieldnames(typeof(value.execution))}(Tuple(getfield(value.execution,key)
        for key in fieldnames(typeof(value.execution))))
    Record=NamedTuple{(:backend,:requested,:methods,:options,:execution),
        Tuple{Symbol,NamedTuple,NamedTuple,NamedTuple,NamedTuple}}
    return Record((:fem,map(record,value.definitions),map(record,value.methods),value.options,execution))
end
