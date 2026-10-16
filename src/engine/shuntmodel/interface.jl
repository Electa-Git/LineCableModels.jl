"""
$(TYPEDEF)

Select a local shunt geometry approximation. Material admittivity remains owned
by the insulation and semiconductor constitutive selections.

Each formula implements `dielectric_groups(selected, design, blueprint, formulas,
design_index)`, which returns the dielectric groups that the formula solves in one design and
the records of every concentric assembly. The `:boundary` formula solves each group with
`internal_capacitance`.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions} <: ShuntModelFormulation
    "Model parameters. The built-in selections accept none."
    parameters::P
    "Boundary discretization and quadrature controls."
    options::O
end

"""
$(TYPEDEF)

Report a recognized failure of the local geometric boundary approximation. Numerical
failure is distinct from invalid physical input and must not trigger Monte
Carlo rejection and resampling.

$(TYPEDFIELDS)
"""
struct BoundarySolveError <: Exception
    "Failure stage or unsupported model assumption."
    category::Symbol
    "Physical location and available numerical diagnostics."
    context::NamedTuple
    "Explanation of the failure."
    message::String
end

function Base.showerror(io::IO, error::BoundarySolveError)
    print(io, "BoundarySolveError(", error.category, "): ", error.message)
    isempty(error.context) || print(io, "; ", error.context)
end
Base.show(io::IO, error::BoundarySolveError) = showerror(io, error)
function Base.summary(io::IO, error::BoundarySolveError)
    print(io, "BoundarySolveError(", error.category, ")")
end
function Base.show(io::IO, ::MIME"text/plain", error::BoundarySolveError)
    get(io, :compact, false) && return summary(io, error)
    TextDisplay.tree(io,
        sprint(summary, error),
        (
            (label = error.message, noun = "details"),
            (label = sprint(show, error.context; context = :limit => true),
                noun = "details")
        );
        noun = "details")
end

"Reference boundary discretization. Accuracy depends on geometry and the requested terminal quantity."
const DEFAULT_RESOLUTION = (wire = 64, order = 32, quadrature = 64, modes = 512)
"Default controls for dimensionless logarithmic-moment integration."
const DEFAULT_INTEGRATION = (rtol = 1e-8, atol = 1e-10, maxevals = 100_000)

"""
Return the dielectric groups that a shunt model solves in one cable design and the records of
its concentric assemblies.
"""
function dielectric_groups end

"Return the available local shunt model identifiers."
formulas(::Type{<:Formula}) = (:default, :equivalent, :boundary)

"""
$(TYPEDSIGNATURES)

Construct a cable-local shunt model. `:default` and `:equivalent` use the equivalent annular
layer of each dielectric interval. `:boundary` resolves supported lossless wire and
tape arrangements before frequency evaluation, or reports an unsupported geometry.

# Keywords

- `parameters=(;)`: the model rejects parameters. Numerical failures propagate.
- `options=(;)`: geometric boundary `resolution=(wire=64, order=32, quadrature=256,
  modes=1024)` and `integration=(rtol=1e-8, atol=1e-10, maxevals=100_000)`. The
  equivalent model does not accept numerical controls.

# Returns

- A concrete shunt model selection.
"""
function Formula{ID}(; parameters::NamedTuple = (;),
        options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    options = FormulationOptions(options)
    isempty(parameters) || throw(ArgumentError("shunt model :$ID accepts no parameters"))
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    normalized = formulation_options(Expression(selected, dielectric_groups), options)
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
end

Formula{:default}(; kwargs...) = Formula{:equivalent}(; kwargs...)

formulation_options(::Expression{<:Formula{:equivalent}}) = FormulationOptions()

function formulation_options(::Expression{<:Formula{:boundary}})
    FormulationOptions((resolution = DEFAULT_RESOLUTION, integration = DEFAULT_INTEGRATION))
end

function formulation_options(::Expression{<:Formula{:boundary}}, ::Val{:resolution},
        defaults::NamedTuple, supplied::NamedTuple)
    isempty(setdiff(keys(supplied), keys(defaults))) ||
        throw(ArgumentError("unknown boundary resolution controls"))
    resolution = merge(defaults, supplied)
    all(x -> x isa Integer && !(x isa Bool) && x > 0, values(resolution)) ||
        throw(ArgumentError("boundary resolution controls must be positive integers"))
    resolution.quadrature >= resolution.order+1 || throw(ArgumentError(
        "boundary quadrature must contain at least order+1 nodes"))
    return map(Int, resolution)
end

function formulation_options(::Expression{<:Formula{:boundary}}, ::Val{:integration},
        defaults::NamedTuple, supplied::NamedTuple)
    isempty(setdiff(keys(supplied), keys(defaults))) ||
        throw(ArgumentError("unknown boundary integration controls"))
    integration = merge(defaults, supplied)
    all(x -> x isa Real && isfinite(x) && x >= 0, (integration.rtol, integration.atol)) &&
    max(integration.rtol, integration.atol) > 0 || throw(ArgumentError(
        "boundary integration requires nonnegative finite tolerances, at least one positive"))
    integration.maxevals isa Integer && !(integration.maxevals isa Bool) &&
    integration.maxevals > 0 ||
        throw(ArgumentError("boundary maxevals must be a positive integer"))
    return (rtol = Float64(integration.rtol), atol = Float64(integration.atol),
        maxevals = Int(integration.maxevals))
end

"""Describe the equivalent annular layer local shunt approximation."""
function description(::Type{<:Formula{:equivalent}}; compact::Bool = false)
    compact ? "equivalent annular layer" :
    "Equivalent annular layer shunt geometry"
end
"""Describe the default equivalent annular layer local shunt approximation."""
function description(::Type{<:Formula{:default}}; compact::Bool = false)
    description(Formula{:equivalent}; compact)
end
"""Describe the lossless wire and tape geometric boundary approximation."""
function description(::Type{<:Formula{:boundary}}; compact::Bool = false)
    compact ? "boundary" :
    "Lossless wire/tape boundary shunt geometry"
end

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(::Val{ID}; kwargs...) where {ID} = Formula{ID}(; kwargs...)
Formula(value::ShuntModelFormulation) = value

function Formula(value::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default && value.equivalent_earth === nothing || throw(ArgumentError(
        "shunt_model does not accept equivalent-earth reductions or ordering"))
    return Formula{ID}(; parameters = value.parameters, options = value.options)
end

formula_id(::Formula{ID}) where {ID} = ID
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
description(value::Formula; compact::Bool = false) = description(typeof(value); compact)
formulation_options(value::Formula) = value.options
function Base.NamedTuple(value::Formula)
    (identifier = formula_id(value),
        parameters = value.parameters, options = value.options.data)
end
Base.pairs(::Type{<:Formula}; quantity = nothing) = pairs((;))
