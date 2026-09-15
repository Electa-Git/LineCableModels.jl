"""
$(TYPEDEF)

Select a local shunt geometry approximation. Material admittivity remains owned
by the insulation and semiconductor constitutive selections.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: NamedTuple} <: ShuntModelFormulation
    "Explicit model fallback policy."
    parameters::P
    "No callable replacements are declared for this model."
    hooks::NamedTuple{(), Tuple{}}
    "Boundary discretization, quadrature, and optional audit controls."
    options::O
end

"Reference boundary discretization; accuracy depends on geometry and the requested terminal quantity."
const DEFAULT_RESOLUTION = (wire = 64, order = 32, quadrature = 256, modes = 1024)
"Default controls for dimensionless logarithmic-moment integration."
const DEFAULT_INTEGRATION = (rtol = 1e-8, atol = 1e-10, maxevals = 100_000)

"Return the available local shunt model identifiers."
formulas() = (:default, :coaxial, :boundary)
formula_id(::Formula{ID}) where {ID} = ID
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
Formula(value::Formula) = value
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

function Formula(value::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default && value.equivalent_earth === nothing || throw(ArgumentError(
        "shunt_model does not accept equivalent-earth reductions or ordering"))
    return Formula(Val(ID); parameters = value.parameters, hooks = value.hooks, options = value.options)
end

"""
$(TYPEDSIGNATURES)

Construct a cable-local shunt model. `:default` and `:coaxial` use annular
dielectric intervals. `:boundary` resolves eligible lossless open-screen
domains before frequency evaluation.

# Keywords

- `parameters=(;)`: Boundary `fallback=:error` (default) or explicitly
  `:coaxial` after an unsupported boundary assumption or numerical failure.
- `hooks=(;)`: Must be empty; partial replacements of the coupled operator are
  not supported.
- `options=(;)`: Boundary `resolution=(wire=64, order=32, quadrature=256,
  modes=1024)`, `integration=(rtol=1e-8, atol=1e-10, maxevals=100_000)`, and
  `audit=false`. The audit recomputes an independent boundary grid and checks
  derivative step refinement. Coaxial models accept no numerical controls.

# Returns

- A concrete shunt model selection.
"""
function Formula(::Val{ID}; parameters::NamedTuple = (;), hooks::NamedTuple = (;),
        options::NamedTuple = (;)) where {ID}
    ID in (:default, :coaxial) || throw(ArgumentError("unknown shunt model :$ID"))
    isempty(hooks) || throw(ArgumentError("shunt_model has no callable overrides"))
    isempty(parameters) && isempty(options) || throw(ArgumentError(
        "coaxial shunt models accept no parameters or numerical controls"))
    return Formula{ID, typeof(parameters), typeof(options)}(parameters, (;), options)
end

function Formula(::Val{:boundary}; parameters::NamedTuple = (;), hooks::NamedTuple = (;),
        options::NamedTuple = (;))
    isempty(hooks) || throw(ArgumentError("shunt_model has no callable overrides"))
    isempty(setdiff(keys(parameters), (:fallback,))) || throw(ArgumentError(
        "boundary shunt parameters accept only fallback"))
    fallback = get(parameters, :fallback, :error)
    fallback in (:error, :coaxial) ||
        throw(ArgumentError("boundary fallback must be :error or :coaxial"))
    isempty(setdiff(keys(options), (:resolution, :integration, :audit))) ||
        throw(ArgumentError(
            "boundary shunt options accept resolution, integration, and audit"))
    resolution = get(options, :resolution, (;))
    resolution isa NamedTuple &&
    isempty(setdiff(keys(resolution), keys(DEFAULT_RESOLUTION))) ||
        throw(ArgumentError("unknown boundary resolution controls"))
    resolution = merge(DEFAULT_RESOLUTION, resolution)
    all(x -> x isa Integer && !(x isa Bool) && x > 0, values(resolution)) ||
        throw(ArgumentError("boundary resolution controls must be positive integers"))
    resolution.quadrature >= resolution.order+1 || throw(ArgumentError(
        "boundary quadrature must contain at least order+1 nodes"))
    integration = get(options, :integration, (;))
    integration isa NamedTuple &&
    isempty(setdiff(keys(integration), keys(DEFAULT_INTEGRATION))) ||
        throw(ArgumentError("unknown boundary integration controls"))
    integration = merge(DEFAULT_INTEGRATION, integration)
    all(x -> x isa Real && isfinite(x) && x >= 0, (integration.rtol, integration.atol)) &&
    max(integration.rtol, integration.atol) > 0 || throw(ArgumentError(
        "boundary integration requires nonnegative finite tolerances, at least one positive"))
    integration.maxevals isa Integer && !(integration.maxevals isa Bool) &&
    integration.maxevals > 0 ||
        throw(ArgumentError("boundary maxevals must be a positive integer"))
    audit = get(options, :audit, false)
    audit isa Bool || throw(ArgumentError("boundary audit must be Bool"))
    normalized = (; resolution,
        integration = (rtol = Float64(integration.rtol),
            atol = Float64(integration.atol), maxevals = Int(integration.maxevals)),
        audit)
    policy = (; fallback)
    return Formula{:boundary, typeof(policy), typeof(normalized)}(policy, (;), normalized)
end

"""Describe the equivalent annular local shunt approximation."""
function description(::Type{<:Formula{:coaxial}}; compact::Bool = false)
    compact ? "coaxial" :
    "Coaxial annular shunt geometry"
end
"""Describe the default equivalent annular local shunt approximation."""
function description(::Type{<:Formula{:default}}; compact::Bool = false)
    description(Formula{:coaxial}; compact)
end
"""Describe the lossless wire/tape boundary approximation."""
function description(::Type{<:Formula{:boundary}}; compact::Bool = false)
    compact ? "boundary" :
    "Lossless wire/tape boundary shunt geometry"
end
description(value::Formula; compact::Bool = false) = description(typeof(value); compact)
Base.pairs(::Type{<:Formula}; quantity = nothing) = pairs((;))
function formulation_options(value::Formula)
    formulation_options(typeof(value),
        (parameters = value.parameters, hooks = value.hooks, options = value.options))
end
function formulation_options(::Type{<:Formula}, retained::NamedTuple)
    formulation_options(FormulaDefinition, retained)
end
function Base.NamedTuple(value::Formula)
    (identifier = formula_id(value),
        parameters = value.parameters, hooks = value.hooks, options = value.options)
end
