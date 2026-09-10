"""
$(TYPEDEF)

Select a scalar electrical-resistivity law evaluated from reference material
properties and prescribed temperature. Reference material values are immutable.

$(TYPEDFIELDS)
"""
struct Formula{ID, R, A <: NamedTuple, H <: NamedTuple, O <: NamedTuple} <: AbstractFormulation
    "Declared scalar constitutive equation."
    binding::R
    "Explicit model parameters."
    parameters::A
    "Concrete callable overrides retained for execution and provenance."
    hooks::H
    "Normalized numerical sections declared by the equation provider."
    options::O
end

"""Return the stable identifier of a temperature-dependent resistivity law."""
formula_id(::Formula{ID}) where {ID} = ID

TextDisplay.@showfields Formula "Formula" selected -> (
    id=formula_id(selected), modified=!isempty(selected.hooks))

"""Evaluate a temperature-dependent electrical resistivity in ohm metres."""
function temperature_resistivity end

"""
$(TYPEDSIGNATURES)

Construct a temperature-dependent resistivity law. A complete replacement
`hooks=(contribution=f,)` has the signature
`f(material, temperature, parameters, options, workspace) -> rho`, with
temperature in °C and resistivity in Ω·m. The replacement declares numerical
defaults through `computation_options(binding, f)`. Unknown fields are rejected.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(selected::Formula) = selected

function Formula(::Val{ID}; parameters::NamedTuple=(;), hooks::NamedTuple=(;),
        options::NamedTuple=(;)) where {ID}
    ID in FORMULAS || throw(ArgumentError("unknown temperature formula :$ID"))
    isempty(parameters) || throw(ArgumentError("temperature formula :$ID has no model parameters"))
    isempty(setdiff(keys(hooks), (:contribution,))) ||
        throw(ArgumentError("unknown temperature-law hooks for :$ID"))
    binding = FormulaMethod(Val(ID), temperature_resistivity)
    selected = get(hooks, :contribution, binding)
    selected === nothing && throw(ArgumentError("a contribution hook must be callable"))
    defaults = haskey(hooks, :contribution) ? computation_options(binding, selected) :
               computation_options(binding)
    normalized = computation_options(binding, defaults, options)
    return Formula{ID, typeof(binding), typeof(parameters), typeof(hooks), typeof(normalized)}(
        binding, parameters, hooks, normalized)
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula(Val(ID); parameters=selection.parameters, hooks=selection.hooks,
        options=selection.options)
end

"""
$(TYPEDSIGNATURES)

Validate evaluated resistivity \\[Ω·m\\] for a reference material and prescribed
temperature \\[°C\\]. Resistivity must be real and positive; conductors require a
finite value. Passive materials admit infinite resistivity. Return `rho`.
"""
function validate(::Union{Nothing, Formula}, material::Material, temperature::Real, rho)
    isfinite(temperature) || throw(DomainError(temperature,
        "constitutive temperature must be finite"))
    rho isa Real && !isnan(rho) && rho > zero(rho) || throw(DomainError(rho,
        "temperature law must return positive real resistivity for $(material.kind) at $temperature °C"))
    material.kind === :conductor && !isfinite(rho) && throw(DomainError(rho,
        "conductor resistivity must be finite at $temperature °C"))
    return rho
end

@inline function (formula::Formula)(material::Material{T}, temperature::T;
        workspace=nothing) where {T <: Real}
    isfinite(temperature) || throw(DomainError(temperature,
        "constitutive temperature must be finite"))
    rho = get(formula.hooks, :contribution, formula.binding)(
        material, temperature, formula.parameters, formula.options, workspace)
    return validate(formula, material, temperature, rho)
end

function (formula::Formula)(material::Material{T}, temperature::Real;
        workspace=nothing) where {T <: Real}
    U = promote_type(T, typeof(float(temperature)))
    return formula(convert(Material{U}, material), convert(U, temperature); workspace)
end

"""Evaluate a selected resistivity law at prescribed temperature in °C; return Ω·m."""
constitutive(formula::Formula, material::Material, temperature::Real) = formula(material, temperature)

"""Retain reference resistivity in Ω·m when no temperature law is selected."""
constitutive(::Nothing, material::Material, temperature::Real) =
    validate(nothing, material, temperature, material.rho)

"""
$(TYPEDSIGNATURES)

Expose the selected equation, physical parameters, callable overrides and numerical
options as a native record. Callables are retained unchanged.
"""
function Base.NamedTuple(value::Formula)
    return (identifier=formula_id(value), binding=value.binding,
        parameters=value.parameters, hooks=value.hooks, options=value.options)
end
