"""
Interface for a selected temperature-dependent resistivity law.
Concrete selections expose model `parameters` and numerical `options` records.
"""
abstract type TemperatureDependentFormulation <: AbstractFormulation end

"""
$(TYPEDEF)

Select a scalar electrical-resistivity law evaluated from reference material
properties and prescribed temperature. Reference material values are immutable.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: NamedTuple} <: TemperatureDependentFormulation
    "Resolved physical/model parameters."
    parameters::P
    "Normalized numerical sections for this equation."
    options::O
end

"""Return the stable identifier of a temperature-dependent resistivity law."""
formula_id(::Formula{ID}) where {ID} = ID

TextDisplay.@showfields Formula "Formula" selected -> (
    id=formula_id(selected),)

"""Evaluate a temperature-dependent electrical resistivity in ohm meters."""
function temperature_resistivity end

"""
$(TYPEDSIGNATURES)

Construct a temperature-dependent resistivity law with model parameters and
numerical controls. Custom formulations extend `temperature_resistivity`
on their own concrete selection type.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(selected::TemperatureDependentFormulation) = selected

function Formula(::Val{ID}; parameters::NamedTuple=(;), options::NamedTuple=(;)) where {ID}
    isempty(parameters) || throw(ArgumentError("formula :$ID has no configurable model parameters"))
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    binding = FormulaMethod(selected, temperature_resistivity)
    normalized = computation_options(binding, options)
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula(Val(ID); parameters=selection.parameters,
        options=selection.options)
end

"""
$(TYPEDSIGNATURES)

Validate evaluated resistivity \\[Ω·m\\] for a reference material and prescribed
temperature \\[°C\\]. Resistivity must be real and positive; conductors require a
finite value. Passive materials admit infinite resistivity. Return `rho`.
"""
function validate(::Union{Nothing, TemperatureDependentFormulation}, material::Material, temperature::Real, rho)
    isfinite(temperature) || throw(DomainError(temperature,
        "constitutive temperature must be finite"))
    rho isa Real && !isnan(rho) && rho > zero(rho) || throw(DomainError(rho,
        "temperature law must return positive real resistivity for $(material.kind) at $temperature °C"))
    material.kind === :conductor && !isfinite(rho) && throw(DomainError(rho,
        "conductor resistivity must be finite at $temperature °C"))
    return rho
end

@inline function (formula::TemperatureDependentFormulation)(material::Material{T}, temperature::T;
        workspace=nothing) where {T <: Real}
    isfinite(temperature) || throw(DomainError(temperature,
        "constitutive temperature must be finite"))
    rho = temperature_resistivity(
        formula, material, temperature, formula.parameters, formula.options, workspace)
    return validate(formula, material, temperature, rho)
end

function (formula::TemperatureDependentFormulation)(material::Material{T}, temperature::Real;
        workspace=nothing) where {T <: Real}
    U = promote_type(T, typeof(float(temperature)))
    return formula(convert(Material{U}, material), convert(U, temperature); workspace)
end

"""Evaluate a selected resistivity law at prescribed temperature in °C; return Ω·m."""
constitutive(formula::TemperatureDependentFormulation, material::Material, temperature::Real) = formula(material, temperature)

"""Retain reference resistivity in Ω·m when no temperature law is selected."""
constitutive(::Nothing, material::Material, temperature::Real) =
    validate(nothing, material, temperature, material.rho)

"""
$(TYPEDSIGNATURES)

Expose the selected identity, model parameters, and numerical options as a native record.
"""
function Base.NamedTuple(value::Formula)
    return (identifier=formula_id(value), parameters=value.parameters, options=value.options)
end

# Identity-only dispatch also describes retained selections without constructors.
import ...Grammar: formulation_options
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)

"""Iterate the independently selectable child slots admitted by this formula family."""
Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((;))
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
formulation_options(value::Formula) = formulation_options(typeof(value), (parameters=value.parameters, options=value.options))
formulation_options(::Type{<:Formula}, retained::NamedTuple) =
    formulation_options(FormulaDefinition, retained)
