"""
$(TYPEDEF)

Select one insulation-impedance formula by its stable identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions} <: InsulationImpedanceFormulation
    "Resolved physical/model parameters."
    parameters::P
    "Normalized numerical sections for this equation."
    options::O
end

"""
Return the stable identifier of an insulation-impedance formula.
"""
formula_id(::Formula{ID}) where {ID} = ID

"""
Evaluate one formula-owned insulation-impedance route.
"""
function insulation_impedance end

"""
$(TYPEDSIGNATURES)

Construct a selected formulation with model parameters and numerical controls.
Custom formulations extend `insulation_impedance` on their own concrete selection type.
Unknown controls fail before numerical evaluation.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(selected::InsulationImpedanceFormulation) = selected

function Formula(::Val{ID}; parameters::NamedTuple=(;), options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    options = options isa NamedTuple ? FormulationOptions(options) : options
    isempty(parameters) || throw(ArgumentError("formula :$ID has no configurable model parameters"))
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    binding = FormulaMethod(selected, insulation_impedance)
    normalized = formulation_options(binding, options)
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
end

@inline function (formula::InsulationImpedanceFormulation)(
        r_in::T,
        r_ex::T,
        mu_r::T,
        s::Complex{T}; workspace = nothing
) where {T <: Real}
    isfinite(r_in) && isfinite(r_ex) && zero(T) <= r_in <= r_ex ||
        throw(DomainError((r_in, r_ex), "insulation radii must satisfy 0 ≤ r_in ≤ r_ex [m]"))
    isfinite(mu_r) && mu_r > zero(T) || throw(DomainError(mu_r,
        "relative insulation permeability must be positive and finite"))
    isfinite(s) || throw(DomainError(s, "jω must be finite"))
    value = insulation_impedance(
        formula, r_in, r_ex, mu_r, s, formula.parameters, formula.options, workspace)
    value isa Number && isfinite(value) || throw(DomainError(value,
        "insulation_impedance must return a finite scalar"))
    return value
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula(Val(ID); parameters = selection.parameters,
        options = selection.options)
end

"""
$(TYPEDSIGNATURES)

Expose the selected identity, model parameters, and numerical options as a native record.
"""
function Base.NamedTuple(value::Formula)
    return (identifier=formula_id(value), parameters=value.parameters, options=value.options.data)
end

# Identity-only dispatch also describes retained selections without constructors.
import ...Grammar: formulation_options
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)

"""Iterate the independently selectable child slots admitted by this formula family."""
Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((;))
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
formulation_options(value::Formula) = value.options
