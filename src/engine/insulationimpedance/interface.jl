"""
$(TYPEDEF)

Select one insulation-impedance formula by its stable identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions} <: InsulationImpedanceFormulation
    "Resolved physical model parameters."
    parameters::P
    "Normalized numerical sections for this equation."
    options::O
end

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
function Formula{ID}(; parameters::NamedTuple=(;), options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    options = options isa NamedTuple ? FormulationOptions(options) : options
    isempty(parameters) || throw(ArgumentError("formula :$ID has no configurable model parameters"))
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    expression = Expression(selected, insulation_impedance)
    normalized = formulation_options(expression, options)
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
end

@inline function (formula::InsulationImpedanceFormulation)(
        r_in::T,
        r_ex::T,
        mu_r::T,
        s::Complex{T}; workspace = nothing
) where {T <: Real}
    functor = Functor(formula, (; r_in, r_ex, mu_r, jω = s, options = formula.options);
        workspace)
    value = Expression(formula, insulation_impedance)(functor, workspace)
    value isa Number && isfinite(value) || throw(DomainError(value,
        "insulation_impedance must return a finite scalar"))
    return value
end

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(::Val{ID}; kwargs...) where {ID} = Formula{ID}(; kwargs...)
Formula(selected::InsulationImpedanceFormulation) = selected

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula{ID}(; parameters = selection.parameters,
        options = selection.options)
end

"""
Return the stable identifier of an insulation-impedance formula.
"""
formula_id(::Formula{ID}) where {ID} = ID
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
# Identity-only dispatch also describes retained selections without constructors.
description(value::Formula; compact::Bool = false) = description(typeof(value); compact)
formulation_options(value::Formula) = value.options

"""
$(TYPEDSIGNATURES)

Expose the selected identity, model parameters, and numerical options as a native record.
"""
function Base.NamedTuple(value::Formula)
    return (identifier=formula_id(value), parameters=value.parameters, options=value.options.data)
end

"""Iterate the independently selectable child slots admitted by this formula family."""
Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((;))
