"""
$(TYPEDEF)

Select one semiconducting-screen constitutive relation by its stable literature
identifier.

Each formula implements `semicon_material(selected, material, frequency,
temperature, parameters, options, workspace)` on its concrete selection type.
It returns the material's frequency-evaluated admittivity [S/m]. Geometry and
radial series aggregation remain common Engine operations.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: NamedTuple} <: SemiconAdmittanceFormulation
    "Resolved physical/model parameters."
    parameters::P
    "Normalized numerical sections for this equation."
    options::O
end

"""
Return the stable identifier of a semicon-admittance formula.
"""
formula_id(::Formula{ID}) where {ID} = ID

"""
Evaluate one formula-owned semiconducting-material constitutive relation.
"""
function semicon_material end

"""
$(TYPEDSIGNATURES)

Construct a selected formulation with model parameters and numerical controls.
Custom formulations extend `semicon_material` on their own concrete selection type.
Unknown controls fail before numerical evaluation.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(selected::SemiconAdmittanceFormulation) = selected

function Formula(::Val{ID}; parameters::NamedTuple=(;), options::NamedTuple=(;)) where {ID}
    isempty(parameters) || throw(ArgumentError("formula :$ID has no configurable model parameters"))
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    binding = FormulaMethod(selected, semicon_material)
    normalized = computation_options(binding, options)
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
end

@inline function (formula::SemiconAdmittanceFormulation)(
        material::Material{T},
        frequency::T,
        temperature::T; workspace = nothing
) where {T <: Real}
    isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
        frequency,
        "semicon constitutive frequency must be positive and finite"
    ))
    isfinite(temperature) || throw(DomainError(
        temperature,
        "semicon constitutive temperature must be finite"
    ))
    value = semicon_material(
        formula, material, frequency, temperature, formula.parameters, formula.options, workspace)
    return validate(formula, T, value)
end

function (formula::SemiconAdmittanceFormulation)(
        material::Material{T},
        frequency::Real,
        temperature::Real; workspace = nothing
) where {T <: Real}
    U = promote_type(
        T,
        typeof(float(frequency)),
        typeof(float(temperature))
    )
    return formula(
        convert(Material{U}, material),
        convert(U, float(frequency)),
        convert(U, float(temperature)); workspace
    )
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
