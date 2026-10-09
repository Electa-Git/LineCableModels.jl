"""
$(TYPEDEF)

Select one semiconducting-screen constitutive relation by its stable literature
identifier.

Each formula implements `semicon_material(selected, functor, workspace)` on its concrete
selection type. The input of `functor` holds the material, the frequency, the temperature
and the options. The method returns the material's frequency-evaluated admittivity [S/m]. Geometry and
radial series aggregation remain common Engine operations.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions} <: SemiconAdmittanceFormulation
    "Resolved physical model parameters."
    parameters::P
    "Normalized numerical sections for this equation."
    options::O
end

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
function Formula{ID}(; parameters::NamedTuple=(;), options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    options = options isa NamedTuple ? FormulationOptions(options) : options
    isempty(parameters) || throw(ArgumentError("formula :$ID has no configurable model parameters"))
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    expression = Expression(selected, semicon_material)
    normalized = formulation_options(expression, options)
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
end

@inline function (formula::SemiconAdmittanceFormulation)(
        material::Material{T},
        frequency::T,
        temperature::T; workspace = nothing
) where {T <: Real}
    functor = Functor(formula,
        (; material, frequency, temperature, options = formula.options); workspace)
    value = Expression(formula, semicon_material)(functor, workspace)
    return convert(Complex{T}, validate(value, formula, T))
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

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(::Val{ID}; kwargs...) where {ID} = Formula{ID}(; kwargs...)
Formula(selected::SemiconAdmittanceFormulation) = selected

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula{ID}(; parameters = selection.parameters,
        options = selection.options)
end

"""
Return the stable identifier of a semicon-admittance formula.
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
