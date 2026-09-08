"""
$(TYPEDEF)

Select one semiconducting-screen constitutive relation by its stable literature
identifier.

Each registered formula has one scalar route with the contract
`route(material, frequency, temperature, parameters, options, workspace) -> Complex`. The route
returns the material's frequency-evaluated complex admittivity. Geometry and
radial series aggregation remain common Engine operations.

$(TYPEDFIELDS)
"""
struct Formula{ID, R, A <: NamedTuple, H <: NamedTuple, O <: NamedTuple} <:
       SemiconAdmittanceFormulation
    "Constitutive route for one semiconducting material and one frequency."
    binding::R
    "Explicit model parameters."
    parameters::A
    "Explicit callable overrides retained for provenance."
    hooks::H
    "Normalized numerical sections for the selected contribution."
    options::O
end

"Return the stable identifier of a semicon-admittance formula."
formula_id(::Formula{ID}) where {ID} = ID

"Evaluate one formula-owned semiconducting-material constitutive relation."
function semicon_material end

"""
$(TYPEDSIGNATURES)

Construct a registered formula with separate model parameters and callable
hooks. `hooks=(contribution=f,)` replaces the complete scalar equation using
the signature `f(material, frequency, temperature, parameters, options, workspace) → complex admittivity [S/m]`. A complete replacement declares its numerical defaults with
`computation_options(binding, replacement)`. Unknown fields fail immediately.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(selected::Formula) = selected

function Formula(::Val{ID}; parameters::NamedTuple = (;), hooks::NamedTuple = (;),
        options::NamedTuple = (;)) where {ID}
    ID in FORMULAS || throw(ArgumentError("unknown formula :$ID"))
    isempty(parameters) ||
        throw(ArgumentError("formula :$ID has no configurable model parameters"))
    isempty(setdiff(keys(hooks), (:contribution,))) ||
        throw(ArgumentError("unknown hooks for :$ID"))
    binding = FormulaMethod(Val(ID), semicon_material)
    selected = get(hooks, :contribution, binding)
    selected === nothing && throw(ArgumentError("a contribution hook must be callable"))
    defaults = haskey(hooks, :contribution) ? computation_options(binding, selected) :
               computation_options(binding)
    normalized = computation_options(binding, defaults, options)
    return Formula{
        ID, typeof(binding), typeof(parameters), typeof(hooks), typeof(normalized)}(
        binding, parameters, hooks, normalized)
end

@inline function (formula::Formula)(
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
    value = get(formula.hooks, :contribution, formula.binding)(
        material, frequency, temperature, formula.parameters, formula.options, workspace)
    value isa Number && isfinite(value) || throw(DomainError(value,
        "semicon_material must return a finite scalar"))
    return value
end

function (formula::Formula)(
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
    return Formula(Val(ID); parameters = selection.parameters, hooks = selection.hooks,
        options = selection.options)
end
