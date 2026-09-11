"""
$(TYPEDEF)

Select one frequency-dependent earth-material relation by its stable formula
identifier.

Each registered formula has one scalar route with the contract
`route(material, frequency, parameters, options, workspace) -> EarthMaterial`. The route and its
parameters participate in the concrete Julia type.

$(TYPEDFIELDS)
"""
struct Formula{ID, R, A <: NamedTuple, H <: NamedTuple, O <: NamedTuple} <:
       AbstractFormulation
    "Declared constitutive equation binding for one soil material."
    binding::R
    "Explicit model parameters."
    parameters::A
    "Explicit callable overrides retained for provenance."
    hooks::H
    "Normalized numerical sections for the selected contribution."
    options::O
end

"Return the stable formula identifier of an earth-property formula."
formula_id(::Formula{ID}) where {ID} = ID

"Evaluate one formula-owned frequency-dependent earth material relation."
function earth_material end

"""
$(TYPEDSIGNATURES)

Construct a registered formula with separate model parameters and callable
hooks. `hooks=(contribution=f,)` replaces the complete scalar equation using
the signature `f(material, frequency, parameters, options, workspace) → EarthMaterial`. A complete replacement declares its numerical defaults with
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
    binding = FormulaMethod(Val(ID), earth_material)
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
        material::EarthMaterial{T}, frequency::T; workspace = nothing
) where {T <: Real}
    isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
        frequency,
        "earth-property evaluation frequency must be positive and finite"
    ))
    evaluated = get(formula.hooks, :contribution, formula.binding)(
        material, frequency, formula.parameters, formula.options, workspace)
    evaluated isa EarthMaterial ||
        throw(ArgumentError("an FrequencyDependent contribution must return EarthMaterial"))
    return evaluated
end

function (formula::Formula)(material::EarthMaterial{T}, frequency::Real; workspace = nothing) where {T <:
                                                                                                     Real}
    U = promote_type(T, typeof(float(frequency)))
    return formula(
        convert(EarthMaterial{U}, material),
        convert(U, float(frequency)); workspace
    )
end

"Pass static earth properties through when no constitutive relation is selected."
constitutive(::Nothing, material::EarthMaterial, ::Real) = material

"Evaluate one registered frequency-dependent earth constitutive relation."
function constitutive(formula::Formula, material::EarthMaterial, frequency::Real)
    formula(material, frequency)
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula(Val(ID); parameters = selection.parameters, hooks = selection.hooks,
        options = selection.options)
end

"""
$(TYPEDSIGNATURES)

Expose the selected equation, physical parameters, callable overrides and numerical
options as a native record. Callables are retained unchanged.
"""
function Base.NamedTuple(value::Formula)
    return (identifier=formula_id(value), binding=value.binding,
        parameters=value.parameters, hooks=value.hooks, options=value.options)
end
