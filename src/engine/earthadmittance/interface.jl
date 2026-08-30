"""
$(TYPEDEF)

Select an earth-admittance recipe by its literature identifier.

Self and mutual routes have the contract `route(functor, pair)`. The Γ route
has the contract `route(jω, permeability, permittivity)` and returns a named
tuple containing `Γ` and `squared`.

$(TYPEDFIELDS)
"""
struct Formula{
    ID,
    R <: NamedTuple,
    A <: NamedTuple
} <: EarthAdmittanceFormulation
    "Leaf formulas selected for the recipe."
    routes::R
    "Physical assumptions of the recipe."
    assumptions::A
end

"""
$(TYPEDEF)

Store the frequency-dependent values shared by one earth-admittance recipe.

The state has no common physical layout. Each formula owns its state and call
methods.

$(TYPEDFIELDS)
"""
struct Functor{ID, R, S}
    "Concrete leaf routes selected by the formula."
    routes::R
    "Formula-owned frequency state."
    state::S
end

"Return the stable literature identifier of an earth-admittance formula."
formula_id(::Formula{ID}) where {ID} = ID

"Return the concrete leaf routes of an earth-admittance formula."
routes(formula::Formula) = formula.routes

"Return the physical assumptions of an earth-admittance formula."
assumptions(formula::Formula) = formula.assumptions

"Return whether a formula accepts an explicit longitudinal propagation constant."
propagation(formula::Formula{ID}) where {ID} = propagation(Val(ID))

"Return the longitudinal propagation constant \\[1/m\\] stored by a functor."
Γ(functor::Functor) = functor.state.Γ

"Return the default leaf routes of a literature formula."
function routes end

"Return the default physical assumptions of a literature formula."
function assumptions end

"Return the longitudinal-propagation support of a literature formula."
function propagation end

"""
$(TYPEDSIGNATURES)

Construct an earth-admittance formula from a literature identifier.

Keyword arguments replace individual leaf routes. Unspecified routes retain
the selected formula's defaults.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

function Formula(::Val{ID}; kwargs...) where {ID}
    defaults = routes(Val(ID))
    overrides = (; kwargs...)
    unknown = setdiff(keys(overrides), keys(defaults))
    isempty(unknown) || throw(ArgumentError(
        "unknown routes for earth-admittance formula :$ID: $(collect(unknown))"
    ))
    selected = merge(defaults, overrides)
    values = assumptions(Val(ID))
    return Formula{ID, typeof(selected), typeof(values)}(selected, values)
end

"""
$(TYPEDSIGNATURES)

Construct a fully specified earth-admittance formula.

This constructor permits experimental recipes without modifying a built-in
route table. The selected call method determines the formula-owned functor
state.
"""
function Formula(
        ::Val{ID}, selected::R, values::A = (;)
) where {ID, R <: NamedTuple, A <: NamedTuple}
    return Formula{ID, R, A}(selected, values)
end

function Formula(identifier::Symbol, selected::NamedTuple, values::NamedTuple = (;))
    Formula(Val(identifier), selected, values)
end

@inline function (formula::Formula)(rho, epsilon, mu, jω, Γ, segments, thickness)
    return formula(rho, epsilon, mu, jω, Γ, segments)
end

@inline function (functor::Functor)(::Val{:self}, pair)
    return functor.routes.self(functor, pair)
end

@inline function (functor::Functor)(::Val{:mutual}, pair)
    return functor.routes.mutual(functor, pair)
end

function _longitudinal(formula::Formula, ::Nothing, jω, permeability, permittivity)
    return formula.routes.Γ(jω, permeability, permittivity)
end

function _longitudinal(
        formula::Formula,
        value::Complex{T},
        jω,
        permeability,
        permittivity
) where {T <: Real}
    if propagation(formula) === Val(:zero) && !iszero(value)
        throw(ArgumentError(
            "earth-admittance formula :$(formula_id(formula)) fixes Γ to zero"
        ))
    end
    return (Γ = value, squared = value^2)
end
