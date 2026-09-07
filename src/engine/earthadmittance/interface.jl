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

"Return the longitudinal propagation constant prescribed by a formula."
function propagation_constant end

"Evaluate one final earth potential-coefficient route."
function earth_potential_coefficient end

"Evaluate an impedance support route owned by an admittance formula."
function earth_impedance end

@inline function earth_potential_coefficient(
        identifier, ::Val{:self}, functor, pair
)
    return earth_potential_coefficient(identifier, Val(:mutual), functor, pair)
end

"""
$(TYPEDSIGNATURES)

Construct an earth-admittance formula from a literature identifier.

Keyword arguments replace individual leaf routes. Unspecified routes retain
the selected formula's defaults.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

"Retain an explicitly selected formula when the placement is resolved."
Formula(selected::EarthAdmittanceFormulation, ::Val) = selected

function Formula(::Val{ID}; kwargs...) where {ID}
    identifier = Val(ID)
    ID in REGISTERED || throw(ArgumentError(
        "unknown earth-admittance formula :$ID"
    ))
    defaults = routes(identifier)
    overrides = (; kwargs...)
    unknown = setdiff(keys(overrides), keys(defaults))
    isempty(unknown) || throw(ArgumentError(
        "unknown routes for earth-admittance formula :$ID: $(collect(unknown))"
    ))
    selected = merge(defaults, overrides)
    values = assumptions(identifier)
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

function Formula(
        ::Val{:default}, selected::R, values::A = (;)
) where {R <: NamedTuple, A <: NamedTuple}
    return Formula{:default, R, A}(selected, values)
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

"""
$(TYPEDSIGNATURES)

Check the selected earth-admittance route against resolved pair geometry without
evaluating its numerical kernel. Custom route overrides are checked as the
selected callables, not as the original author's default leaves.

# Arguments

- `formula`: Resolved formula; resolve `:default` before this check.
- `pair`: Resolved earth-return interaction with lengths in \\[m\\].

# Returns

- The same `formula`.

# Errors

- Throws an actionable native exception when the selected route cannot
  represent the pair. Numerical integration is not performed by validation.
"""
function validate(formula::Formula, pair::EarthPair)
    selected = pair.row == pair.column ? formula.routes.self : formula.routes.mutual
    validate(pair, selected, formula)
    return formula
end

function validate(formula::Formula{:default}, pair::EarthPair)
    throw(ArgumentError(
        "resolve earth-admittance :default from the problem before validating pair geometry"))
end

function validate(
        pair::EarthPair, route::FormulaMethod{ID, typeof(earth_potential_coefficient)}, formula
) where {ID}
    throw(ArgumentError(
        "earth-admittance route :$ID needs an owner-defined validate(pair, route, formula) method"))
end

"""
$(TYPEDSIGNATURES)

Check the layer inventory consumed by a resolved earth-admittance formula.
Both static layer tuples and frequency-resolved property vectors contain air
first, followed by the earth layers. Formula-specific counts are checked beside
the corresponding equations.

# Arguments

- `formula`: Resolved formula.
- `layer_count`: Number of media, including air.

# Returns

- The same `formula`.

# Errors

- Throws `DimensionMismatch` when the layer count is incompatible.
"""
function validate(formula::Formula, layer_count::Integer)
    layer_count >= 2 || throw(DimensionMismatch(
        "an earth-admittance formula requires air and at least one earth layer"))
    return formula
end

validate(formula::Formula, layers::Union{Tuple, AbstractVector}) =
    validate(formula, length(layers))

function validate(
        formula::Formula, layers::Union{Tuple, AbstractVector},
        thickness::Union{Nothing, AbstractVector}
)
    validate(formula, layers)
    thickness === nothing || length(thickness) == length(layers) ||
        throw(DimensionMismatch("earth-layer thickness and material vectors must align"))
    return formula
end

function validate(
        formula::Formula, resistivity::AbstractVector,
        permittivity::AbstractVector, permeability::AbstractVector,
        thickness::Union{Nothing, AbstractVector}
)
    length(resistivity) == length(permittivity) == length(permeability) ||
        throw(DimensionMismatch("earth-property vectors must have equal lengths"))
    validate(formula, resistivity, thickness)
    return formula
end
