"""
$(TYPEDEF)

Select one semiconducting-screen constitutive relation by its stable literature
identifier.

Each registered formula has one scalar route with the contract
`route(material, frequency, temperature, assumptions) -> Material`. Geometry
and radial series aggregation remain common Engine operations.

$(TYPEDFIELDS)
"""
struct Formula{ID, R, A <: NamedTuple} <: SemiconAdmittanceFormulation
    "Constitutive route for one semiconducting material and one frequency."
    route::R
    "Physical and numerical assumptions of the selected formula."
    assumptions::A
end

"Generic zero-field route tag shared by built-in semicon relations."
struct Functor{ID} end

"Return the stable identifier of a semicon-admittance formula."
formula_id(::Formula{ID}) where {ID} = ID

"Return the assumptions of a semicon-admittance formula."
assumptions(formula::Formula) = formula.assumptions

"Return the default assumptions of a registered semicon-admittance formula."
function assumptions end

"""
$(TYPEDSIGNATURES)

Construct a semicon constitutive relation from its stable identifier.

# Arguments

- `identifier`: Registered literature identifier.

# Keywords

- `route`: Optional complete constitutive route. The route receives a
  [`Material`](@ref), frequency \\[Hz\\], operating temperature \\[°C\\], and the
  formula assumption tuple.
- `kwargs`: Formula-specific assumption overrides.

# Returns

- A concrete semicon-admittance formula.
"""
function Formula(identifier::Symbol; route = nothing, kwargs...)
    Formula(Val(identifier); route, kwargs...)
end

function Formula(::Val{ID}; route = nothing, kwargs...) where {ID}
    tag = Val(ID)
    applicable(assumptions, tag) || throw(ArgumentError(
        "unknown semicon-admittance formula :$ID"
    ))
    defaults = assumptions(tag)
    overrides = (; kwargs...)
    unknown = setdiff(keys(overrides), keys(defaults))
    isempty(unknown) || throw(ArgumentError(
        "unknown assumptions for semicon-admittance formula :$ID: $(collect(unknown))"
    ))
    selected = merge(defaults, overrides)
    selected_route = route === nothing ? Functor{ID}() : route
    return Formula{ID, typeof(selected_route), typeof(selected)}(
        selected_route,
        selected
    )
end

"""
$(TYPEDSIGNATURES)

Construct a fully specified semicon constitutive relation for an external
route.
"""
function Formula(
        identifier::Symbol,
        route,
        values::NamedTuple = (;)
)
    return Formula(Val(identifier), route, values)
end

function Formula(
        ::Val{ID},
        route::R,
        values::A = (;)
) where {ID, R, A <: NamedTuple}
    return Formula{ID, R, A}(route, values)
end

@inline function (formula::Formula)(
        material::Material{T},
        frequency::T,
        temperature::T
) where {T <: Real}
    isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
        frequency,
        "semicon constitutive frequency must be positive and finite"
    ))
    isfinite(temperature) || throw(DomainError(
        temperature,
        "semicon constitutive temperature must be finite"
    ))
    return formula.route(material, frequency, temperature, formula.assumptions)
end

function (formula::Formula)(
        material::Material{T},
        frequency::Real,
        temperature::Real
) where {T <: Real}
    U = promote_type(
        T,
        typeof(float(frequency)),
        typeof(float(temperature))
    )
    return formula(
        convert(Material{U}, material),
        convert(U, float(frequency)),
        convert(U, float(temperature))
    )
end

"Evaluate one registered semicon constitutive relation."
function constitutive(
        formula::Formula,
        material::Material,
        frequency::Real,
        temperature::Real
)
    formula(material, frequency, temperature)
end
