"""
$(TYPEDEF)

Select one frequency-dependent earth-material relation by its stable literature
identifier.

Each registered formula has one scalar route with the contract
`route(material, frequency, assumptions) -> EarthMaterial`. The route and its
assumptions participate in the concrete Julia type.

$(TYPEDFIELDS)
"""
struct Formula{ID, R, A <: NamedTuple} <: AbstractFormulation
    "Scalar constitutive route for one soil material and one frequency."
    route::R
    "Numerical and physical assumptions of the relation."
    assumptions::A
end

"Generic zero-field route tag shared by built-in earth-property formulas."
struct Functor{ID} end

"Return the stable literature identifier of an earth-property formula."
formula_id(::Formula{ID}) where {ID} = ID

"Return the assumptions of an earth-property formula."
assumptions(formula::Formula) = formula.assumptions

"Return the default assumptions of a registered earth-property formula."
function assumptions end

"""
$(TYPEDSIGNATURES)

Construct a registered earth-property formula from its literature identifier.

Keyword arguments replace named formula assumptions. Unknown identifiers and
unknown assumptions are rejected before a line-parameter computation begins.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

function Formula(::Val{ID}; route = nothing, kwargs...) where {ID}
    tag = Val(ID)
    applicable(assumptions, tag) || throw(ArgumentError(
        "unknown earth-property formula :$ID"
    ))
    defaults = assumptions(tag)
    overrides = (; kwargs...)
    unknown = setdiff(keys(overrides), keys(defaults))
    isempty(unknown) || throw(ArgumentError(
        "unknown assumptions for earth-property formula :$ID: $(collect(unknown))"
    ))
    values = merge(defaults, overrides)
    selected = route === nothing ? Functor{ID}() : route
    return Formula{ID, typeof(selected), typeof(values)}(selected, values)
end

"""
$(TYPEDSIGNATURES)

Construct a fully specified experimental earth-property relation.

The caller supplies the scalar route and its assumption tuple. This permits an
external formula without adding it to the built-in discovery directory.
"""
function Formula(
        ::Val{ID}, route::R, values::A = (;)
) where {ID, R, A <: NamedTuple}
    return Formula{ID, R, A}(route, values)
end

function Formula(identifier::Symbol, route, values::NamedTuple = (;))
    return Formula(Val(identifier), route, values)
end

@inline function (formula::Formula)(
        material::EarthMaterial{T}, frequency::T
) where {T <: Real}
    isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
        frequency,
        "earth-property evaluation frequency must be positive and finite"
    ))
    return formula.route(material, frequency, formula.assumptions)
end

function (formula::Formula)(material::EarthMaterial{T}, frequency::Real) where {T <: Real}
    U = promote_type(T, typeof(float(frequency)))
    return formula(
        convert(EarthMaterial{U}, material),
        convert(U, float(frequency))
    )
end

"Pass static earth properties through when no constitutive relation is selected."
constitutive(::Nothing, material::EarthMaterial, ::Real) = material

"Evaluate one registered frequency-dependent earth constitutive relation."
constitutive(formula::Formula, material::EarthMaterial, frequency::Real) =
    formula(material, frequency)

"Return vacuum permittivity represented in the scalar type of `value` \\[F/m\\]."
vacuum_permittivity(value) = one(value) * 88541878128 * (one(value) * 10)^(-22)
