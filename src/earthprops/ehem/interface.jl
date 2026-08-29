"Abstract equivalent homogeneous-earth rule."
abstract type AbstractRule <: AbstractFormulation end

"Abstract ordering of material frequency dependence and EHEM reduction."
abstract type AbstractSequence <: AbstractFormulation end

"""
$(TYPEDEF)

Select one equivalent homogeneous-earth rule by its stable literature
identifier.

The rule receives the physical pair layout, evaluated layer properties, the
static earth geometry, the interaction pair, and frequency. It returns one
artificial homogeneous [`EarthMaterial`](@ref). Formula assumptions participate
in the concrete Julia type.

$(TYPEDFIELDS)
"""
struct Formula{ID, R, A <: NamedTuple} <: AbstractRule
    "Callable EHEM route."
    route::R
    "Numerical and reconstruction assumptions."
    assumptions::A
end

"Zero-field route tag shared by built-in EHEM formulas."
struct Functor{ID} end

"""
$(TYPEDEF)

Apply material frequency dependence to every physical layer before the EHEM
rule constructs an equivalent material.

$(TYPEDFIELDS)
"""
struct AfterFD{R <: AbstractRule} <: AbstractSequence
    "Equivalent homogeneous-earth rule."
    rule::R
end

"""
$(TYPEDEF)

Construct an equivalent material from static layers before applying the
selected material frequency dependence to that artificial material.

$(TYPEDFIELDS)
"""
struct BeforeFD{R <: AbstractRule} <: AbstractSequence
    "Equivalent homogeneous-earth rule."
    rule::R
end

"Return the rule stored by an EHEM composition."
rule(sequence::AbstractSequence) = sequence.rule

"Return the stable literature identifier of an EHEM formula."
formula_id(::Formula{ID}) where {ID} = ID

"Return the assumptions of an EHEM formula."
assumptions(formula::Formula) = formula.assumptions

"Return the default assumptions of a registered EHEM formula."
function assumptions end

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

function Formula(::Val{ID}; route = nothing, kwargs...) where {ID}
    tag = Val(ID)
    applicable(assumptions, tag) || throw(ArgumentError(
        "unknown EHEM formula :$ID"
    ))
    defaults = assumptions(tag)
    overrides = (; kwargs...)
    unknown = setdiff(keys(overrides), keys(defaults))
    isempty(unknown) || throw(ArgumentError(
        "unknown assumptions for EHEM formula :$ID: $(collect(unknown))"
    ))
    values = merge(defaults, overrides)
    selected = route === nothing ? Functor{ID}() : route
    return Formula{ID, typeof(selected), typeof(values)}(selected, values)
end

function Formula(
        ::Val{ID}, route::R, values::A = (;)
) where {ID, R, A <: NamedTuple}
    return Formula{ID, R, A}(route, values)
end

Formula(identifier::Symbol, route, values::NamedTuple = (;)) =
    Formula(Val(identifier), route, values)

AfterFD(identifier::Symbol; kwargs...) = AfterFD(Formula(identifier; kwargs...))
BeforeFD(identifier::Symbol; kwargs...) = BeforeFD(Formula(identifier; kwargs...))

description(sequence::AfterFD) = "$(description(sequence.rule)) after layerwise FD"
description(sequence::BeforeFD) = "$(description(sequence.rule)) before layerwise FD"

function _check(
        rho::AbstractVector,
        eps_r::AbstractVector,
        mu_r::AbstractVector,
        model::EarthModel
)
    count = length(model.layers)
    length(rho) == length(eps_r) == length(mu_r) == count ||
        throw(DimensionMismatch(
            "EHEM property vectors must align with the $count earth-model layers"
        ))
    return nothing
end

function _index(model::EarthModel, layer::Int)
    selected = layer == -1 ? length(model.layers) : layer
    selected in 2:length(model.layers) || throw(BoundsError(model.layers, selected))
    return selected
end

@inline function _material(rho, eps_r, mu_r, model::EarthModel, layer::Int)
    selected = _index(model, layer)
    return EarthMaterial(rho[selected], eps_r[selected], mu_r[selected])
end

function _horizontal(model::EarthModel, identifier::Symbol)
    model.vertical_layers && throw(ArgumentError(
        "EHEM formula :$identifier supports horizontal earth layers only"
    ))
    return nothing
end

function _nonmagnetic(mu_r, identifier::Symbol)
    all(isone, @view mu_r[2:end]) || throw(ArgumentError(
        "EHEM formula :$identifier assumes relative permeability equal to one"
    ))
    return nothing
end

@inline function (formula::Formula)(
        layout::Val,
        rho::AbstractVector,
        eps_r::AbstractVector,
        mu_r::AbstractVector,
        model::EarthModel,
        pair,
        frequency::Real
)
    _check(rho, eps_r, mu_r, model)
    isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
        frequency,
        "EHEM evaluation frequency must be positive and finite"
    ))
    return formula.route(
        layout, rho, eps_r, mu_r, model, pair, frequency, formula.assumptions
    )
end

function (::Functor{ID})(::Val{Layout}, args...) where {ID, Layout}
    throw(ArgumentError(
        "EHEM formula :$ID does not support :$Layout conductor pairs"
    ))
end
