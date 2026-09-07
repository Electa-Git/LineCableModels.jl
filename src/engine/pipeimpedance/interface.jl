"""
$(TYPEDEF)

Store two circular cable boundaries relative to a common pipe axis.

$(TYPEDFIELDS)
"""
struct Pair{T <: Real}
    "First coaxial-unit index."
    row::Int
    "Second coaxial-unit index."
    column::Int
    "Two Cartesian axis positions relative to the pipe axis [m]."
    positions::Tuple{Tuple{T,T},Tuple{T,T}}
    "Outer radii of the two local coaxial boundaries [m]."
    radii::Tuple{T,T}
end

"""
$(TYPEDEF)

Select a circular-pipe formulation and its numerical sum controls.

$(TYPEDFIELDS)
"""
struct Formula{ID,R <: NamedTuple,A <: NamedTuple} <: PipeImpedanceFormulation
    "Self and mutual cavity routes."
    routes::R
    "Physical and numerical assumptions."
    assumptions::A
end

"""
$(TYPEDEF)

Store one pipe's material and frequency-dependent surface response.

$(TYPEDFIELDS)
"""
struct Functor{ID,R,S}
    "Selected self and mutual cavity routes."
    routes::R
    "Formula-owned numerical values."
    state::S
end

formula_id(::Formula{ID}) where {ID} = ID
routes(formula::Formula) = formula.routes
assumptions(formula::Formula) = formula.assumptions
function routes end
function assumptions end
function pipe_impedance end

Formula(id::Symbol; kwargs...) = Formula(Val(id); kwargs...)
Formula(::Val{:default}; kwargs...) = Formula(Val(DEFAULT); kwargs...)

function Formula(::Val{ID}; rtol::Real=1e-10, max_terms::Integer=1024,
        wall::Symbol=:finite,
        proximity::Symbol=(ID===:Kane1995 || (ID===:DaSilva2006 && wall===:infinite)) ? :Kane1995 : :none,
        kwargs...) where {ID}
    ID in FORMULAS || throw(ArgumentError("unknown pipe-impedance formula :$ID"))
    isfinite(rtol) && zero(rtol)<rtol<one(rtol) ||
        throw(DomainError(rtol, "pipe sum tolerance must lie between zero and one"))
    max_terms>0 || throw(DomainError(max_terms, "pipe sum term limit must be positive"))
    wall in (:finite,:infinite) || throw(ArgumentError("pipe wall must be :finite or :infinite"))
    wall===:infinite && !(ID in (:DaSilva2006,:Hoidalen2013)) &&
        throw(ArgumentError("this source does not supply an infinite-wall selection"))
    proximity in (:none,:Kane1995,:Hoidalen2013) ||
        throw(ArgumentError("unknown core proximity selection"))
    defaults=routes(Val(ID))
    overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) ||
        throw(ArgumentError("unknown pipe-impedance routes"))
    selected=merge(defaults,overrides)
    values=merge(assumptions(Val(ID)),(;rtol,max_terms=Int(max_terms),wall,proximity))
    return Formula{ID,typeof(selected),typeof(values)}(selected,values)
end

@inline (functor::Functor)(::Val{:self},pair::Pair) = functor.routes.self(functor,pair)
@inline (functor::Functor)(::Val{:mutual},pair::Pair) = functor.routes.mutual(functor,pair)
@inline (functor::Functor)(::Val{:inner}) = functor.state.inner
@inline function (functor::Functor)(::Val{:outer})
    functor.state.outer===nothing && throw(ArgumentError("an infinite wall supplies no finite outer-surface term"))
    return functor.state.outer
end
@inline function (functor::Functor)(::Val{:mutual})
    functor.state.transfer===nothing && throw(ArgumentError("an infinite wall supplies no finite through-wall transfer term"))
    return functor.state.transfer
end
