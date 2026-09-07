"""
$(TYPEDEF)

Select one pipe-interior potential-coefficient formula.

$(TYPEDFIELDS)
"""
struct Formula{ID,R <: NamedTuple,A <: NamedTuple} <: PipeAdmittanceFormulation
    "Self and mutual potential-coefficient routes."
    routes::R
    "Formula assumptions."
    assumptions::A
end

"""
$(TYPEDEF)

Store a pipe's radius and homogeneous dielectric permittivity.

$(TYPEDFIELDS)
"""
struct Functor{ID,R,T}
    "Selected potential-coefficient routes."
    routes::R
    "Inner pipe radius [m]."
    radius::T
    "Absolute lossless permittivity [F/m]."
    epsilon::T
end

formula_id(::Formula{ID}) where {ID}=ID
routes(formula::Formula)=formula.routes
assumptions(formula::Formula)=formula.assumptions
function routes end
function assumptions end
function pipe_potential_coefficient end

Formula(id::Symbol;kwargs...)=Formula(Val(id);kwargs...)
Formula(::Val{:default};kwargs...)=Formula(Val(DEFAULT);kwargs...)
function Formula(::Val{ID};kwargs...) where {ID}
    ID in FORMULAS || throw(ArgumentError("unknown pipe potential formula :$ID"))
    defaults=routes(Val(ID)); overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) ||
        throw(ArgumentError("unknown pipe potential routes"))
    selected=merge(defaults,overrides); values=assumptions(Val(ID))
    return Formula{ID,typeof(selected),typeof(values)}(selected,values)
end

@inline (f::Functor)(::Val{:self},pair::PipeImpedance.Pair)=f.routes.self(f,pair)
@inline (f::Functor)(::Val{:mutual},pair::PipeImpedance.Pair)=f.routes.mutual(f,pair)
