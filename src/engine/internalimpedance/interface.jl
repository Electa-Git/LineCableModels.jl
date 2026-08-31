"""
$(TYPEDEF)

Select an internal-impedance recipe by its literature identifier.

The route tuple contains leaf functions with the contract `route(state)`. The
assumption tuple contains formula parameters that do not identify an
interaction. Both tuples participate in the concrete Julia type.

$(TYPEDFIELDS)
"""
struct Formula{
    ID,
    R <: NamedTuple,
    A <: NamedTuple
} <: InternalImpedanceFormulation
    "Leaf formulas for the supported electromagnetic interactions."
    routes::R
    "Physical assumptions of the selected recipe."
    assumptions::A
end

"""
$(TYPEDEF)

Store the values shared by the leaf interactions of one formula call.

The state has no common physical layout: every formula owns its state and call
methods.

$(TYPEDFIELDS)
"""
struct Functor{ID, R, S}
    "Concrete leaf routes selected by the formula."
    routes::R
    "Formula-owned shared numerical state."
    state::S
end

"Return the stable literature identifier of an internal-impedance formula."
formula_id(::Formula{ID}) where {ID} = ID

"Return the concrete leaf routes of an internal-impedance formula."
routes(formula::Formula) = formula.routes

"Return the physical assumptions of an internal-impedance formula."
assumptions(formula::Formula) = formula.assumptions

"Return the default leaf routes of a literature formula."
function routes end

"Return the default physical assumptions of a literature formula."
function assumptions end

"""
$(TYPEDSIGNATURES)

Construct an internal-impedance formula from a literature identifier.

Keyword arguments replace individual leaf routes. Unspecified routes retain
the selected formula's defaults.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

Formula(::Val{:default}; kwargs...) = Formula(Val(DEFAULT); kwargs...)

function Formula(::Val{ID}; kwargs...) where {ID}
    defaults = routes(Val(ID))
    overrides = (; kwargs...)
    unknown = setdiff(keys(overrides), keys(defaults))
    isempty(unknown) || throw(ArgumentError(
        "unknown routes for internal-impedance formula :$ID: $(collect(unknown))"
    ))
    selected = merge(defaults, overrides)
    values = assumptions(Val(ID))
    return Formula{ID, typeof(selected), typeof(values)}(selected, values)
end

"""
$(TYPEDSIGNATURES)

Construct a fully specified internal-impedance formula.

This constructor permits experimental recipes without modifying a built-in
default route table. The caller must provide the corresponding formula call
method.
"""
function Formula(
        ::Val{ID}, selected::R, values::A = (;)
) where {ID, R <: NamedTuple, A <: NamedTuple}
    return Formula{ID, R, A}(selected, values)
end

function Formula(
        ::Val{:default}, selected::R, values::A = (;)
) where {R <: NamedTuple, A <: NamedTuple}
    return Formula(Val(DEFAULT), selected, values)
end

function Formula(identifier::Symbol, selected::NamedTuple, values::NamedTuple = (;))
    Formula(Val(identifier), selected, values)
end
