"""
$(TYPEDEF)

Select one insulation-impedance formula by its stable identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, R, A <: NamedTuple} <: InsulationImpedanceFormulation
    "Single series-impedance route selected for the formula."
    route::R
    "Physical and numerical assumptions of the selected formula."
    assumptions::A
end

"Return the stable identifier of an insulation-impedance formula."
formula_id(::Formula{ID}) where {ID} = ID

"Return the assumptions of an insulation-impedance formula."
assumptions(formula::Formula) = formula.assumptions

"Return the default assumptions of a registered insulation-impedance formula."
function assumptions end

"Evaluate one formula-owned insulation-impedance route."
function insulation_impedance end

"""
$(TYPEDSIGNATURES)

Construct an insulation-impedance formula from its stable identifier.

# Arguments

- `identifier`: Registered formula identifier.

# Keywords

- `route`: Optional complete layer-impedance route. The route receives inner
  radius \\[m\\], outer radius \\[m\\], relative permeability
  \\[dimensionless\\], complex angular frequency \\[rad/s\\], and the formula
  assumption tuple.
- `kwargs`: Formula-specific assumption overrides.

# Returns

- A concrete insulation-impedance formula.
"""
Formula(identifier::Symbol; route = nothing, kwargs...) =
    Formula(Val(identifier); route, kwargs...)

Formula(::Val{:default}; route = nothing, kwargs...) =
    Formula(Val(DEFAULT); route, kwargs...)

function Formula(::Val{ID}; route = nothing, kwargs...) where {ID}
    tag = Val(ID)
    applicable(assumptions, tag) || throw(ArgumentError(
        "unknown insulation-impedance formula :$ID"
    ))
    defaults = assumptions(tag)
    overrides = (; kwargs...)
    unknown = setdiff(keys(overrides), keys(defaults))
    isempty(unknown) || throw(ArgumentError(
        "unknown assumptions for insulation-impedance formula :$ID: $(collect(unknown))"
    ))
    selected = merge(defaults, overrides)
    selected_route = route === nothing ?
                     FormulaMethod(tag, insulation_impedance) : route
    return Formula{ID, typeof(selected_route), typeof(selected)}(
        selected_route,
        selected
    )
end

"""
$(TYPEDSIGNATURES)

Construct a fully specified insulation-impedance formula for an external
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

function Formula(
        ::Val{:default},
        route::R,
        values::A = (;)
) where {R, A <: NamedTuple}
    return Formula(Val(DEFAULT), route, values)
end

@inline function (formula::Formula)(
        r_in::T,
        r_ex::T,
        mu_r::T,
        s::Complex{T}
) where {T <: Real}
    return formula.route(r_in, r_ex, mu_r, s, formula.assumptions)
end
