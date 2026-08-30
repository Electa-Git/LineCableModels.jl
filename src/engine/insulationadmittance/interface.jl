"""
$(TYPEDEF)

Select one insulation-admittance formula by its stable identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, R, A <: NamedTuple} <: InsulationAdmittanceFormulation
    "Single layer-potential route selected for the formula."
    route::R
    "Physical and numerical assumptions of the selected formula."
    assumptions::A
end

"Generic zero-field route tag shared by built-in insulation-admittance formulas."
struct Functor{ID} end

"Return the stable identifier of an insulation-admittance formula."
formula_id(::Formula{ID}) where {ID} = ID

"Return the assumptions of an insulation-admittance formula."
assumptions(formula::Formula) = formula.assumptions

"Return the default assumptions of a registered insulation-admittance formula."
function assumptions end

"""
$(TYPEDSIGNATURES)

Construct an insulation-admittance formula from its stable identifier.

# Arguments

- `identifier`: Registered formula identifier.

# Keywords

- `route`: Optional complete layer-potential route. The route receives inner
  radius \\[m\\], outer radius \\[m\\], resistivity \\[Ω·m\\], relative
  permittivity \\[dimensionless\\], complex angular frequency \\[rad/s\\], and
  the formula assumption tuple.
- `kwargs`: Formula-specific assumption overrides.

# Returns

- A concrete insulation-admittance formula.
"""
Formula(identifier::Symbol; route = nothing, kwargs...) =
    Formula(Val(identifier); route, kwargs...)

function Formula(::Val{ID}; route = nothing, kwargs...) where {ID}
    tag = Val(ID)
    applicable(assumptions, tag) || throw(ArgumentError(
        "unknown insulation-admittance formula :$ID"
    ))
    defaults = assumptions(tag)
    overrides = (; kwargs...)
    unknown = setdiff(keys(overrides), keys(defaults))
    isempty(unknown) || throw(ArgumentError(
        "unknown assumptions for insulation-admittance formula :$ID: $(collect(unknown))"
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

Construct a fully specified insulation-admittance formula for an external
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
        r_in::T,
        r_ex::T,
        rho::T,
        eps_r::T,
        s::Complex{T}
) where {T <: Real}
    return formula.route(r_in, r_ex, rho, eps_r, s, formula.assumptions)
end

@inline function potential_coefficient(
        formula::Formula,
        workspace,
        component::Int,
        s::Complex{T}
) where {T <: Real}
    coefficient = zero(Complex{T})
    @inbounds for layer in workspace.insulator_layer_ranges[component]
        coefficient += formula(
            workspace.r_ins_layer_in[layer],
            workspace.r_ins_layer_ext[layer],
            workspace.rho_ins_layer[layer],
            workspace.eps_ins_layer[layer],
            s
        )
    end
    return coefficient
end
