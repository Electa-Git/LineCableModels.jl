"""
$(TYPEDEF)

Select one modal-decomposition route by its stable identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, R, A <: NamedTuple}
    "Single decomposition route returning [`ModalOperators`](@ref)."
    route::R
    "Numerical and physical assumptions of the route."
    assumptions::A
end

"Generic zero-field route tag shared by built-in modal formulas."
struct Functor{ID} end

"Return the stable identifier of a modal-transformation formula."
formula_id(::Formula{ID}) where {ID} = ID

"Return the assumptions of a modal-transformation formula."
assumptions(formula::Formula) = formula.assumptions

"Return the default assumptions of a registered modal formula."
function assumptions end

function Formula(identifier::Symbol; route = nothing, kwargs...)
    return Formula(Val(identifier); route, kwargs...)
end

function Formula(::Val{ID}; route = nothing, kwargs...) where {ID}
    tag = Val(ID)
    applicable(assumptions, tag) || throw(
        ArgumentError("unknown modal-transformation formula :$ID")
    )
    defaults = assumptions(tag)
    overrides = (; kwargs...)
    unknown = setdiff(keys(overrides), keys(defaults))
    isempty(unknown) || throw(ArgumentError(
        "unknown assumptions for modal-transformation formula :$ID: $(collect(unknown))"
    ))
    selected = merge(defaults, overrides)
    selected_route = route === nothing ? Functor{ID}() : route
    return Formula{ID, typeof(selected_route), typeof(selected)}(
        selected_route,
        selected
    )
end

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

(formula::Formula)(parameters) = formula.route(parameters, formula.assumptions)

"""
$(TYPEDEF)

Select the one registered route used by a modal-transformation computation.

$(TYPEDFIELDS)
"""
struct ModalTransformationFormulation{F <: Formula} <: AbstractFormulation
    "Selected modal-decomposition formula."
    formula::F

    function ModalTransformationFormulation(formula::F) where {F <: Formula}
        haskey(formula.assumptions, :tolerance) || throw(ArgumentError(
            "modal-transformation formulas must define a :tolerance assumption"
        ))
        tolerance = formula.assumptions.tolerance
        tolerance isa Real || throw(ArgumentError(
            "modal-transformation tolerance must be real"
        ))
        isfinite(tolerance) && tolerance >= zero(tolerance) || throw(DomainError(
            tolerance,
            "modal-transformation tolerance must be finite and nonnegative"
        ))
        return new{F}(formula)
    end
end

function ModalTransformationFormulation(identifier::Symbol; kwargs...)
    return ModalTransformationFormulation(Formula(identifier; kwargs...))
end

function ModalTransformationFormulation(
        selection::FormulaSpec{ID, Order}
) where {ID, Order}
    Order === :default || throw(ArgumentError(
        "formula order is only valid for equivalent_earth; got :$Order for modal transformation"
    ))
    return ModalTransformationFormulation(
        Formula(Val(ID); selection.overrides...)
    )
end

formula_id(formulation::ModalTransformationFormulation) = formula_id(formulation.formula)
assumptions(formulation::ModalTransformationFormulation) = assumptions(formulation.formula)
description(formulation::ModalTransformationFormulation) = description(formulation.formula)
