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

"Return the stable identifier of a modal-transformation formula."
formula_id(::Formula{ID}) where {ID} = ID

"Return the assumptions of a modal-transformation formula."
assumptions(formula::Formula) = formula.assumptions

"Return the default assumptions of a registered modal formula."
function assumptions end

"Construct phase-to-modal operators for one registered transformation."
function modal_operators end

function Formula(identifier::Symbol; route = nothing, kwargs...)
    return Formula(Val(identifier); route, kwargs...)
end

function Formula(::Val{:default}; route = nothing, kwargs...)
    Formula(Val(DEFAULT); route, kwargs...)
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
    selected_route = route === nothing ? FormulaMethod(tag, modal_operators) : route
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

function Formula(
        ::Val{:default},
        route::R,
        values::A = (;)
) where {R, A <: NamedTuple}
    return Formula(Val(DEFAULT), route, values)
end

(formula::Formula)(parameters) = formula.route(parameters, formula.assumptions)

"""
$(TYPEDEF)

Select the one registered route used by a modal-transformation computation.

The zero-argument constructor selects the module's `:default` route,
`:Chrysochos2014`.

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

function ModalTransformationFormulation()
    return ModalTransformationFormulation(Formula(Val(DEFAULT)))
end

function _modal_formulation(identifier::Symbol, overrides::NamedTuple)
    return ModalTransformationFormulation(Formula(identifier; overrides...))
end

function _modal_formulation(
        selection::FormulaSpec{ID, Order},
        overrides::NamedTuple
) where {ID, Order}
    isempty(overrides) || throw(ArgumentError(
        "formula(...) selections already contain their modal assumptions"
    ))
    Order === :default || throw(ArgumentError(
        "formula order is only valid for equivalent_earth; got :$Order for modal transformation"
    ))
    return ModalTransformationFormulation(
        Formula(Val(ID); selection.overrides...)
    )
end

function _modal_formulation(formula::Formula, overrides::NamedTuple)
    isempty(overrides) || throw(ArgumentError(
        "completed modal formulas cannot receive additional assumptions"
    ))
    return ModalTransformationFormulation(formula)
end

"""
$(TYPEDSIGNATURES)

Select one or more completed modal-transformation formulations.

A scalar symbol, `FormulaSpec`, or completed [`Formula`](@ref) returns one
[`ModalTransformationFormulation`](@ref). An explicit
[`Grid`](@ref LineCableModels.ParametricBuilder.Grid) or
[`Gridspace`](@ref LineCableModels.ParametricBuilder.Gridspace) returns a
`Gridspace{ModalTransformationFormulation}`.
Formula-specific assumptions vary by placing complete `formula(...)`
selections in the finite source.
"""
function ModalTransformationFormulation(
        selection;
        combine::Symbol = :product,
        kwargs...
)
    return _construction(
        ModalTransformationFormulation,
        _modal_formulation,
        (selection, (; kwargs...));
        combine
    )
end

formula_id(formulation::ModalTransformationFormulation) = formula_id(formulation.formula)
assumptions(formulation::ModalTransformationFormulation) = assumptions(formulation.formula)
description(formulation::ModalTransformationFormulation) = description(formulation.formula)
