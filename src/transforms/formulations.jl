"""
$(TYPEDEF)

Select one modal-decomposition route by its stable identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, R, A <: NamedTuple, H <: NamedTuple, O <: NamedTuple}
    "Declared equation binding returning [`ModalOperators`](@ref)."
    binding::R
    "Explicit physical/model parameters."
    parameters::A
    "Explicit callable overrides retained for provenance."
    hooks::H
    "Normalized numerical sections for the selected modal algorithm."
    options::O
end

"Return the stable identifier of a modal-transformation formula."
formula_id(::Formula{ID}) where {ID} = ID

"Construct phase-to-modal operators for one registered transformation."
function modal_operators end

"""
$(TYPEDSIGNATURES)

Construct a registered formula with separate model parameters and callable
hooks. `hooks=(contribution=f,)` replaces the modal decomposition using
its documented argument and result contract. Unknown fields fail immediately.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(selected::Formula) = selected

function Formula(::Val{ID}; parameters::NamedTuple = (;),
        hooks::NamedTuple = (;), options::NamedTuple = (;)) where {ID}
    ID in FORMULAS || throw(ArgumentError("unknown formula :$ID"))
    isempty(parameters) || throw(ArgumentError("modal :$ID has no physical parameters"))
    isempty(setdiff(keys(hooks), (:contribution,))) ||
        throw(ArgumentError("unknown hooks for :$ID"))
    binding = FormulaMethod(Val(ID), modal_operators)
    selected = get(hooks, :contribution, binding)
    selected === nothing && throw(ArgumentError("a contribution hook must be callable"))
    defaults = haskey(hooks, :contribution) ? computation_options(binding, selected) :
               computation_options(binding)
    normalized = computation_options(binding, defaults, options)
    return Formula{
        ID, typeof(binding), typeof(parameters), typeof(hooks), typeof(normalized)}(
        binding, parameters, hooks, normalized)
end

function (formula::Formula)(parameters; workspace = (fallback_frequencies = Int[],))
    maps = get(formula.hooks, :contribution, formula.binding)(
        parameters, formula.parameters, formula.options, workspace)
    maps isa ModalOperators ||
        throw(ArgumentError("modal contribution must return ModalOperators"))
    return _check_operators(maps, parameters)
end

"""
$(TYPEDEF)

Select the one registered route used by a modal-transformation computation.

The zero-argument constructor selects the package's `:default`
Levenberg–Marquardt modal-tracking route.

$(TYPEDFIELDS)
"""
struct ModalTransformationFormulation{F <: Formula} <: AbstractFormulation
    "Selected modal-decomposition formula."
    formula::F
end

function ModalTransformationFormulation()
    return ModalTransformationFormulation(Formula(:default))
end

function _modal_formulation(identifier::Symbol, overrides::NamedTuple)
    return ModalTransformationFormulation(Formula(identifier; overrides...))
end

function _modal_formulation(
        selection::FormulaDefinition{ID, Order},
        overrides::NamedTuple
) where {ID, Order}
    isempty(overrides) || throw(ArgumentError(
        "formula(...) selections already contain their modal parameters and hooks"
    ))
    Order === :default || throw(ArgumentError(
        "formula order is only valid for equivalent_earth; got :$Order for modal transformation"
    ))
    return ModalTransformationFormulation(
        Formula(selection)
    )
end

function _modal_formulation(formula::Formula, overrides::NamedTuple)
    isempty(overrides) || throw(ArgumentError(
        "completed modal formulas cannot receive additional parameters or hooks"
    ))
    return ModalTransformationFormulation(formula)
end

"""
$(TYPEDSIGNATURES)

Select one or more completed modal-transformation formulations.

A scalar symbol, `FormulaDefinition`, or completed [`Formula`](@ref) returns one
[`ModalTransformationFormulation`](@ref). An explicit
[`Grid`](@ref LineCableModels.ParametricBuilder.Grid) or
[`Gridspace`](@ref LineCableModels.ParametricBuilder.Gridspace) returns a
`Gridspace{ModalTransformationFormulation}`.
Formula-specific parameters and hooks vary by placing complete `formula(...)`
selections in the finite source.
"""
function ModalTransformationFormulation(
        selection;
        combine::Symbol = :product,
        kwargs...
)
    return parameterize(
        ModalTransformationFormulation,
        _modal_formulation,
        (selection, (; kwargs...));
        combine
    )
end

formula_id(formulation::ModalTransformationFormulation) = formula_id(formulation.formula)
description(formulation::ModalTransformationFormulation) = description(formulation.formula)

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing ||
        throw(ArgumentError("modal formulas cannot consume equivalent_earth"))
    return Formula(Val(ID); parameters = selection.parameters, hooks = selection.hooks,
        options = selection.options)
end
