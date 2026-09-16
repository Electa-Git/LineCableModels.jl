"""
$(TYPEDEF)

Select a modal decomposition with model parameters and normalized numerical
controls. The selected concrete type owns `modal_operators`.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions} <: AbstractFormulation
    "Resolved model parameters."
    parameters::P
    "Normalized numerical sections for the modal algorithm."
    options::O
end

formula_id(::Formula{ID}) where {ID} = ID
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)
Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((;))

"""Construct phase-to-modal operators for the selected concrete formulation."""
function modal_operators end

"""
$(TYPEDSIGNATURES)

Construct a modal decomposition with model parameters and numerical controls.
Custom formulations implement `modal_operators` on their own concrete type.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(selected::Formula) = selected

function Formula(::Val{ID}; parameters::NamedTuple=(;), options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    options = options isa NamedTuple ? FormulationOptions(options) : options
    isempty(parameters) || throw(ArgumentError("modal :$ID has no physical parameters"))
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    normalized = formulation_options(FormulaMethod(selected, modal_operators), options)
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
end

function (selected::Formula)(parameters; workspace=(fallback_frequencies=Int[],))
    maps = modal_operators(selected, parameters, selected.parameters, selected.options, workspace)
    maps isa ModalOperators || throw(ArgumentError("modal_operators must return ModalOperators"))
    return _check_operators(maps, parameters)
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "modal formulas cannot consume equivalent_earth"))
    return Formula(Val(ID); parameters=selection.parameters, options=selection.options)
end

"""Expose a selected modal equation and its model and numerical controls."""
Base.NamedTuple(value::Formula) = (identifier=formula_id(value),
    parameters=value.parameters, options=value.options.data)
formulation_options(value::Formula) = value.options

"""
$(TYPEDEF)

Own a modal computation's requested declaration and resolved equation.
The requested declaration is retained before default resolution and numerical
normalization, using the same provenance contract as line-parameter actions.

$(TYPEDFIELDS)
"""
struct ModalTransformationFormulation{F <: AbstractFormulation, D} <: AbstractFormulation
    "Resolved modal decomposition."
    formula::F
    "Requested selection before default resolution and control normalization."
    definition::D
end

ModalTransformationFormulation() = ModalTransformationFormulation(:default)

function _modal_formulation(identifier::Symbol, controls::NamedTuple)
    selection = formula(identifier; controls...)
    return ModalTransformationFormulation(Formula(selection), selection)
end

function _modal_formulation(selection::FormulaDefinition, controls::NamedTuple)
    isempty(controls) || throw(ArgumentError(
        "formula(...) already contains its modal parameters and numerical controls"))
    return ModalTransformationFormulation(Formula(selection), selection)
end

function _modal_formulation(selected::AbstractFormulation, controls::NamedTuple)
    isempty(controls) || throw(ArgumentError(
        "a completed modal formulation cannot receive additional controls"))
    return ModalTransformationFormulation(selected, selected)
end

"""
$(TYPEDSIGNATURES)

Select one or more modal computations. A symbol, formula declaration, or
completed user-owned formulation selects one equation. `Grid` and `Gridspace`
vary complete selections with product or zip composition.
"""
function ModalTransformationFormulation(selection; combine::Symbol=:product, kwargs...)
    return parameterize(ModalTransformationFormulation, _modal_formulation,
        (selection, (; kwargs...)); combine)
end

formula_id(::Type{<:ModalTransformationFormulation}) = :modal
formula_id(::ModalTransformationFormulation) = :modal
description(::Type{<:ModalTransformationFormulation}; compact::Bool=false) = "modal"
description(::ModalTransformationFormulation; compact::Bool=false) = "modal"
description(::Type{ModalTransformationFormulation}, ::Val{:transformation}) = "modal operators"
Base.pairs(::Type{ModalTransformationFormulation}; quantity=nothing) = pairs((transformation=Formula,))
formulation_options(::ModalTransformationFormulation) = FormulationOptions()

function Base.pairs(value::ModalTransformationFormulation; quantity=nothing)
    return pairs(ModalTransformationFormulation,
        (methods=(transformation=value.formula,),
         requested=(transformation=value.definition,), options=(;)); quantity)
end

function Base.pairs(::Type{ModalTransformationFormulation}, retained::NamedTuple; quantity=nothing)
    return pairs(LineParametersFormulation, retained; quantity, owner=ModalTransformationFormulation)
end

function Base.NamedTuple(value::ModalTransformationFormulation)
    return (backend=:modal, requested=(transformation=NamedTuple(value.definition),),
        methods=(transformation=NamedTuple(value.formula),), options=(;))
end
