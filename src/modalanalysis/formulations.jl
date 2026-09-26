"""
$(TYPEDEF)

Select a modal decomposition with model parameters and normalized numerical
controls. The modal workspace supplies common numerical storage. The selected
formula implements `decompose!` and allocates additional operation-specific
work through `initialize_buffers`.

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

"""
$(TYPEDSIGNATURES)

Construct a modal decomposition with model parameters and numerical controls.
Custom formulations extend `initialize_buffers` and `decompose!`.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(selected::Formula) = selected

function Formula(::Val{ID}; parameters::NamedTuple=(;), options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    options = options isa NamedTuple ? FormulationOptions(options) : options
    isempty(parameters) || throw(ArgumentError("modal :$ID has no physical parameters"))
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    normalized = formulation_options(FormulaMethod(selected, decompose!), options)
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
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
The requested formula and controls are retained before default resolution and
normalization. The resolved formula and effective controls are stored separately.

$(TYPEDFIELDS)
"""
struct ModalAnalysisFormulation{F <: AbstractFormulation, D} <: AbstractFormulation
    "Resolved modal decomposition."
    formula::F
    "Requested selection before default resolution and control normalization."
    definition::D
end

ModalAnalysisFormulation() = ModalAnalysisFormulation(:default)

function _modal_formulation(identifier::Symbol, controls::NamedTuple)
    selection = formula(identifier; controls...)
    return ModalAnalysisFormulation(Formula(selection), selection)
end

function _modal_formulation(selection::FormulaDefinition, controls::NamedTuple)
    isempty(controls) || throw(ArgumentError(
        "formula(...) already contains its modal parameters and numerical controls"))
    return ModalAnalysisFormulation(Formula(selection), selection)
end

function _modal_formulation(selected::AbstractFormulation, controls::NamedTuple)
    isempty(controls) || throw(ArgumentError(
        "a completed modal formulation cannot receive additional controls"))
    return ModalAnalysisFormulation(selected, selected)
end

function _modal_formulation(selected::ModalAnalysisFormulation, controls::NamedTuple)
    isempty(controls) || throw(ArgumentError(
        "a completed modal formulation cannot receive additional controls"))
    return selected
end

"""
$(TYPEDSIGNATURES)

Select one or more modal computations. A symbol, formula declaration, or
completed user-owned formulation selects one equation. `Grid` and `Gridspace`
vary complete selections with product or zip composition.
"""
function ModalAnalysisFormulation(selection; combine::Symbol=:product, kwargs...)
    return parameterize(ModalAnalysisFormulation, _modal_formulation,
        (selection, (; kwargs...)); combine)
end

formula_id(::Type{<:ModalAnalysisFormulation}) = :modal
formula_id(::ModalAnalysisFormulation) = :modal
description(::Type{<:ModalAnalysisFormulation}; compact::Bool=false) = "modal"
description(::ModalAnalysisFormulation; compact::Bool=false) = "modal"
description(::Type{ModalAnalysisFormulation}, ::Val{:transformation}; compact::Bool=false) = "modal operators"
description(::Type{ModalAnalysisFormulation},selected::AbstractFormulation;
    compact::Bool=false,quantity=nothing) =
    applicable(description,selected) ? description(selected;compact) :
    string(formula_id(selected))
description(::Type{ModalAnalysisFormulation},selected::Pair{<:AbstractFormulation,<:NamedTuple};
    compact::Bool=false,quantity=nothing) =
    description(ModalAnalysisFormulation,first(selected);compact,quantity)
Base.pairs(::Type{ModalAnalysisFormulation}; quantity=nothing) = pairs((transformation=Formula,))
formulation_options(::ModalAnalysisFormulation) = FormulationOptions()

function Base.pairs(value::ModalAnalysisFormulation; quantity=nothing)
    return pairs(ModalAnalysisFormulation,
        (methods=(transformation=value.formula,),
         requested=(transformation=value.definition,), options=(;)); quantity)
end

function Base.pairs(::Type{ModalAnalysisFormulation}, retained::NamedTuple; quantity=nothing)
    return pairs(LineParametersFormulation, retained; quantity, owner=ModalAnalysisFormulation)
end

function Base.NamedTuple(value::ModalAnalysisFormulation)
    return (backend=:modal, requested=(transformation=NamedTuple(value.definition),),
        methods=(transformation=NamedTuple(value.formula),), options=(;))
end
