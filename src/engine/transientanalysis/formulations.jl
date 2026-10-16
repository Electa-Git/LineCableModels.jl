"""
$(TYPEDEF)

Select a transient response with model parameters and normalized numerical
controls. The transient workspace supplies common numerical storage. The selected
formula implements `respond!` and allocates additional operation-specific
arrays through `initialize_buffers`.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions} <: AbstractFormulation
    "Resolved model parameters."
    parameters::P
    "Normalized numerical sections for the transient algorithm."
    options::O
end

"""
$(TYPEDSIGNATURES)

Construct a transient response with model parameters and numerical controls.
Custom formulations extend `initialize_buffers` and `respond!`.
"""
function Formula{ID}(; parameters::NamedTuple=(;), options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    ID in formulas(Formula) || throw(ArgumentError("unknown transient formula :$ID"))
    options = FormulationOptions(options)
    isempty(parameters) || throw(ArgumentError("transient :$ID has no physical parameters"))
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    normalized = formulation_options(Expression(selected, respond!), options)
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
end

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(::Val{ID}; kwargs...) where {ID} = Formula{ID}(; kwargs...)
Formula(selected::Formula) = selected

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "transient formulas cannot consume equivalent_earth"))
    return Formula{ID}(; parameters=selection.parameters, options=selection.options)
end

formula_id(::Formula{ID}) where {ID} = ID
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)
formulation_options(value::Formula) = value.options

"""Expose a selected transient equation and its model and numerical controls."""
Base.NamedTuple(value::Formula) = (identifier=formula_id(value),
    parameters=value.parameters, options=value.options.data)

Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((;))

"""
$(TYPEDEF)

Own a transient computation's requested declaration and resolved equation.
The requested formula and controls are retained before default resolution and
normalization. The resolved formula and effective controls are stored separately.

$(TYPEDFIELDS)
"""
struct TransientFormulation{F <: AbstractFormulation, D} <: AbstractFormulation
    "Resolved transient response."
    formula::F
    "Requested selection before default resolution and control normalization."
    definition::D
end

function _transient_formulation(identifier::Symbol, options::NamedTuple)
    selection = formula(identifier; options=options)
    return TransientFormulation(Formula(selection), selection)
end

function _transient_formulation(selection::FormulaDefinition, options::NamedTuple)
    isempty(options) || throw(ArgumentError(
        "formula(...) already contains its transient parameters and numerical controls"))
    return TransientFormulation(Formula(selection), selection)
end

function _transient_formulation(selected::AbstractFormulation, options::NamedTuple)
    isempty(options) || throw(ArgumentError(
        "a completed transient formulation cannot receive additional controls"))
    return TransientFormulation(selected, selected)
end

function _transient_formulation(selected::TransientFormulation, options::NamedTuple)
    isempty(options) || throw(ArgumentError(
        "a completed transient formulation cannot receive additional controls"))
    return selected
end

"""
$(TYPEDSIGNATURES)

Select one or more transient computations. A symbol, formula declaration, or
completed user-owned formulation selects one equation. `Grid` and `Gridspace`
vary complete selections with product or zip composition.
"""
function TransientFormulation(selection=:default; options::NamedTuple=(;), combine::Symbol=:product)
    return parameterize(TransientFormulation, _transient_formulation,
        (selection, options); combine)
end

formula_id(::Type{<:TransientFormulation}) = :transient
formula_id(::TransientFormulation) = :transient
description(::Type{<:TransientFormulation}; compact::Bool=false) = "transient"
description(::TransientFormulation; compact::Bool=false) = "transient"
description(::Type{TransientFormulation}, ::Val{:response}; compact::Bool=false) = "transient response"
description(::Type{TransientFormulation},selected::AbstractFormulation;
    compact::Bool=false,quantity=nothing) = description(selected;compact)
description(::Type{TransientFormulation},selected::Pair{<:AbstractFormulation,<:NamedTuple};
    compact::Bool=false,quantity=nothing) =
    description(TransientFormulation,first(selected);compact,quantity)
Base.pairs(::Type{TransientFormulation}; quantity=nothing) = pairs((response=Formula,))
formulation_options(::TransientFormulation) = FormulationOptions()

function Base.pairs(value::TransientFormulation; quantity=nothing)
    return pairs(TransientFormulation,
        (formulas=(response=value.formula,),
         requested=(response=value.definition,), options=(;)); quantity)
end

function Base.pairs(::Type{TransientFormulation}, retained::NamedTuple; quantity=nothing)
    return pairs(LineParametersFormulation, retained; quantity, owner=TransientFormulation)
end

function Base.NamedTuple(value::TransientFormulation)
    return (backend=:transient, requested=(response=NamedTuple(value.definition),),
        formulas=(response=NamedTuple(value.formula),), options=(;))
end
