"""
$(TYPEDEF)

Select a line model with model parameters and normalized numerical
controls. The line-model workspace supplies common numerical storage. The selected
formula implements `fit` and allocates arrays for its operations through `initialize_buffers`.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions} <: AbstractFormulation
    "Resolved model parameters."
    parameters::P
    "Normalized numerical sections for the line-model algorithm."
    options::O
end

"""
$(TYPEDSIGNATURES)

Construct a line model with model parameters and numerical controls.
Custom formulations extend `initialize_buffers` and `fit`.
"""
function Formula{ID}(; parameters::NamedTuple=(;), options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    ID in formulas(Formula) || throw(ArgumentError("unknown line-model formula :$ID"))
    options = FormulationOptions(options)
    isempty(parameters) || throw(ArgumentError("line-model :$ID has no physical parameters"))
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    normalized = formulation_options(Expression(selected, fit), options)
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
end

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(::Val{ID}; kwargs...) where {ID} = Formula{ID}(; kwargs...)
Formula(selected::Formula) = selected

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "line-model formulas cannot consume equivalent_earth"))
    return Formula{ID}(; parameters=selection.parameters, options=selection.options)
end

formula_id(::Formula{ID}) where {ID} = ID
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)
formulation_options(value::Formula) = value.options

"""Expose a selected line-model equation and its model and numerical controls."""
Base.NamedTuple(value::Formula) = (identifier=formula_id(value),
    parameters=value.parameters, options=value.options.data)

Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((;))

"""
$(TYPEDEF)

Own a line-model computation's requested declaration and resolved equation.
The definition preserves the requested formula and controls before default
resolution and normalization. The formula contains the resolved equation and controls.

$(TYPEDFIELDS)
"""
struct LineModelFormulation{F <: AbstractFormulation, D} <: AbstractFormulation
    "Resolved line-model fitting."
    formula::F
    "Requested selection before default resolution and control normalization."
    definition::D
end

function _line_model_formulation(identifier::Symbol, options::NamedTuple)
    selection = formula(identifier; options=options)
    return LineModelFormulation(Formula(selection), selection)
end

function _line_model_formulation(selection::FormulaDefinition, options::NamedTuple)
    isempty(options) || throw(ArgumentError(
        "formula(...) already contains its line-model parameters and numerical controls"))
    return LineModelFormulation(Formula(selection), selection)
end

function _line_model_formulation(selected::AbstractFormulation, options::NamedTuple)
    isempty(options) || throw(ArgumentError(
        "a completed line-model formulation cannot receive additional controls"))
    return LineModelFormulation(selected, selected)
end

function _line_model_formulation(selected::LineModelFormulation, options::NamedTuple)
    isempty(options) || throw(ArgumentError(
        "a completed line-model formulation cannot receive additional controls"))
    return selected
end

"""
$(TYPEDSIGNATURES)

Select one or more line-model computations. A symbol, formula declaration, or
completed user-owned formulation selects one equation. `Grid` and `Gridspace`
vary complete selections with product or zip composition.
"""
function LineModelFormulation(selection=:default; options::NamedTuple=(;), combine::Symbol=:product)
    return parameterize(LineModelFormulation, _line_model_formulation,
        (selection, options); combine)
end

formula_id(::Type{<:LineModelFormulation}) = :linemodel

# The constructor does not lift a Grid inside the options, so the formula would receive it as a
# value. These checks refuse it before a fit fails on it.
function validate(section::NamedTuple, context::Type{<:Formula})
    foreach(value -> validate(value, context), section)
    return section
end
function validate(value::Union{AbstractGrid, Gridspace}, ::Type{<:Formula})
    throw(ArgumentError("Grids inside line-model options do not lift (got $(length(value)) " *
        "alternatives); build one formulation for each set of options"))
end
validate(value, ::Type{<:Formula}) = value
formula_id(::LineModelFormulation) = :linemodel
description(::Type{<:LineModelFormulation}; compact::Bool=false) = "line-model"
description(::LineModelFormulation; compact::Bool=false) = "line-model"
description(::Type{LineModelFormulation}, ::Val{:fitting}; compact::Bool=false) = "line-model fitting"
description(::Type{LineModelFormulation},selected::AbstractFormulation;
    compact::Bool=false,quantity=nothing) = description(selected;compact)
description(::Type{LineModelFormulation},selected::Pair{<:AbstractFormulation,<:NamedTuple};
    compact::Bool=false,quantity=nothing) =
    description(LineModelFormulation,first(selected);compact,quantity)
Base.pairs(::Type{LineModelFormulation}; quantity=nothing) = pairs((fitting=Formula,))
formulation_options(::LineModelFormulation) = FormulationOptions()

function Base.pairs(value::LineModelFormulation; quantity=nothing)
    return pairs(LineModelFormulation,
        (formulas=(fitting=value.formula,),
         requested=(fitting=value.definition,), options=(;)); quantity)
end

function Base.pairs(::Type{LineModelFormulation}, retained::NamedTuple; quantity=nothing)
    return pairs(LineParametersFormulation, retained; quantity, owner=LineModelFormulation)
end

function Base.NamedTuple(value::LineModelFormulation)
    return (backend=:linemodel, requested=(fitting=NamedTuple(value.definition),),
        formulas=(fitting=NamedTuple(value.formula),), options=(;))
end
