"""
Declare formulation-option defaults for one expression.
"""
function formulation_options(expression::Expression)
    throw(ArgumentError("missing formulation-option defaults for $expression"))
end

function formulation_options(expression::Expression, supplied::FormulationOptions)
    return formulation_options(expression, formulation_options(expression), supplied)
end

"""
$(TYPEDSIGNATURES)

Normalize supplied formulation options against defaults declared by the actual
selected expression. A family with multiple cases projects supplied options to
each consuming expression before calling this constructor. Empty defaults exclude options. Each option's dispatched
normalizer defines its value type, including scalar physical choices and structured
numerical controls.
"""
function formulation_options(expression::Expression, defaults::FormulationOptions, supplied::FormulationOptions)
    default_data, supplied_data = defaults.data, supplied.data
    unknown = filter(name -> !haskey(default_data, name), keys(supplied_data))
    isempty(unknown) || throw(ArgumentError(
        "unused formulation options $(Tuple(unknown)) for $expression"))
    sections = map(keys(default_data)) do name
        default = getproperty(default_data, name)
        explicit = get(supplied_data, name, default isa NamedTuple ? (;) : default)
        formulation_options(expression, Val(name), default, explicit)
    end
    return FormulationOptions(NamedTuple{keys(default_data)}(sections))
end

function formulation_options(expression::Expression, ::Val{Section}, defaults,
        supplied) where {Section}
    throw(ArgumentError("no formulation-option constructor for :$Section of $expression with $(typeof(supplied))"))
end

"""
$(TYPEDSIGNATURES)

Project the options supplied to `formula` onto the expressions it declares. Each supplied
section must be consumed by at least one expression. Each distinct expression receives the
sections its defaults declare, normalized by its own constructor. Return the distinct
expressions in order of first appearance and their `FormulationOptions`, as
`(expressions, options)`.
"""
function formulation_options(formula::AbstractFormulation,
        expressions::Union{Tuple, AbstractVector{<:Expression}})
    supplied = formula.options.data
    identities = unique(expressions)
    defaults = map(formulation_options, identities)
    admitted = union((keys(value.data) for value in defaults)...)
    unknown = setdiff(keys(supplied), admitted)
    isempty(unknown) || throw(ArgumentError(
        "unused formulation options $(Tuple(unknown)) for :$(formula_id(formula))"))
    normalized = map(eachindex(identities)) do index
        declared = defaults[index]
        names = Tuple(intersect(keys(supplied), keys(declared.data)))
        formulation_options(identities[index], declared, FormulationOptions(supplied[names]))
    end
    return (expressions = identities, options = normalized)
end

"""
    formulas(family)

Return the identifiers that a formula family registers, in registration order. `family` is
the family's `Formula` type, such as `EarthImpedance.Formula`. Each family adds one
method on its own `Formula` type.
"""
function formulas end

"""
    bindings(formula, interactions)

Bind each interaction to the expression that `formula` declares for it, together with that
expression's normalized formulation options. Return one record per interaction, in order.
A formula family extends this function for its own formulas and interactions.
"""
function bindings end
import ..LineCableModels: description, formula_id

"""
$(TYPEDSIGNATURES)

Check that the operation of `expression` has a method for its formula and selectors, before
any evaluation. An engine expression takes at least its functor and its workspace after the
selectors, of any type. A backend variant, which takes the backend in their place, does not
count. Return `expression`.

# Errors

- Throws `ArgumentError` naming the formula and the selectors that have no expression.
"""
function validate(expression::Expression)
    F = typeof(expression.selection)
    selectors = map(typeof, expression.arguments)
    hasmethod(expression.method, Tuple{F, selectors..., Any, Any}) && return expression
    isempty(methods(expression.method, Tuple{F, selectors..., Any, Any, Vararg{Any}})) ||
        return expression
    value(::Val{X}) where {X} = X
    throw(ArgumentError("formula :$(formula_id(expression.selection)) has no expression for " *
        join(map(selector -> repr(value(selector)), expression.arguments), ", ")))
end

"""Read a declaration's actual formulation-owned inputs without resolving them."""
formulation_options(value::FormulaDefinition) = value.options
formula_id(source::Pair{<:AbstractFormulation,<:NamedTuple}) = formula_id(first(source))
description(source::Pair{<:AbstractFormulation,<:NamedTuple}; compact::Bool=false) = description(first(source);compact)

"""
$(TYPEDSIGNATURES)

Return the ordered, owner-scoped formula selections and controls relevant to
`quantity`. Pass `nothing` to retain the complete formulation. The formulation's owning
`pairs` methods list them. Descriptions are not identities.
Return `missing` if any selected identity is unavailable. No formula is evaluated.
"""
function formula_id(source::AbstractFormulation, quantity)
    selections = Tuple((scope, formula_id(value), value isa Pair ? last(value) : (;))
        for (scope, value) in pairs(source; quantity))
    return any(selection -> ismissing(selection[2]), selections) ? missing : selections
end
formula_id(::Missing, quantity) = missing

"""Describe a formulation's composition for one existing physical request."""
description(source::AbstractFormulation,request;compact::Bool=true) =
    only(description([source];roles=[:none],quantity=request,compact))

"""Describe a selection in its consuming owner's scientific context."""
description(owner::Type,selected;compact::Bool=false,quantity=nothing) = description(selected;compact)


"""Describe an owner-scoped route using the selected leaf's own description."""
description(scope::Tuple{Type,Tuple},selected;kwargs...) = description(first(scope),last(scope),selected;kwargs...)
function description(owner::Type,route::Tuple, selected; compact::Bool=true, settings::Bool=false)
    name=isempty(route) ? description(owner;compact) : description(owner,Val(first(route)))
    length(route)>1 && (name *= "("*join(string.(Base.tail(route)),",")*")")
    controls=selected isa Pair ? last(selected) : (;)
    value=settings ? "" : isempty(route) ? description(selected;compact) : description(owner,selected;compact)
    isempty(controls) || (value *= (isempty(value) ? "" : " ") *
        (isempty(route) ? sprint(show,controls;context=:compact=>compact) :
            description(FormulaDefinition,controls;compact)))
    return name*"="*value
end

"""Render only the explicit controls admitted by a formula declaration."""
function description(::Type{FormulaDefinition},controls::NamedTuple;compact::Bool=true)
    return "("*join([string(key)*"="*description(Val(key),value;compact)
        for (key,value) in pairs(controls)],", ")*")"
end
description(::Union{Val{:parameters},Val{:options}},value;compact::Bool=true) =
    sprint(show,value;context=:compact=>compact)
