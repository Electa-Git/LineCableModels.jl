"""
Declare numerical defaults for one complete equation binding.
"""
function computation_options(binding::FormulaMethod)
    throw(ArgumentError("missing numerical-default declaration for $binding"))
end

"""
Declare the numerical requirements of a complete contribution replacement.
"""
function computation_options(binding::FormulaMethod, replacement)
    throw(ArgumentError(
        "complete contribution override $(typeof(replacement)) must declare computation_options(binding, replacement) for $binding"))
end

function computation_options(binding::FormulaMethod, supplied::NamedTuple)
    return computation_options(binding, computation_options(binding), supplied)
end

"""
$(TYPEDSIGNATURES)

Normalize supplied numerical sections against defaults declared by the actual
equation provider. A family with multiple cases projects supplied sections to
each consuming binding before calling this constructor. Empty defaults admit
no numerical options; they do not declare equation availability.
"""
function computation_options(binding::FormulaMethod, defaults::NamedTuple, supplied::NamedTuple)
    unknown = setdiff(keys(supplied), keys(defaults))
    isempty(unknown) || throw(ArgumentError(
        "unused numerical sections $(Tuple(unknown)) for $binding"))
    sections = map(keys(defaults)) do name
        default = getproperty(defaults, name)
        explicit = get(supplied, name, (;))
        default isa NamedTuple && explicit isa NamedTuple || throw(ArgumentError(
            "numerical section :$name must be a NamedTuple"))
        computation_options(binding, Val(name), default, explicit)
    end
    return NamedTuple{keys(defaults)}(sections)
end

function computation_options(binding::FormulaMethod, ::Val{Section}, defaults::NamedTuple,
        supplied::NamedTuple) where {Section}
    throw(ArgumentError("no numerical constructor for section :$Section of $binding"))
end
import ..LineCableModels: FormulaDefinition, description, formula_id

"""Read explicit formula controls without validating or evaluating them."""
function formulation_options(::Type{FormulaDefinition}, retained::NamedTuple)
    return (; (key => retained[key] for key in (:parameters,:hooks,:options,:equivalent_earth)
        if haskey(retained,key) && retained[key] !== nothing &&
            !(retained[key] isa NamedTuple && isempty(retained[key])))...)
end
formulation_options(source::Pair{<:Type,<:NamedTuple}) = formulation_options(first(source),last(source),Val(:retained))
formulation_options(owner::Type,retained::NamedTuple,::Val{:retained}) = formulation_options(owner,retained)
formulation_options(source::Pair{<:AbstractFormulation,<:NamedTuple}) = last(source)
formulation_options(value::FormulaDefinition) = formulation_options(FormulaDefinition,
    (parameters=value.parameters,hooks=value.hooks,options=value.options,equivalent_earth=value.equivalent_earth))
formulation_options(::Symbol) = (;)
formulation_options(value::NamedTuple) = map(formulation_options,value)
formulation_options(::Nothing) = (;)
formulation_options(::Missing) = (;)
formula_id(source::Pair{<:AbstractFormulation,<:NamedTuple}) = formula_id(first(source))
description(source::Pair{<:AbstractFormulation,<:NamedTuple}; compact::Bool=false) = description(first(source);compact)

"""
$(TYPEDSIGNATURES)

Return the ordered, owner-scoped formula selections and controls relevant to
`quantity`. Pass `nothing` to retain the complete formulation. Native and saved
formulations use the same owning `pairs` methods; descriptions are not identities.
Return `missing` if any selected identity is unavailable. No formula is evaluated.
"""
function formula_id(source::Union{AbstractFormulation,Pair{<:Type,<:NamedTuple}}, quantity)
    selections = Tuple((scope, formula_id(value), formulation_options(value))
        for (scope, value) in pairs((source isa Pair ? Tuple(source) : (source,))...; quantity))
    return any(selection -> ismissing(selection[2]), selections) ? missing : selections
end
formula_id(::Missing, quantity) = missing

"""Describe a formulation's composition for one existing physical request."""
description(source::AbstractFormulation,request;compact::Bool=true) =
    only(description([source];roles=[:none],quantity=request,compact))
description(source::Pair{<:Type,<:NamedTuple},request;compact::Bool=true) =
    only(description([source];roles=[:none],quantity=request,compact))

"""Describe a selection in its consuming owner's scientific context."""
description(owner::Type,selected;compact::Bool=false) = description(selected;compact)


"""Describe an owner-scoped route using the selected leaf's own description."""
description(scope::Tuple{Type,Tuple},selected;kwargs...) = description(first(scope),last(scope),selected;kwargs...)
function description(owner::Type,route::Tuple, selected; compact::Bool=true, settings::Bool=false)
    name=isempty(route) ? description(owner;compact) : description(owner,Val(first(route)))
    length(route)>1 && (name *= "("*join(string.(Base.tail(route)),",")*")")
    controls=formulation_options(selected)
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
description(::Union{Val{:parameters},Val{:hooks},Val{:options}},value;compact::Bool=true) =
    sprint(show,value;context=:compact=>compact)
