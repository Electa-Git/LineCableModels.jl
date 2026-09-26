"""
Declare formulation-option defaults for one complete equation binding.
"""
function formulation_options(binding::FormulaMethod)
    throw(ArgumentError("missing formulation-option defaults for $binding"))
end

function formulation_options(binding::FormulaMethod, supplied::FormulationOptions)
    return formulation_options(binding, formulation_options(binding), supplied)
end

"""
$(TYPEDSIGNATURES)

Normalize supplied formulation options against defaults declared by the actual
selected equation. A family with multiple cases projects supplied options to
each consuming binding before calling this constructor. Empty defaults admit
no options; they do not declare equation availability. Each option's dispatched
normalizer owns its value type, including scalar physical choices and structured
numerical controls.
"""
function formulation_options(binding::FormulaMethod, defaults::FormulationOptions, supplied::FormulationOptions)
    default_data, supplied_data = defaults.data, supplied.data
    unknown = filter(name -> !haskey(default_data, name), keys(supplied_data))
    isempty(unknown) || throw(ArgumentError(
        "unused formulation options $(Tuple(unknown)) for $binding"))
    sections = map(keys(default_data)) do name
        default = getproperty(default_data, name)
        explicit = get(supplied_data, name, default isa NamedTuple ? (;) : default)
        formulation_options(binding, Val(name), default, explicit)
    end
    return FormulationOptions(NamedTuple{keys(default_data)}(sections))
end

function formulation_options(binding::FormulaMethod, ::Val{Section}, defaults,
        supplied) where {Section}
    throw(ArgumentError("no formulation-option constructor for :$Section of $binding with $(typeof(supplied))"))
end
import ..LineCableModels: description, formula_id

"""Read a declaration's actual formulation-owned inputs without resolving them."""
formulation_options(value::FormulaDefinition) = value.options
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
    selections = Tuple((scope, formula_id(value), value isa Pair ? last(value) : (;))
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
