"Declare numerical defaults for one complete equation binding."
function computation_options(binding::FormulaMethod)
    throw(ArgumentError("missing numerical-default declaration for $binding"))
end

"Declare the numerical requirements of a complete contribution replacement."
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
