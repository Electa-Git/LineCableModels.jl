# Native renderer attributes follow the same semantic groups as the legends.
function _series_attributes(attributes, count::Integer)
    attributes === nothing && return fill((;), count)
    return collect(validate(attributes, Makie.Plot, count))
end

# One NamedTuple of native plot attributes per series, in a tuple or a vector.
function validate(series_attributes::Union{Tuple, AbstractVector}, ::Type{Makie.Plot},
        count::Integer)
    length(series_attributes) == count &&
        all(value -> value isa NamedTuple, series_attributes) || throw(ArgumentError(
        "series_attributes must contain one NamedTuple per legend group"))
    return series_attributes
end
function validate(series_attributes, ::Type{Makie.Plot}, count::Integer)
    throw(ArgumentError("series_attributes must be a tuple or vector of NamedTuples, " *
        "one per legend group; got $(typeof(series_attributes))"))
end
