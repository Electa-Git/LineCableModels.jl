# Native renderer attributes follow the same semantic groups as the legends.
function _series_attributes(attributes, count::Integer)
    attributes === nothing && return ntuple(_ -> (;), count)
    attributes isa NamedTuple && return ntuple(_ -> attributes, count)
    attributes isa Union{Tuple, AbstractVector} && length(attributes) == count &&
        all(value -> value isa NamedTuple, attributes) || throw(ArgumentError(
            "series_attributes must be a NamedTuple for all series, or one NamedTuple per legend group"))
    return Tuple(attributes)
end
