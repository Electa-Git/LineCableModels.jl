# One legend label per series, in a tuple or a vector: a nonempty string or rich text.
function validate(series_labels::Union{Tuple, AbstractVector}, ::Type{Legend}, count::Int)
    length(series_labels) == count || throw(DimensionMismatch(
        "series_labels must contain one label for each LineParameters result",
    ))
    all(
        label -> label isa typeof(Makie.rich("")) ||
                 label isa AbstractString && !isempty(strip(label)), series_labels) || throw(
        ArgumentError("series labels must be nonempty strings"),
    )
    return series_labels
end
function validate(series_labels, ::Type{Legend}, count::Int)
    throw(ArgumentError(
        "series_labels must be a tuple or vector of labels; got $(typeof(series_labels))"))
end
