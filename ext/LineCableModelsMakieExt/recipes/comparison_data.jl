function _comparison_labels(labels, count::Int)
    labels isa Tuple || labels isa AbstractVector ||
        throw(ArgumentError(
            "series_labels must be a tuple or vector of labels",
        ))
    length(labels) == count || throw(DimensionMismatch(
        "series_labels must contain one label for each LineParameters result",
    ))
    all(label -> label isa Makie.RichText || label isa AbstractString && !isempty(strip(label)), labels) || throw(
        ArgumentError("series labels must be nonempty strings"),
    )
    return Tuple(labels)
end
