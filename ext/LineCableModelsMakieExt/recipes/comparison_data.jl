const _LineParametersComparison = Tuple{
    LineParameters, LineParameters, Vararg{LineParameters}}

function _comparison_labels(labels, count::Int)
    labels isa Tuple || labels isa AbstractVector ||
        throw(ArgumentError(
            "series_labels must be a tuple or vector of labels",
        ))
    length(labels) == count || throw(DimensionMismatch(
        "series_labels must contain one label for each LineParameters result",
    ))
    all(label -> label isa AbstractString && !isempty(strip(label)), labels) || throw(
        ArgumentError("series labels must be nonempty strings"),
    )
    return Tuple(String(label) for label in labels)
end

function _validate_comparison_inputs(parameters::_LineParametersComparison)
    reference = first(parameters)
    nconductors(reference) > 0 || throw(ArgumentError(
        "line-parameter comparison requires at least one conductor",
    ))
    isempty(frequencies(reference)) && throw(ArgumentError(
        "line-parameter comparison requires at least one frequency",
    ))
    for candidate in Base.tail(parameters)
        size(Z(reference))[1:2] == size(Z(candidate))[1:2] || throw(DimensionMismatch(
            "all compared Z tensors must have identical matrix dimensions",
        ))
        size(Y(reference))[1:2] == size(Y(candidate))[1:2] || throw(DimensionMismatch(
            "all compared Y tensors must have identical matrix dimensions",
        ))
        basis(reference) === basis(candidate) || throw(ArgumentError(
            "all compared line parameters must use the same basis",
        ))
        domain(reference) === domain(candidate) || throw(ArgumentError(
            "all compared line parameters must use the same domain",
        ))
    end
    return parameters
end

function _prepare_line_comparison(
        parameters::_LineParametersComparison;
        requests,
        series_labels,
        freq_unit = :base,
        length_unit = :kilo,
        quantity_units = nothing,
        clip::Bool = true
)
    _validate_comparison_inputs(parameters)
    labels = _comparison_labels(series_labels, length(parameters))
    published = map(parameters) do parameter
        _prepare_line_observations(
            parameter;
            requests,
            freq_unit,
            length_unit,
            quantity_units,
            clip
        )
    end
    return (; published, labels)
end
