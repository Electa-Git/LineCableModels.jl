import Colors: HSV, RGB

"""
$(TYPEDEF)

Prepare one matrix-grid page for each selected scientific request when comparing
two or more [`LineParameters`](@ref) results.
"""
struct LineParametersBenchmarkPlotDefinition <: PlotBuilder.AbstractPlotDefinition end

const _LineParametersBenchmarkTuple = Tuple{
    LineParameters, LineParameters, Vararg{LineParameters}}

function _comparison_labels(labels, count::Int)
    labels isa Tuple || throw(ArgumentError("legend must be a tuple of labels"))
    length(labels) == count || throw(DimensionMismatch(
        "legend must contain one label for each LineParameters result",
    ))
    all(label -> label isa AbstractString && !isempty(strip(label)), labels) || throw(
        ArgumentError("legend labels must be nonempty strings"),
    )
    return Tuple(String(label) for label in labels)
end

function _comparison_color(index::Int)
    hue = mod(210.0 + 137.50776405003785 * (index - 1), 360.0)
    return RGB(HSV(hue, 0.72, 0.78))
end

function _validate_comparison_inputs(parameters::_LineParametersBenchmarkTuple)
    reference = first(parameters)
    nconductors(reference) > 0 || throw(ArgumentError(
        "line-parameter comparison requires at least one conductor",
    ))
    isempty(frequencies(reference)) && throw(ArgumentError(
        "line-parameter comparison requires at least one frequency",
    ))
    for candidate in Base.tail(parameters)
        size(Z(reference)) == size(Z(candidate)) || throw(DimensionMismatch(
            "all compared Z tensors must have identical dimensions",
        ))
        size(Y(reference)) == size(Y(candidate)) || throw(DimensionMismatch(
            "all compared Y tensors must have identical dimensions",
        ))
        frequencies(reference) == frequencies(candidate) || throw(ArgumentError(
            "all compared frequency vectors must match exactly and in order",
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

function _comparison_input_defaults(::_LineParametersBenchmarkTuple)
    return (;
        requests = (),
        legend = nothing,
        title = nothing,
        labels = nothing,
        freq_unit = :base,
        length_unit = :kilo,
        quantity_units = nothing,
        xscale = :linear,
        yscale = :linear
    )
end

function PlotBuilder.dispatch_on(::Type{LineParametersBenchmarkPlotDefinition})
    return _LineParametersBenchmarkTuple
end

function PlotBuilder.input_defaults(
        ::Type{LineParametersBenchmarkPlotDefinition},
        parameters::_LineParametersBenchmarkTuple
)
    return _comparison_input_defaults(parameters)
end

function PlotBuilder.renderer_defaults(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::_LineParametersBenchmarkTuple
)
    return (; fig_size = (1200, 800))
end

function PlotBuilder.resolve(
        ::Type{LineParametersBenchmarkPlotDefinition},
        parameters,
        request::NamedTuple
)
    input = request.input
    requests = input.requests
    requests isa Tuple || throw(ArgumentError("requests must be a tuple"))
    isempty(requests) && throw(ArgumentError(
        "at least one explicit observable request is required",
    ))
    all(item -> item isa Tuple, requests) || throw(ArgumentError(
        "line comparisons accept explicit observable request tuples",
    ))
    identities = map(_line_request_identity, requests)
    length(unique(identities)) == length(identities) || throw(ArgumentError(
        "line comparisons do not accept duplicate scientific quantities",
    ))
    validate_observables(
        first(parameters),
        NamedTuple{_line_request_names(requests)}(requests)
    )
    foreach(_line_request_indices, requests)
    labels = _comparison_labels(input.legend, length(parameters))
    input.xscale in (:linear, :log10) || throw(
        ArgumentError("xscale must be :linear or :log10"),
    )
    input.yscale in (:linear, :log10) || throw(
        ArgumentError("yscale must be :linear or :log10"),
    )
    request.renderer.fig_size isa Tuple{Int, Int} || throw(
        ArgumentError("fig_size must be a tuple of two integers"),
    )
    all(>(0), request.renderer.fig_size) || throw(
        ArgumentError("fig_size dimensions must be positive"),
    )
    input.title === nothing || input.title isa AbstractString || throw(
        ArgumentError("title must be a string or nothing"),
    )
    input.labels === nothing || input.labels isa Tuple || throw(
        ArgumentError("labels must be a tuple or nothing"),
    )
    input.labels === nothing || length(input.labels) == length(requests) || throw(
        DimensionMismatch("labels must contain one entry per normalized request"),
    )
    input.labels === nothing || all(label -> label isa AbstractString, input.labels) ||
        throw(ArgumentError("labels must contain strings"))
    colors = Tuple(_comparison_color(index) for index in eachindex(parameters))
    return merge(
        request,
        (; input = merge(input, (; requests, legend = labels, colors)))
    )
end

function _comparison_page(configuration, published, request_index, page_index::Int)
    scientific_request = configuration.input.requests[request_index]
    family = _line_request_family(scientific_request)
    parent = _family_parent(family)
    first_observation = first(published).observations[request_index]
    rows, columns, _ = first(published).coordinates[request_index]
    family_symbol = Units.symbol(first_observation.quantity)
    positions = Tuple((local_row, local_column)
    for local_row in eachindex(rows)
    for local_column in eachindex(columns))
    titles = Tuple(
        "$family_symbol[$row,$column] · $(configuration.input.labels === nothing ? Units.label(first_observation.quantity) : configuration.input.labels[request_index])"
    for row in rows
    for column in columns)
    xscales = map(positions) do _
        admitted = _axis_scales(first(published).frequency.values)
        configuration.input.xscale in admitted || throw(DomainError(
            first(published).frequency.values,
            "logarithmic frequency axes require positive finite data and uncertainty bounds"
        ))
        return admitted
    end
    yscales = map(positions) do (local_row, local_column)
        panel_values = collect(Iterators.flatten(
            view(source.observations[request_index].values, local_row, local_column, :)
        for source in published))
        admitted = _axis_scales(panel_values)
        configuration.input.yscale in admitted || throw(DomainError(
            panel_values,
            "logarithmic ordinate axes require positive finite data and uncertainty bounds"
        ))
        return admitted
    end
    attributes = Tuple(
        (;
            xlabelvisible = local_row == length(rows),
            xticklabelsvisible = local_row == length(rows),
            xticksvisible = local_row == length(rows),
            ylabelvisible = local_column == 1,
            yticklabelsvisible = local_column == 1,
            yticksvisible = local_column == 1
        )
    for local_row in eachindex(rows)
    for local_column in eachindex(columns))
    title = configuration.input.title === nothing ?
            "$(Units.label(first_observation.quantity)) comparison" :
            String(configuration.input.title)
    payload = LineDashboardPayload(
        first(published).frequency,
        (scientific_request,),
        Tuple(source.observations[request_index] for source in published),
        (first(published).coordinates[request_index],),
        positions,
        titles,
        xscales,
        yscales,
        attributes,
        configuration.input.legend,
        configuration.input.colors,
        (;
            xscale = configuration.input.xscale,
            yscale = configuration.input.yscale,
            panels = nothing
        )
    )
    return PlotBuilder.PlotPage(
        title,
        configuration.renderer.fig_size,
        (; page = page_index, request = scientific_request),
        payload;
        legend = PlotBuilder.LegendDefinition(),
        export_definition = PlotBuilder.ExportDefinition(
            theme = configuration.renderer.export_theme,
            name = title,
            open_file = configuration.renderer.open_export
        )
    )
end

function PlotBuilder.fetch(
        ::Type{LineParametersBenchmarkPlotDefinition},
        parameters,
        request::NamedTuple
)
    parameters = _validate_comparison_inputs(parameters)
    internal_input = merge(request.input, (; frequencies = nothing, legend = nothing))
    published = Tuple(
        _publish_line_source(parameters[index], internal_input, request.input.requests)
    for index in eachindex(parameters)
    )
    length(first(published).frequency.values) <= 1 && return PlotBuilder.PlotPage[]
    return PlotBuilder.PlotPage[_comparison_page(request, published, request_index, page_index)
                                for (page_index, request_index) in enumerate(eachindex(request.input.requests))]
end
