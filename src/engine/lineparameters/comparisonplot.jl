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
        quantities = (),
        legend = nothing,
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
    requests = line_requests(first(parameters), input.quantities)
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
    colors = Tuple(_comparison_color(index) for index in eachindex(parameters))
    return merge(
        request,
        (; input = merge(input, (; requests, legend = labels, colors)))
    )
end

function _comparison_curves(published, request_index, row, column, labels, colors)
    return map(eachindex(published)) do result_index
        observation = published[result_index].observations[request_index]
        LineCurve(
            collect(view(observation.values, row, column, :)),
            labels[result_index],
            Symbol("line_parameters_$result_index"),
            (;
                color = colors[result_index],
                linestyle = :solid,
                linewidth = 2
            )
        )
    end
end

function _comparison_page(configuration, published, request_index, page_index::Int)
    scientific_request = configuration.input.requests[request_index]
    parent = line_parent(scientific_request)
    first_observation = first(published).observations[request_index]
    count = size(first_observation.values, 1)
    family_symbol = Units.symbol(Units.quantity(parent))
    panels = Tuple(
        begin
            curves = _comparison_curves(
                published,
                request_index,
                row,
                column,
                configuration.input.legend,
                configuration.input.colors
            )
            y_observation = _axis_observation(first_observation, curves)
            xscales = _axis_scales(first(published).frequency.values)
            yscales = _axis_scales(y_observation.values)
            configuration.input.xscale in xscales || throw(DomainError(
                first(published).frequency.values,
                "logarithmic frequency axes require positive finite data and uncertainty bounds"
            ))
            configuration.input.yscale in yscales || throw(DomainError(
                y_observation.values,
                "logarithmic ordinate axes require positive finite data and uncertainty bounds"
            ))
            LinePanelPayload(
                scientific_request,
                (row, column),
                "$family_symbol[$row,$column] · $(Units.label(first_observation.quantity))",
                first(published).frequency,
                y_observation,
                curves,
                configuration.input.xscale,
                configuration.input.yscale,
                xscales,
                yscales,
                (;
                    xlabelvisible = row == count,
                    xticklabelsvisible = row == count,
                    xticksvisible = row == count,
                    ylabelvisible = column == 1,
                    yticklabelsvisible = column == 1,
                    yticksvisible = column == 1
                )
            )
        end
    for row in 1:count
    for column in 1:count
    )
    title = "$(Units.label(first_observation.quantity)) comparison"
    payload = LinePagePayload(panels, nothing)
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
    internal_input = merge(request.input, (; frequencies = nothing, con = nothing))
    published = Tuple(
        _publish_line_source(parameters[index], internal_input, request.input.requests)
    for index in eachindex(parameters)
    )
    length(first(published).frequency.values) <= 1 && return PlotBuilder.PlotPage[]
    return PlotBuilder.PlotPage[_comparison_page(request, published, request_index, page_index)
                                for (page_index, request_index) in enumerate(eachindex(request.input.requests))]
end
