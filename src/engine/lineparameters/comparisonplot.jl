import Colors: HSV, RGB

"""
$(TYPEDEF)

Build one matrix-grid page for each selected benchmark quantity when comparing two or
more [`LineParameters`](@ref) results. Each page places matrix term ``(i,j)``
at grid position ``(i,j)`` and overlays one solid line per result.
"""
struct LineParametersBenchmarkPlotDefinition <: PlotBuilder.AbstractPlotDefinition end

const _LineParametersBenchmarkTuple = Tuple{
    LineParameters, LineParameters, Vararg{LineParameters}}

function _comparison_page_key(component::Symbol)
    component in _SERIES_COMPONENTS && return LinePageKey{:series, component}()
    component in _SHUNT_COMPONENTS && return LinePageKey{:shunt, component}()
    throw(ArgumentError("unsupported line-parameter component :$component"))
end

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

function _comparison_input_defaults(parameters::_LineParametersBenchmarkTuple)
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

function PlotBuilder.input_kwargs(::Type{LineParametersBenchmarkPlotDefinition})
    return (
        :quantities,
        :legend,
        :freq_unit,
        :length_unit,
        :quantity_units,
        :xscale,
        :yscale
    )
end

PlotBuilder.renderer_kwargs(::Type{LineParametersBenchmarkPlotDefinition}) = (:fig_size,)

function PlotBuilder.input_defaults(
        ::Type{LineParametersBenchmarkPlotDefinition},
        parameters::_LineParametersBenchmarkTuple
)
    return _comparison_input_defaults(parameters)
end

function PlotBuilder.renderer_defaults(
        ::Type{LineParametersBenchmarkPlotDefinition},
        parameters::_LineParametersBenchmarkTuple
)
    return (; fig_size = (1200, 800))
end

function PlotBuilder.resolve(
        ::Type{LineParametersBenchmarkPlotDefinition},
        recipe::PlotBuilder.PlotRecipe
)
    parameters = recipe.object
    input = recipe.input
    components = _resolve_line_components(first(parameters), input.quantities)
    labels = _comparison_labels(input.legend, length(parameters))
    input.xscale in (:linear, :log10) || throw(
        ArgumentError("xscale must be :linear or :log10"),
    )
    input.yscale in (:linear, :log10) || throw(
        ArgumentError("yscale must be :linear or :log10"),
    )
    recipe.renderer.fig_size isa Tuple{Int, Int} || throw(
        ArgumentError("fig_size must be a tuple of two integers"),
    )
    all(>(0), recipe.renderer.fig_size) || throw(
        ArgumentError("fig_size dimensions must be positive"),
    )
    colors = Tuple(_comparison_color(index) for index in eachindex(parameters))
    return PlotBuilder.PlotRecipe(
        LineParametersBenchmarkPlotDefinition,
        parameters,
        merge(input, (; components, legend = labels, colors)),
        recipe.renderer
    )
end

function PlotBuilder.fetch(
        ::Type{LineParametersBenchmarkPlotDefinition},
        recipe::PlotBuilder.PlotRecipe
)
    parameters = _validate_comparison_inputs(recipe.object)
    internal_input = merge(recipe.input, (; frequencies = nothing, con = nothing))
    published = Tuple(
        _publish_line_source(parameters[index], internal_input, recipe.input.components)
    for index in eachindex(parameters)
    )
    frequency = first(published).frequency
    input = merge(recipe.input, (; frequency, published))
    return PlotBuilder.PlotRecipe(
        LineParametersBenchmarkPlotDefinition,
        parameters,
        input,
        recipe.renderer
    )
end

function PlotBuilder._recipe_variant(
        ::Type{LineParametersBenchmarkPlotDefinition},
        recipe::PlotBuilder.PlotRecipe
)
    return Val(:comparison)
end

function PlotBuilder._composition(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe
)
    return Val(:panels)
end

function PlotBuilder._page_keys(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        ::Val{:panels},
        recipe::PlotBuilder.PlotRecipe
)
    length(recipe.input.frequency.values) <= 1 && return ()
    return Tuple(_comparison_page_key(component) for component in recipe.input.components)
end

function PlotBuilder._view_keys(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        ::Val{:panels},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey
)
    first_payload = first(Base.values(first(recipe.input.published).component_payloads))
    count = size(first_payload.values, 1)
    return Tuple((row, column) for row in 1:count for column in 1:count)
end

function PlotBuilder._series_keys(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        ::Val{:panels},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int}
)
    return eachindex(recipe.input.published)
end

function _comparison_values(
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        result_index::Int
)
    component = _line_component(page_key)
    return getproperty(
        recipe.input.published[result_index].component_payloads,
        component
    ).values
end

function PlotBuilder.axis_payload(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        ::Val{:x},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int}
)
    return recipe.input.frequency
end

function PlotBuilder.axis_payload(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int}
)
    return getproperty(
        first(recipe.input.published).component_payloads,
        _line_component(page_key)
    )
end

function PlotBuilder.axis_scale(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        ::Val{:x},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int}
)
    return recipe.input.xscale
end

function PlotBuilder.axis_scale(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int}
)
    return recipe.input.yscale
end

function PlotBuilder.axis_scales(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        ::Val{dim},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int},
        series::Vector{PlotBuilder.SeriesSpec}
) where {dim}
    supports_log = dim === :x ?
                   _supports_log_values(recipe.input.frequency.values) :
                   _supports_log(series, dim)
    return supports_log ? (:linear, :log10) : (:linear,)
end

function PlotBuilder.axis_exponent(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        ::Val{dim},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int},
        series::Vector{PlotBuilder.SeriesSpec}
) where {dim}
    return _finite_exponent(
        (dim === :x ? item.xdata : item.ydata for item in series)
    )
end

function PlotBuilder.series_values(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        ::Val{:x},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int},
        result_index::Int
)
    return recipe.input.frequency.values
end

function PlotBuilder.series_values(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int},
        result_index::Int
)
    values = _comparison_values(recipe, page_key, result_index)
    row, column = view_key
    return collect(view(values, row, column, :))
end

function PlotBuilder.legend_label(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int},
        result_index::Int
)
    return recipe.input.legend[result_index]
end

function PlotBuilder.series_group(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int},
        result_index::Int
)
    return Symbol("line_parameters_$result_index")
end

function PlotBuilder.series_attributes(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int},
        result_index::Int
)
    return (;
        color = recipe.input.colors[result_index],
        linestyle = :solid,
        linewidth = 2
    )
end

function PlotBuilder.default_title(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        ::Nothing
)
    quantity = PlotBuilder.axis_payload(
        LineParametersBenchmarkPlotDefinition,
        Val(:comparison),
        Val(:y),
        recipe,
        page_key,
        (1, 1)
    ).quantity
    return "$(Units.label(quantity)) comparison"
end

function PlotBuilder.default_title(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey{K},
        view_key::Tuple{Int, Int}
) where {K}
    symbol = K === :series ? "Z" : "Y"
    quantity = PlotBuilder.axis_payload(
        LineParametersBenchmarkPlotDefinition,
        Val(:comparison),
        Val(:y),
        recipe,
        page_key,
        view_key
    ).quantity
    label = Units.label(quantity)
    row, column = view_key
    return "$symbol[$row,$column] · $label"
end

function PlotBuilder.view_key(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int}
)
    row, column = view_key
    return (; component = _line_component(page_key), row, column)
end

function PlotBuilder.view_placement(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int}
)
    row, column = view_key
    return PlotBuilder.PlacementSpec(:canvas, PlotBuilder.GridArea(row, column))
end

function PlotBuilder.view_attributes(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key::Tuple{Int, Int}
)
    row, column = view_key
    payload = first(Base.values(first(recipe.input.published).component_payloads))
    count = size(payload.values, 1)
    return (;
        xlabelvisible = row == count,
        xticklabelsvisible = row == count,
        xticksvisible = row == count,
        ylabelvisible = column == 1,
        yticklabelsvisible = column == 1,
        yticksvisible = column == 1
    )
end

function _comparison_layout()
    root = PlotBuilder.GridSpec(
        :root;
        rows = PlotBuilder.AbstractTrackSize[
            PlotBuilder.FixedTrack(36),
            PlotBuilder.RelativeTrack(),
            PlotBuilder.FixedTrack(20)
        ],
        columns = PlotBuilder.AbstractTrackSize[
            PlotBuilder.RelativeTrack(),
            PlotBuilder.ContentTrack()
        ],
        rowgap = 6,
        columngap = 12,
        padding = (20, 20, 28, 28)
    )
    slots = [
        PlotBuilder.SlotSpec(
            :toolbar, :root, PlotBuilder.GridArea(1, 1:2);
            halign = :left, valign = :bottom
        ),
        PlotBuilder.SlotSpec(
            :canvas, :root, PlotBuilder.GridArea(2, 1);
            halign = :stretch, valign = :stretch
        ),
        PlotBuilder.SlotSpec(
            :legend, :root, PlotBuilder.GridArea(2, 2);
            halign = :left, valign = :top
        ),
        PlotBuilder.SlotSpec(
            :status, :root, PlotBuilder.GridArea(3, 1:2);
            halign = :left, valign = :center
        )
    ]
    return PlotBuilder.LayoutSpec(:line_parameters_comparison, [root], slots)
end

function PlotBuilder.layout_spec(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey
)
    return _comparison_layout()
end

function PlotBuilder.default_figsize(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey
)
    return recipe.renderer.fig_size
end

function PlotBuilder.control_spec(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey
)
    return PlotBuilder.ControlSpec(reset = true, export_svg = true, slot = :toolbar)
end

function PlotBuilder.legend_spec(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey
)
    return PlotBuilder.LegendSpec(
        enabled = true,
        interactive = true,
        slot = :legend,
        overflow = :ellipsis
    )
end

function PlotBuilder.colorbar_specs(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey
)
    return PlotBuilder.ColorbarSpec[]
end

function PlotBuilder.status_spec(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey
)
    return PlotBuilder.StatusSpec(enabled = true, initial = "Ready.", slot = :status)
end

function PlotBuilder.page_identity(
        ::Type{LineParametersBenchmarkPlotDefinition},
        ::Val{:comparison},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey
)
    payload = first(Base.values(first(recipe.input.published).component_payloads))
    return (;
        component = _line_component(page_key),
        results = length(recipe.input.published),
        conductors = size(payload.values, 1)
    )
end
