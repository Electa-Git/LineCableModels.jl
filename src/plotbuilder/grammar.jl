const COMMON_RENDERER_KWARGS = (:export_theme, :open_export, :layout)

"""
    dispatch_on(::Type{S})

Return the domain type accepted by plot specification `S`.
"""
dispatch_on(::Type{S}) where {S <: AbstractPlotDefinition} = Any

"""
    input_kwargs(::Type{S})

Return semantic keyword names accepted by plot specification `S`.
"""
input_kwargs(::Type{S}) where {S <: AbstractPlotDefinition} = ()

"""
    renderer_kwargs(::Type{S})

Return recipe-specific renderer keyword names accepted by plot specification `S`.
"""
renderer_kwargs(::Type{S}) where {S <: AbstractPlotDefinition} = ()

"""
    input_defaults(::Type{S}, object)

Return defaults for every name declared by `input_kwargs(S)`.
"""
input_defaults(::Type{S}, object) where {S <: AbstractPlotDefinition} = (;)

"""
    renderer_defaults(::Type{S}, object)

Return defaults for every name declared by `renderer_kwargs(S)`.
"""
renderer_defaults(::Type{S}, object) where {S <: AbstractPlotDefinition} = (;)

function _symbol_tuple(::Type{S}, values, accessor::Symbol) where {S <:
                                                                   AbstractPlotDefinition}
    values isa Tuple || throw(
        ArgumentError("$accessor($S) must return a tuple of symbols"),
    )
    all(value -> value isa Symbol, values) || throw(
        ArgumentError("$accessor($S) must return a tuple of symbols"),
    )
    length(unique(values)) == length(values) || throw(
        ArgumentError("$accessor($S) cannot declare duplicate keywords"),
    )
    return values
end

function _select_kwargs(kwargs::NamedTuple, names::Tuple)
    selected = Tuple(pair for pair in pairs(kwargs) if first(pair) in names)
    return NamedTuple(selected)
end

function _validate_defaults(::Type{S}, defaults::NamedTuple, names::Tuple,
        accessor::Symbol) where {
        S <: AbstractPlotDefinition,
}
    actual = Tuple(keys(defaults))
    Set(actual) == Set(names) || throw(
        ArgumentError(
        "$accessor($S) must define exactly the declared keywords; " *
        "declared $(collect(names)), received $(collect(actual))"
    ),
    )
    return defaults
end

"""
    parse_kwargs(::Type{S}, object, kwargs)

Validate and split caller input into the semantic and renderer options declared
by a recipe. Unsupported keywords are errors.
"""
function parse_kwargs(
        ::Type{S},
        object,
        kwargs::NamedTuple
) where {S <: AbstractPlotDefinition}
    input_names = _symbol_tuple(S, input_kwargs(S), :input_kwargs)
    renderer_names = _symbol_tuple(S, renderer_kwargs(S), :renderer_kwargs)
    collisions = intersect(input_names, renderer_names)
    isempty(collisions) || throw(
        ArgumentError("$S declares keywords in both input and renderer options: $(join(collisions, ", "))"),
    )
    common_collisions = intersect((input_names..., renderer_names...), COMMON_RENDERER_KWARGS)
    isempty(common_collisions) || throw(
        ArgumentError("$S redeclares common renderer keywords: $(join(common_collisions, ", "))"),
    )
    allowed_renderer = (COMMON_RENDERER_KWARGS..., renderer_names...)
    allowed = (input_names..., allowed_renderer...)
    unsupported = Tuple(name for name in keys(kwargs) if name ∉ allowed)
    isempty(unsupported) || throw(
        ArgumentError("unsupported plot keyword(s) for $S: $(join(unsupported, ", "))"),
    )

    declared_input = input_defaults(S, object)
    declared_input isa NamedTuple || throw(
        ArgumentError("input_defaults($S) must return a NamedTuple"),
    )
    declared_renderer = renderer_defaults(S, object)
    declared_renderer isa NamedTuple || throw(
        ArgumentError("renderer_defaults($S) must return a NamedTuple"),
    )
    _validate_defaults(S, declared_input, input_names, :input_defaults)
    _validate_defaults(S, declared_renderer, renderer_names, :renderer_defaults)

    input = merge(declared_input, _select_kwargs(kwargs, input_names))
    renderer = merge(
        (; export_theme = :default, open_export = true, layout = nothing),
        declared_renderer,
        _select_kwargs(kwargs, allowed_renderer)
    )
    _validate_export_theme(renderer.export_theme)
    renderer.open_export isa Bool || throw(ArgumentError("open_export must be Bool"))
    renderer.layout isa Union{Nothing, Symbol, LayoutSpec} || throw(
        ArgumentError("layout must be nothing, a preset symbol, or LayoutSpec"),
    )
    return PlotRecipe(object, input, renderer)
end

function parse_kwargs(::Type{S}, object; kwargs...) where {S <: AbstractPlotDefinition}
    parse_kwargs(S, object, (; kwargs...))
end

"""
    resolve_input(::Type{S}, recipe)

Validate and enrich parsed recipe input before materialization.
"""
resolve_input(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition} = recipe

"""
    recipe_mode(::Type{S}, recipe)

Return the value-dispatched plotting mode for a resolved recipe.
"""
function recipe_mode(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition}
    Val(:default)
end

"""
    grouping_mode(::Type{S}, mode, recipe)

Return the value-dispatched grouping mode for a resolved recipe.
"""
function grouping_mode(
        ::Type{S},
        mode::Val,
        recipe::PlotRecipe
) where {S <: AbstractPlotDefinition}
    Val(:overlay)
end

"""
    page_facets(::Type{S}, mode, recipe)

Return semantic page facets for recipes using `:faceted_pages`.
"""
function page_facets(
        ::Type{S},
        mode::Val,
        recipe::PlotRecipe
) where {S <: AbstractPlotDefinition}
    (nothing,)
end

"""
    group_facets(::Type{S}, mode, recipe, page_key)

Return semantic series facets for the selected recipe page.
"""
function group_facets(
        ::Type{S},
        mode::Val,
        recipe::PlotRecipe,
        page_key
) where {S <: AbstractPlotDefinition}
    (nothing,)
end

function page_keys(::Type{S}, mode::Val, ::Val{:overlay},
        recipe::PlotRecipe) where {
        S <: AbstractPlotDefinition,
}
    (nothing,)
end
function page_keys(::Type{S}, mode::Val, ::Val{:panels},
        recipe::PlotRecipe) where {
        S <: AbstractPlotDefinition,
}
    (nothing,)
end
function page_keys(::Type{S}, mode::Val, ::Val{:pages}, recipe::PlotRecipe) where {
        S <: AbstractPlotDefinition,
}
    group_facets(S, mode, recipe, nothing)
end
function page_keys(::Type{S}, mode::Val, ::Val{:faceted_pages},
        recipe::PlotRecipe) where {
        S <: AbstractPlotDefinition,
}
    page_facets(S, mode, recipe)
end
function page_keys(::Type{S}, mode::Val, ::Val{:empty}, recipe::PlotRecipe) where {
        S <: AbstractPlotDefinition,
}
    (nothing,)
end
function page_keys(
        ::Type{S},
        mode::Val,
        grouping::Val,
        recipe::PlotRecipe
) where {S <: AbstractPlotDefinition}
    grouping_name = only(typeof(grouping).parameters)
    throw(
        ArgumentError(
        "unsupported grouping mode :$grouping_name for $S; " *
        "specialize PlotBuilder.page_keys for this Val mode"
    ),
    )
end

function view_keys(::Type{S}, mode::Val, ::Val{:overlay}, recipe::PlotRecipe,
        page_key) where {
        S <: AbstractPlotDefinition,
}
    (nothing,)
end
function view_keys(::Type{S}, mode::Val, ::Val{:panels}, recipe::PlotRecipe,
        page_key) where {
        S <: AbstractPlotDefinition,
}
    group_facets(S, mode, recipe, page_key)
end
function view_keys(::Type{S}, mode::Val, ::Val{:pages}, recipe::PlotRecipe,
        page_key) where {
        S <: AbstractPlotDefinition,
}
    (nothing,)
end
function view_keys(::Type{S}, mode::Val, ::Val{:faceted_pages},
        recipe::PlotRecipe, page_key) where {
        S <: AbstractPlotDefinition,
}
    (nothing,)
end
function view_keys(::Type{S}, mode::Val, ::Val{:empty}, recipe::PlotRecipe,
        page_key) where {
        S <: AbstractPlotDefinition,
}
    ()
end

function series_keys(::Type{S}, mode::Val, ::Val{:overlay}, recipe::PlotRecipe,
        page_key, view_key) where {S <: AbstractPlotDefinition}
    group_facets(S, mode, recipe, page_key)
end
function series_keys(::Type{S}, mode::Val, ::Val{:panels}, recipe::PlotRecipe,
        page_key, view_key) where {S <: AbstractPlotDefinition}
    (view_key,)
end
function series_keys(::Type{S}, mode::Val, ::Val{:pages}, recipe::PlotRecipe,
        page_key, view_key) where {S <: AbstractPlotDefinition}
    (page_key,)
end
function series_keys(::Type{S}, mode::Val, ::Val{:faceted_pages}, recipe::PlotRecipe,
        page_key, view_key) where {S <: AbstractPlotDefinition}
    group_facets(S, mode, recipe, page_key)
end

"""
    geom_axes(::Type{S}, mode, recipe, page_key, view_key)

Return the dimensions used by a recipe view.
"""
function geom_axes(::Type{S}, mode::Val, recipe::PlotRecipe,
        page_key, view_key) where {
        S <: AbstractPlotDefinition,
}
    (:x, :y)
end

"""
    axis_quantity(::Type{S}, dim, recipe)

Return the semantic quantity tag for one axis.
"""
function axis_quantity(::Type{S}, dim::Val, recipe::PlotRecipe) where {S <:
                                                                       AbstractPlotDefinition}
    QuantityTag{:unknown}()
end
function axis_quantity(
        ::Type{S}, mode::Val, dim::Val, recipe::PlotRecipe,
        page_key, view_key
) where {S <: AbstractPlotDefinition}
    axis_quantity(S, dim, recipe)
end

"""
    axis_unit(::Type{S}, dim, quantity, recipe)

Return the display units for one axis quantity.
"""
function axis_unit(::Type{S}, dim::Val, quantity::QuantityTag,
        recipe::PlotRecipe) where {
        S <: AbstractPlotDefinition,
}
    display_unit(quantity)
end
function axis_unit(
        ::Type{S}, mode::Val, dim::Val, quantity::QuantityTag,
        recipe::PlotRecipe, page_key, view_key
) where {S <: AbstractPlotDefinition}
    axis_unit(S, dim, quantity, recipe)
end

"""
    axis_label(::Type{S}, dim, quantity, unit, recipe)

Return the displayed label for one axis.
"""
function axis_label(
        ::Type{S}, dim::Val, quantity::QuantityTag, unit::Units,
        recipe::PlotRecipe
) where {S <: AbstractPlotDefinition}
    quantity_label = get_label(quantity)
    unit_label = get_label(unit)
    return isempty(unit_label) ? quantity_label : "$quantity_label [$unit_label]"
end
function axis_label(
        ::Type{S}, mode::Val, dim::Val, quantity::QuantityTag, unit::Units,
        recipe::PlotRecipe, page_key, view_key
) where {S <: AbstractPlotDefinition}
    axis_label(S, dim, quantity, unit, recipe)
end

"""
    axis_scale(::Type{S}, dim, recipe)

Return the initial scale for one axis.
"""
function axis_scale(::Type{S}, dim::Val, recipe::PlotRecipe) where {S <:
                                                                    AbstractPlotDefinition}
    :linear
end
function axis_scale(
        ::Type{S}, mode::Val, dim::Val, recipe::PlotRecipe,
        page_key, view_key
) where {S <: AbstractPlotDefinition}
    axis_scale(S, dim, recipe)
end

"""
    axis_scales(::Type{S}, dim, recipe, series)

Return the scales available to one fully resolved axis.
"""
function axis_scales(
        ::Type{S}, dim::Val, recipe::PlotRecipe,
        series::Vector{SeriesSpec}
) where {S <: AbstractPlotDefinition}
    (axis_scale(S, dim, recipe),)
end
function axis_scales(
        ::Type{S}, mode::Val, dim::Val, recipe::PlotRecipe, page_key, view_key,
        series::Vector{SeriesSpec}
) where {S <: AbstractPlotDefinition}
    axis_scales(S, dim, recipe, series)
end

"""
    axis_exponent(::Type{S}, dim, recipe, series)

Return the base-ten display exponent for linear ticks on one axis.
"""
function axis_exponent(
        ::Type{S}, dim::Val, recipe::PlotRecipe,
        series::Vector{SeriesSpec}
) where {S <: AbstractPlotDefinition}
    0
end
function axis_exponent(
        ::Type{S}, mode::Val, dim::Val, recipe::PlotRecipe, page_key, view_key,
        series::Vector{SeriesSpec}
) where {S <: AbstractPlotDefinition}
    axis_exponent(S, dim, recipe, series)
end

"""
    axis_attributes(::Type{S}, dim, recipe)

Return visual renderer attributes for one axis.
"""
function axis_attributes(
        ::Type{S}, dim::Val, recipe::PlotRecipe
) where {S <: AbstractPlotDefinition}
    (;)
end
function axis_attributes(
        ::Type{S}, mode::Val, dim::Val, recipe::PlotRecipe,
        page_key, view_key
) where {S <: AbstractPlotDefinition}
    axis_attributes(S, dim, recipe)
end

function _make_axis(
        ::Type{S}, mode::Val, ::Val{dim}, recipe::PlotRecipe, page_key,
        view_key
) where {S <: AbstractPlotDefinition, dim}
    quantity = axis_quantity(S, mode, Val(dim), recipe, page_key, view_key)
    unit = axis_unit(S, mode, Val(dim), quantity, recipe, page_key, view_key)
    label = axis_label(S, mode, Val(dim), quantity, unit, recipe, page_key, view_key)
    scale = axis_scale(S, mode, Val(dim), recipe, page_key, view_key)
    return AxisSpec(
        dim,
        quantity,
        unit,
        label,
        scale;
        attributes = axis_attributes(S, mode, Val(dim), recipe, page_key, view_key)
    )
end

"""
    make_axes(::Type{S}, mode, recipe, page_key, view_key)

Construct backend-neutral axes for one view.
"""
function make_axes(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key,
        view_key
) where {S <: AbstractPlotDefinition}
    dims = geom_axes(S, mode, recipe, page_key, view_key)
    all(dim -> dim in (:x, :y, :z), dims) || throw(
        ArgumentError("geom_axes($S) may only contain :x, :y, and :z"),
    )
    length(unique(dims)) == length(dims) || throw(
        ArgumentError("geom_axes($S) cannot contain duplicate dimensions"),
    )
    xaxis = :x in dims ? _make_axis(S, mode, Val(:x), recipe, page_key, view_key) : nothing
    yaxis = :y in dims ? _make_axis(S, mode, Val(:y), recipe, page_key, view_key) : nothing
    zaxis = :z in dims ? _make_axis(S, mode, Val(:z), recipe, page_key, view_key) : nothing
    return (; xaxis, yaxis, zaxis)
end

"""
    plot_kind(::Type{S}, recipe, series_key)

Return the primitive symbol used by one semantic series facet.
"""
function plot_kind(::Type{S}, recipe::PlotRecipe, series_key) where {
        S <: AbstractPlotDefinition,
}
    :line
end
function plot_kind(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key, view_key,
        series_key
) where {S <: AbstractPlotDefinition}
    plot_kind(S, recipe, series_key)
end

"""
    series_data(::Type{S}, dim, recipe, series_key)

Return data for axis `dim` and one semantic series facet.
"""
function series_data(::Type{S}, dim::Val, recipe::PlotRecipe, series_key) where {
        S <: AbstractPlotDefinition,
}
    nothing
end
function series_data(
        ::Type{S}, mode::Val, dim::Val, recipe::PlotRecipe,
        page_key, view_key, series_key
) where {S <: AbstractPlotDefinition}
    series_data(S, dim, recipe, series_key)
end

"""
    legend_label(::Type{S}, recipe, series_key)

Return the legend label for one semantic series facet.
"""
function legend_label(::Type{S}, recipe::PlotRecipe, series_key) where {S <:
                                                                        AbstractPlotDefinition}
    nothing
end
function legend_label(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key, view_key,
        series_key
) where {S <: AbstractPlotDefinition}
    legend_label(S, recipe, series_key)
end

"""
    series_group(::Type{S}, recipe, series_key)

Return the visibility-group symbol for one semantic series facet.
"""
function series_group(::Type{S}, recipe::PlotRecipe, series_key) where {S <:
                                                                        AbstractPlotDefinition}
    nothing
end
function series_group(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key, view_key,
        series_key
) where {S <: AbstractPlotDefinition}
    series_group(S, recipe, series_key)
end

"""
    series_visible(::Type{S}, recipe, series_key)

Return the initial visibility of one semantic series facet.
"""
function series_visible(::Type{S}, recipe::PlotRecipe, series_key) where {S <:
                                                                          AbstractPlotDefinition}
    true
end
function series_visible(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key, view_key,
        series_key
) where {S <: AbstractPlotDefinition}
    series_visible(S, recipe, series_key)
end

"""
    series_attributes(::Type{S}, recipe, series_key)

Return backend-neutral visual attributes for one series facet.
"""
function series_attributes(::Type{S}, recipe::PlotRecipe, series_key) where {
        S <: AbstractPlotDefinition,
}
    (;)
end
function series_attributes(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key, view_key,
        series_key
) where {S <: AbstractPlotDefinition}
    series_attributes(S, recipe, series_key)
end

"""
    make_series(::Type{S}, mode, grouping, recipe, page_key, view_key, axes)

Construct backend-neutral primitive specifications for one view.
"""
function make_series(
        ::Type{S}, mode::Val, grouping::Val, recipe::PlotRecipe,
        page_key, view_key, axes::NamedTuple
) where {S <: AbstractPlotDefinition}
    series = SeriesSpec[]
    for series_key in series_keys(S, mode, grouping, recipe, page_key, view_key)
        push!(
            series,
            SeriesSpec(
                plot_kind(S, mode, recipe, page_key, view_key, series_key),
                series_data(S, mode, Val(:x), recipe, page_key, view_key, series_key),
                series_data(S, mode, Val(:y), recipe, page_key, view_key, series_key),
                series_data(S, mode, Val(:z), recipe, page_key, view_key, series_key),
                legend_label(S, mode, recipe, page_key, view_key, series_key);
                group = series_group(S, mode, recipe, page_key, view_key, series_key),
                visible = series_visible(S, mode, recipe, page_key, view_key, series_key),
                attributes = series_attributes(
                    S, mode, recipe, page_key, view_key, series_key)
            )
        )
    end
    return series
end

function _decorate_axis(
        axis::Nothing, ::Type{S}, mode::Val, dim::Val, recipe::PlotRecipe,
        page_key, view_key, series::Vector{SeriesSpec}
) where {S <: AbstractPlotDefinition}
    nothing
end
function _decorate_axis(
        axis::AxisSpec, ::Type{S}, mode::Val, dim::Val, recipe::PlotRecipe,
        page_key, view_key, series::Vector{SeriesSpec}
) where {S <: AbstractPlotDefinition}
    scales = axis_scales(S, mode, dim, recipe, page_key, view_key, series)
    exponent = axis_exponent(S, mode, dim, recipe, page_key, view_key, series)
    return AxisSpec(
        axis.dim,
        axis.quantity,
        axis.units,
        axis.label,
        axis.scale;
        allowed_scales = scales,
        exponent,
        attributes = axis.attributes
    )
end

"""
    default_title(::Type{S}, recipe)

Return the title for a semantic page or view facet.
"""
default_title(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition} = ""
function default_title(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key,
        view_key
) where {S <: AbstractPlotDefinition}
    default_title(S, recipe)
end

"""
    view_key(::Type{S}, recipe)

Return the semantic identity for one view.
"""
view_key(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition} = (;)
function view_key(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key, key
) where {S <: AbstractPlotDefinition}
    view_key(S, recipe)
end

"""
    view_placement(::Type{S}, recipe)

Return the named-slot placement for one view.
"""
function view_placement(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition}
    PlacementSpec()
end
function view_placement(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key, key
) where {S <: AbstractPlotDefinition}
    view_placement(S, recipe)
end

"""
    view_aspect(::Type{S}, recipe)

Return the aspect declaration for one view.
"""
view_aspect(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition} = nothing
function view_aspect(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key, key
) where {S <: AbstractPlotDefinition}
    view_aspect(S, recipe)
end

"""
    view_limits(::Type{S}, recipe)

Return explicit axis limits for one view, or `nothing`.
"""
view_limits(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition} = nothing
function view_limits(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key, key
) where {S <: AbstractPlotDefinition}
    view_limits(S, recipe)
end

"""
    view_attributes(::Type{S}, recipe)

Return visual renderer attributes for one view.
"""
view_attributes(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition} = (;)
function view_attributes(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key, key
) where {S <: AbstractPlotDefinition}
    view_attributes(S, recipe)
end

"""
    make_views(::Type{S}, mode, grouping, recipe, page_key)

Construct backend-neutral views for one semantic page facet.
"""
function make_views(
        ::Type{S}, mode::Val, grouping::Val, recipe::PlotRecipe,
        page_key
) where {S <: AbstractPlotDefinition}
    views = ViewSpec[]
    for key in view_keys(S, mode, grouping, recipe, page_key)
        axes = make_axes(S, mode, recipe, page_key, key)
        series = make_series(S, mode, grouping, recipe, page_key, key, axes)
        xaxis = _decorate_axis(axes.xaxis, S, mode, Val(:x), recipe, page_key, key, series)
        yaxis = _decorate_axis(axes.yaxis, S, mode, Val(:y), recipe, page_key, key, series)
        zaxis = _decorate_axis(axes.zaxis, S, mode, Val(:z), recipe, page_key, key, series)
        push!(
            views,
            ViewSpec(
                xaxis,
                yaxis,
                zaxis,
                default_title(S, mode, recipe, page_key, key),
                series,
                view_key(S, mode, recipe, page_key, key);
                placement = view_placement(S, mode, recipe, page_key, key),
                aspect = view_aspect(S, mode, recipe, page_key, key),
                limits = view_limits(S, mode, recipe, page_key, key),
                attributes = view_attributes(S, mode, recipe, page_key, key)
            )
        )
    end
    return views
end

"""
    default_figsize(::Type{S}, recipe)

Return the default page width and height in pixels for a recipe.
"""
function default_figsize(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition}
    (800, 400)
end
function default_figsize(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key
) where {S <: AbstractPlotDefinition}
    default_figsize(S, recipe)
end

function _standard_layout(name::Symbol)
    root = GridSpec(
        :root;
        rows = AbstractTrackSize[FixedTrack(36), RelativeTrack(), FixedTrack(20)],
        columns = AbstractTrackSize[RelativeTrack(), ContentTrack()],
        rowgap = 6,
        columngap = 12,
        padding = (20, 20, 28, 28)
    )
    side = GridSpec(
        :side;
        parent = :root,
        area = GridArea(2, 2),
        rows = AbstractTrackSize[RelativeTrack(), ContentTrack()],
        columns = AbstractTrackSize[ContentTrack()],
        rowgap = 4
    )
    slots = [
        SlotSpec(:toolbar, :root, GridArea(1, 1:2); halign = :left, valign = :bottom),
        SlotSpec(:canvas, :root, GridArea(2, 1); halign = :stretch, valign = :stretch),
        SlotSpec(:status, :root, GridArea(3, 1:2); halign = :left, valign = :center),
        SlotSpec(:legend, :side, GridArea(1, 1); halign = :left, valign = :top),
        SlotSpec(:colorbars, :side, GridArea(2, 1); halign = :left, valign = :top)
    ]
    return LayoutSpec(name, [root, side], slots)
end

"""
    layout_preset(::Val{name}, view_count)

Construct a built-in named layout preset for `view_count` views.

# Errors

- `ArgumentError` when `name` is not `:single`, `:grid`, `:preview`, or
  `:material_scale`.
"""
layout_preset(::Val{:single}, view_count::Integer) = _standard_layout(:single)
layout_preset(::Val{:grid}, view_count::Integer) = _standard_layout(:grid)
layout_preset(::Val{:preview}, view_count::Integer) = _standard_layout(:preview)
function layout_preset(::Val{:material_scale}, view_count::Integer)
    root = GridSpec(
        :root;
        rows = AbstractTrackSize[FixedTrack(36), RelativeTrack(), FixedTrack(20)],
        columns = AbstractTrackSize[ContentTrack()],
        rowgap = 6,
        padding = (20, 20, 28, 28)
    )
    side = GridSpec(
        :side;
        parent = :root,
        area = GridArea(2, 1),
        rows = AbstractTrackSize[RelativeTrack()],
        columns = AbstractTrackSize[ContentTrack()]
    )
    slots = [
        SlotSpec(:toolbar, :root, GridArea(1, 1); halign = :left, valign = :bottom),
        SlotSpec(:status, :root, GridArea(3, 1); halign = :left, valign = :center),
        SlotSpec(:colorbars, :side, GridArea(1, 1); halign = :left, valign = :center)
    ]
    return LayoutSpec(:material_scale, [root, side], slots)
end
function layout_preset(::Val{name}, view_count::Integer) where {name}
    throw(ArgumentError("unknown PlotBuilder layout preset :$name"))
end

"""
    layout_spec(::Type{S}, recipe)

Return a layout preset symbol or complete `LayoutSpec` for a recipe page.
"""
layout_spec(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition} = :single
function layout_spec(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key
) where {S <: AbstractPlotDefinition}
    layout_spec(S, recipe)
end

function _resolve_layout(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key,
        view_count::Integer
) where {S <: AbstractPlotDefinition}
    selected = recipe.renderer.layout
    selected === nothing && (selected = layout_spec(S, mode, recipe, page_key))
    selected isa LayoutSpec && return selected
    selected isa Symbol || throw(
        ArgumentError("layout_spec($S) must return a preset symbol or LayoutSpec"),
    )
    return layout_preset(Val(selected), view_count)
end

"""
    page_identity(::Type{S}, recipe, page_key)

Return the semantic `NamedTuple` identity for one page.
"""
function page_identity(
        ::Type{S}, recipe::PlotRecipe, page_key
) where {S <: AbstractPlotDefinition}
    page_key === nothing && return (;)
    page_key isa NamedTuple && return page_key
    return (; facet = page_key)
end
function page_identity(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key
) where {S <: AbstractPlotDefinition}
    page_identity(S, recipe, page_key)
end

"""
    control_spec(::Type{S}, recipe)

Return the typed interactive-control declaration for a recipe page.
"""
function control_spec(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition}
    ControlSpec()
end
function control_spec(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key
) where {S <: AbstractPlotDefinition}
    control_spec(S, recipe)
end

"""
    legend_spec(::Type{S}, recipe)

Return the typed legend declaration for a recipe page.
"""
function legend_spec(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition}
    LegendSpec()
end
function legend_spec(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key
) where {S <: AbstractPlotDefinition}
    legend_spec(S, recipe)
end

"""
    colorbar_specs(::Type{S}, recipe)

Return typed colorbar declarations for a recipe page.
"""
function colorbar_specs(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition}
    ColorbarSpec[]
end
function colorbar_specs(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key
) where {S <: AbstractPlotDefinition}
    colorbar_specs(S, recipe)
end

"""
    status_spec(::Type{S}, recipe)

Return the typed status-line declaration for a recipe page.
"""
function status_spec(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition}
    StatusSpec()
end
function status_spec(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key
) where {S <: AbstractPlotDefinition}
    status_spec(S, recipe)
end

"""
    export_spec(::Type{S}, recipe, title)

Return the typed SVG export declaration for a recipe page.
"""
function export_spec(
        ::Type{S}, recipe::PlotRecipe, title::AbstractString
) where {S <: AbstractPlotDefinition}
    ExportSpec(
        theme = recipe.renderer.export_theme,
        name = title,
        open_file = recipe.renderer.open_export
    )
end
function export_spec(
        ::Type{S}, mode::Val, recipe::PlotRecipe, page_key,
        title::AbstractString
) where {S <: AbstractPlotDefinition}
    export_spec(S, recipe, title)
end

"""
    make_pages(::Type{S}, mode, grouping, recipe)

Construct every backend-neutral page for a resolved recipe.
"""
function make_pages(
        ::Type{S}, mode::Val, grouping::Val, recipe::PlotRecipe
) where {S <: AbstractPlotDefinition}
    pages = PageSpec[]
    for page_key in page_keys(S, mode, grouping, recipe)
        views = make_views(S, mode, grouping, recipe, page_key)
        title = default_title(S, mode, recipe, page_key, nothing)
        push!(
            pages,
            PageSpec(
                title,
                default_figsize(S, mode, recipe, page_key),
                page_identity(S, mode, recipe, page_key),
                _resolve_layout(S, mode, recipe, page_key, length(views)),
                views;
                controls = control_spec(S, mode, recipe, page_key),
                legend = legend_spec(S, mode, recipe, page_key),
                colorbars = colorbar_specs(S, mode, recipe, page_key),
                status = status_spec(S, mode, recipe, page_key),
                export_spec = export_spec(S, mode, recipe, page_key, title)
            )
        )
    end
    return pages
end

"""
    make_render(::Type{S}, object; kwargs...)

Materialize a domain object through the uniform PlotBuilder grammar. Plot
specifications specialize accessors; they do not replace this rendering
sequence.
"""
function make_render(::Type{S}, object; kwargs...) where {S <: AbstractPlotDefinition}
    expected = dispatch_on(S)
    object isa expected || throw(
        ArgumentError("$S accepts $expected, not $(typeof(object))"),
    )
    recipe = resolve_input(S, parse_kwargs(S, object; kwargs...))
    recipe isa PlotRecipe || throw(
        ArgumentError("resolve_input($S) must return PlotRecipe"),
    )
    mode = recipe_mode(S, recipe)
    mode isa Val || throw(ArgumentError("recipe_mode($S) must return Val(mode)"))
    grouping = grouping_mode(S, mode, recipe)
    grouping isa Val || throw(
        ArgumentError("grouping_mode($S) must return Val(mode)"),
    )
    return RenderSpec(S, make_pages(S, mode, grouping, recipe))
end
