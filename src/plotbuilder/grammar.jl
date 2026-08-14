const COMMON_RENDERER_KWARGS = (:export_theme, :open_export)

"""Return the domain type accepted by plot specification `S`."""
dispatch_on(::Type{S}) where {S <: AbstractPlotSpec} = Any

"""Return semantic keyword names accepted by plot specification `S`."""
input_kwargs(::Type{S}) where {S <: AbstractPlotSpec} = ()

"""Return renderer keyword names accepted by plot specification `S`."""
renderer_kwargs(::Type{S}) where {S <: AbstractPlotSpec} = ()

"""Return semantic keyword defaults for `object`."""
input_defaults(::Type{S}, object) where {S <: AbstractPlotSpec} = (;)

"""Return renderer keyword defaults for `object`."""
renderer_defaults(::Type{S}, object) where {S <: AbstractPlotSpec} = (;)

function _symbol_tuple(::Type{S}, values, accessor::Symbol) where {S <: AbstractPlotSpec}
    values isa Tuple || throw(
        ArgumentError("$accessor($S) must return a tuple of symbols"),
    )
    all(value -> value isa Symbol, values) || throw(
        ArgumentError("$accessor($S) must return a tuple of symbols"),
    )
    return values
end

function _select_kwargs(kwargs::NamedTuple, names::Tuple)
    selected = Tuple(pair for pair in pairs(kwargs) if first(pair) in names)
    return NamedTuple(selected)
end

"""
    parse_kwargs(::Type{S}, object, kwargs)

Split user input into semantic and renderer options declared by the recipe.
Unsupported keywords are reported and ignored, matching the established plot
entry-point behavior.
"""
function parse_kwargs(
        ::Type{S},
        object,
        kwargs::NamedTuple
) where {S <: AbstractPlotSpec}
    input_names = _symbol_tuple(S, input_kwargs(S), :input_kwargs)
    renderer_names = _symbol_tuple(S, renderer_kwargs(S), :renderer_kwargs)
    allowed_renderer = (COMMON_RENDERER_KWARGS..., renderer_names...)
    allowed = (input_names..., allowed_renderer...)
    for name in keys(kwargs)
        name in allowed || @warn "unsupported plot keyword" specification=S keyword=name
    end

    input = merge(input_defaults(S, object), _select_kwargs(kwargs, input_names))
    renderer = merge(
        (; export_theme = :default, open_export = true),
        renderer_defaults(S, object),
        _select_kwargs(kwargs, allowed_renderer)
    )
    _validate_export_theme(renderer.export_theme)
    renderer.open_export isa Bool || throw(ArgumentError("open_export must be Bool"))
    return (; object, input, renderer)
end

function parse_kwargs(::Type{S}, object; kwargs...) where {S <: AbstractPlotSpec}
    parse_kwargs(S, object, (; kwargs...))
end

"""Validate and enrich parsed recipe input before materialization."""
resolve_input(::Type{S}, recipe::NamedTuple) where {S <: AbstractPlotSpec} = recipe

"""Return the value-dispatched plotting mode for a resolved recipe."""
recipe_mode(::Type{S}, recipe::NamedTuple) where {S <: AbstractPlotSpec} = Val(:default)

"""Return the value-dispatched grouping mode for a resolved recipe."""
function grouping_mode(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple
) where {S <: AbstractPlotSpec}
    Val(:overlay)
end

"""Return semantic page facets for recipes using `:faceted_pages`."""
function page_facets(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple
) where {S <: AbstractPlotSpec}
    (nothing,)
end

"""Return semantic series facets for the selected recipe page."""
function group_facets(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key
) where {S <: AbstractPlotSpec}
    (nothing,)
end

function page_keys(::Type{S}, mode::Val, ::Val{:overlay}, recipe::NamedTuple) where {
        S <: AbstractPlotSpec,
}
    (nothing,)
end
function page_keys(::Type{S}, mode::Val, ::Val{:panels}, recipe::NamedTuple) where {
        S <: AbstractPlotSpec,
}
    (nothing,)
end
function page_keys(::Type{S}, mode::Val, ::Val{:pages}, recipe::NamedTuple) where {
        S <: AbstractPlotSpec,
}
    group_facets(S, mode, recipe, nothing)
end
function page_keys(::Type{S}, mode::Val, ::Val{:faceted_pages},
        recipe::NamedTuple) where {
        S <: AbstractPlotSpec,
}
    page_facets(S, mode, recipe)
end
function page_keys(::Type{S}, mode::Val, ::Val{:empty}, recipe::NamedTuple) where {
        S <: AbstractPlotSpec,
}
    (nothing,)
end
function page_keys(
        ::Type{S},
        mode::Val,
        grouping::Val,
        recipe::NamedTuple
) where {S <: AbstractPlotSpec}
    grouping_name = only(typeof(grouping).parameters)
    throw(
        ArgumentError(
        "unsupported grouping mode :$grouping_name for $S; " *
        "specialize PlotBuilder.page_keys for this Val mode"
    ),
    )
end

function view_keys(
        ::Type{S},
        mode::Val,
        ::Val{:overlay},
        recipe::NamedTuple,
        page_key
) where {S <: AbstractPlotSpec}
    (nothing,)
end
function view_keys(
        ::Type{S},
        mode::Val,
        ::Val{:panels},
        recipe::NamedTuple,
        page_key
) where {S <: AbstractPlotSpec}
    group_facets(S, mode, recipe, page_key)
end
function view_keys(
        ::Type{S},
        mode::Val,
        ::Val{:pages},
        recipe::NamedTuple,
        page_key
) where {S <: AbstractPlotSpec}
    (nothing,)
end
function view_keys(
        ::Type{S},
        mode::Val,
        ::Val{:faceted_pages},
        recipe::NamedTuple,
        page_key
) where {S <: AbstractPlotSpec}
    (nothing,)
end
function view_keys(
        ::Type{S},
        mode::Val,
        ::Val{:empty},
        recipe::NamedTuple,
        page_key
) where {S <: AbstractPlotSpec}
    ()
end

function series_keys(
        ::Type{S},
        mode::Val,
        ::Val{:overlay},
        recipe::NamedTuple,
        page_key,
        view_key
) where {S <: AbstractPlotSpec}
    group_facets(S, mode, recipe, page_key)
end
function series_keys(
        ::Type{S},
        mode::Val,
        ::Val{:panels},
        recipe::NamedTuple,
        page_key,
        view_key
) where {S <: AbstractPlotSpec}
    (view_key,)
end
function series_keys(
        ::Type{S},
        mode::Val,
        ::Val{:pages},
        recipe::NamedTuple,
        page_key,
        view_key
) where {S <: AbstractPlotSpec}
    (page_key,)
end
function series_keys(
        ::Type{S},
        mode::Val,
        ::Val{:faceted_pages},
        recipe::NamedTuple,
        page_key,
        view_key
) where {S <: AbstractPlotSpec}
    group_facets(S, mode, recipe, page_key)
end

"""Return geometric axes used by a recipe view."""
function geom_axes(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key
) where {S <: AbstractPlotSpec}
    (:x, :y)
end

function axis_quantity(
        ::Type{S},
        mode::Val,
        ::Val{dim},
        recipe::NamedTuple,
        page_key,
        view_key
) where {S <: AbstractPlotSpec, dim}
    return QuantityTag{:unknown}()
end

function axis_unit(
        ::Type{S},
        mode::Val,
        ::Val{dim},
        quantity::QuantityTag,
        recipe::NamedTuple,
        page_key,
        view_key
) where {S <: AbstractPlotSpec, dim}
    return display_unit(quantity)
end

function axis_label(
        ::Type{S},
        mode::Val,
        ::Val{dim},
        quantity::QuantityTag,
        unit::Units,
        recipe::NamedTuple,
        page_key,
        view_key
) where {S <: AbstractPlotSpec, dim}
    quantity_label = get_label(quantity)
    unit_label = get_label(unit)
    return isempty(unit_label) ? quantity_label : "$quantity_label [$unit_label]"
end

function axis_scale(
        ::Type{S},
        mode::Val,
        ::Val{dim},
        recipe::NamedTuple,
        page_key,
        view_key
) where {S <: AbstractPlotSpec, dim}
    return :linear
end

function _make_axis(
        ::Type{S},
        mode::Val,
        ::Val{dim},
        recipe::NamedTuple,
        page_key,
        view_key
) where {S <: AbstractPlotSpec, dim}
    quantity = axis_quantity(S, mode, Val(dim), recipe, page_key, view_key)
    unit = axis_unit(S, mode, Val(dim), quantity, recipe, page_key, view_key)
    label = axis_label(S, mode, Val(dim), quantity, unit, recipe, page_key, view_key)
    scale = axis_scale(S, mode, Val(dim), recipe, page_key, view_key)
    return AxisSpec(dim, quantity, unit, label, scale)
end

function make_axes(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key
) where {S <: AbstractPlotSpec}
    dims = geom_axes(S, mode, recipe, page_key, view_key)
    all(dim -> dim in (:x, :y, :z), dims) || throw(
        ArgumentError("geom_axes($S) may only contain :x, :y, and :z"),
    )
    xaxis = :x in dims ? _make_axis(S, mode, Val(:x), recipe, page_key, view_key) : nothing
    yaxis = :y in dims ? _make_axis(S, mode, Val(:y), recipe, page_key, view_key) : nothing
    zaxis = :z in dims ? _make_axis(S, mode, Val(:z), recipe, page_key, view_key) : nothing
    return (; xaxis, yaxis, zaxis)
end

"""Return the primitive used by one semantic series facet."""
function plot_kind(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key,
        series_key
) where {S <: AbstractPlotSpec}
    :line
end

"""Return data for axis `dim` and one semantic series facet."""
function series_data(
        ::Type{S},
        mode::Val,
        ::Val{dim},
        recipe::NamedTuple,
        page_key,
        view_key,
        series_key
) where {S <: AbstractPlotSpec, dim}
    return nothing
end

"""Return the legend label for one semantic series facet."""
function legend_label(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key,
        series_key
) where {S <: AbstractPlotSpec}
    nothing
end

"""Return backend-neutral primitive attributes for one series facet."""
function series_attributes(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key,
        series_key
) where {S <: AbstractPlotSpec}
    (;)
end

function make_series(
        ::Type{S},
        mode::Val,
        grouping::Val,
        recipe::NamedTuple,
        page_key,
        view_key,
        axes::NamedTuple
) where {S <: AbstractPlotSpec}
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
                attributes = series_attributes(
                    S,
                    mode,
                    recipe,
                    page_key,
                    view_key,
                    series_key
                )
            )
        )
    end
    return series
end

"""Return the title for a semantic page or view facet."""
function default_title(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key
) where {S <: AbstractPlotSpec}
    ""
end

function view_key(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        key
) where {S <: AbstractPlotSpec}
    (;)
end

function view_attributes(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        key
) where {S <: AbstractPlotSpec}
    (;)
end

function make_views(
        ::Type{S},
        mode::Val,
        grouping::Val,
        recipe::NamedTuple,
        page_key
) where {S <: AbstractPlotSpec}
    views = ViewSpec[]
    for key in view_keys(S, mode, grouping, recipe, page_key)
        axes = make_axes(S, mode, recipe, page_key, key)
        series = make_series(S, mode, grouping, recipe, page_key, key, axes)
        push!(
            views,
            ViewSpec(
                axes.xaxis,
                axes.yaxis,
                axes.zaxis,
                default_title(S, mode, recipe, page_key, key),
                series,
                view_key(S, mode, recipe, page_key, key);
                attributes = view_attributes(S, mode, recipe, page_key, key)
            )
        )
    end
    return views
end

"""Return the default page size for a recipe."""
function default_figsize(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key
) where {S <: AbstractPlotSpec}
    (800, 400)
end

"""Return the renderer layout symbol for a recipe page."""
function figure_layout(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key
) where {S <: AbstractPlotSpec}
    :single
end

"""Return axes whose logarithmic controls are available on a recipe page."""
function enable_logscale(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key
) where {S <: AbstractPlotSpec}
    ()
end

"""Return recipe-specific page attributes."""
function page_kwargs(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        views::Vector{ViewSpec}
) where {S <: AbstractPlotSpec}
    (;)
end

function _common_page_kwargs(
        ::Type{S},
        mode::Val,
        recipe::NamedTuple,
        page_key
) where {S <: AbstractPlotSpec}
    log_axes = enable_logscale(S, mode, recipe, page_key)
    controls = control_definitions(xlog = :x in log_axes, ylog = :y in log_axes)
    return (;
        export_theme = recipe.renderer.export_theme,
        open_export = recipe.renderer.open_export,
        controls
    )
end

function make_pages(
        ::Type{S},
        mode::Val,
        grouping::Val,
        recipe::NamedTuple
) where {S <: AbstractPlotSpec}
    pages = PageSpec[]
    for page_key in page_keys(S, mode, grouping, recipe)
        views = make_views(S, mode, grouping, recipe, page_key)
        title = default_title(S, mode, recipe, page_key, nothing)
        kwargs = merge(
            _common_page_kwargs(S, mode, recipe, page_key),
            page_kwargs(S, mode, recipe, page_key, views)
        )
        push!(
            pages,
            PageSpec(
                title,
                default_figsize(S, mode, recipe, page_key),
                figure_layout(S, mode, recipe, page_key),
                views,
                kwargs
            )
        )
    end
    return pages
end

"""
    make_render(::Type{S}, object; kwargs...)

Materialize a domain object through the uniform PlotBuilder grammar. Plot
specifications specialize accessors; they do not replace this pipeline.
"""
function make_render(::Type{S}, object; kwargs...) where {S <: AbstractPlotSpec}
    expected = dispatch_on(S)
    object isa expected || throw(
        ArgumentError("$S accepts $expected, not $(typeof(object))"),
    )
    parsed = parse_kwargs(S, object; kwargs...)
    recipe = resolve_input(S, parsed)
    mode = recipe_mode(S, recipe)
    mode isa Val || throw(ArgumentError("recipe_mode($S) must return Val(mode)"))
    grouping = grouping_mode(S, mode, recipe)
    grouping isa Val || throw(
        ArgumentError("grouping_mode($S) must return Val(mode)"),
    )
    return RenderSpec(S, make_pages(S, mode, grouping, recipe))
end
