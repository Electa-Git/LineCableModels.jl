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
    series_values(::Type{S}, dim, recipe, series_key)

Return data for axis `dim` and one semantic series facet.
"""
function series_values(::Type{S}, dim::Val, recipe::PlotRecipe, series_key) where {
        S <: AbstractPlotDefinition,
}
    nothing
end
function series_values(
        ::Type{S}, mode::Val, dim::Val, recipe::PlotRecipe,
        page_key, view_key, series_key
) where {S <: AbstractPlotDefinition}
    series_values(S, dim, recipe, series_key)
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

Return renderer-independent visual attributes for one series facet.
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

Construct renderer-independent primitive declarations for one view.
"""
function make_series(
        ::Type{S}, mode::Val, grouping::Val, recipe::PlotRecipe,
        page_key, view_key, axes::NamedTuple
) where {S <: AbstractPlotDefinition}
    series = SeriesSpec[]
    for series_key in _series_keys(S, mode, grouping, recipe, page_key, view_key)
        push!(
            series,
            SeriesSpec(
                plot_kind(S, mode, recipe, page_key, view_key, series_key),
                series_values(S, mode, Val(:x), recipe, page_key, view_key, series_key),
                series_values(S, mode, Val(:y), recipe, page_key, view_key, series_key),
                series_values(S, mode, Val(:z), recipe, page_key, view_key, series_key),
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

"Construct the series for every axis bundle in one stage."
function make_series(
        ::Type{S}, mode::Val, grouping::Val, recipe::PlotRecipe,
        axis_entries::AbstractVector
) where {S <: AbstractPlotDefinition}
    return map(axis_entries) do entry
        merge(entry, (;
            series = make_series(
                S,
                mode,
                grouping,
                recipe,
                entry.page_key,
                entry.view_key,
                entry.axes
            ),
        ))
    end
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
