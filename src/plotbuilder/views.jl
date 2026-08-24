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
    make_views(::Type{S}, mode, grouping, recipe, series_entries)

Construct renderer-independent views from the completed series stage.
"""
function make_views(
        ::Type{S}, mode::Val, grouping::Val, recipe::PlotRecipe,
        series_entries::AbstractVector
) where {S <: AbstractPlotDefinition}
    views = NamedTuple[]
    for entry in series_entries
        page_key = entry.page_key
        key = entry.view_key
        axes = entry.axes
        series = entry.series
        xaxis = _decorate_axis(
            axes.xaxis, S, mode, Val(:x), recipe, page_key, key, series)
        yaxis = _decorate_axis(
            axes.yaxis, S, mode, Val(:y), recipe, page_key, key, series)
        zaxis = _decorate_axis(
            axes.zaxis, S, mode, Val(:z), recipe, page_key, key, series)
        push!(
            views,
            (;
                page_key,
                view = ViewSpec(
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
        )
    end
    return views
end
