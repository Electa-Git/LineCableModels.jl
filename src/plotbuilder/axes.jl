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
    axis_payload(::Type{S}, dim, recipe)

Return the published or synthetic payload used to define one axis.
"""
function axis_payload(::Type{S}, dim::Val, recipe::PlotRecipe) where {S <:
                                                                      AbstractPlotDefinition}
    return (
        values = nothing,
        quantity = QuantityTag{:dimensionless}(),
        unit = units(:base, :dimensionless)
    )
end
function axis_payload(
        ::Type{S}, mode::Val, dim::Val, recipe::PlotRecipe,
        page_key, view_key
) where {S <: AbstractPlotDefinition}
    axis_payload(S, dim, recipe)
end

"""
    axis_label(::Type{S}, dim, quantity, unit, recipe)

Return the displayed label for one axis.
"""
function axis_label(
        ::Type{S}, dim::Val, quantity::QuantityTag, unit::UnitExpr,
        recipe::PlotRecipe
) where {S <: AbstractPlotDefinition}
    quantity_label = label(quantity)
    unit_label = label(unit)
    return isempty(unit_label) ? quantity_label : "$quantity_label [$unit_label]"
end
function axis_label(
        ::Type{S}, mode::Val, dim::Val, quantity::QuantityTag, unit::UnitExpr,
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
    payload = axis_payload(S, mode, Val(dim), recipe, page_key, view_key)
    keys(payload) == (:values, :quantity, :unit) || throw(
        ArgumentError("axis payloads must contain values, quantity, and unit"),
    )
    quantity = payload.quantity
    unit = payload.unit
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

Construct renderer-independent axes for one view.
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

"Construct the axes for every page/view facet in one stage."
function make_axes(
        ::Type{S}, mode::Val, grouping::Val, recipe::PlotRecipe
) where {S <: AbstractPlotDefinition}
    entries = NamedTuple[]
    for page_key in _page_keys(S, mode, grouping, recipe)
        for view_key in _view_keys(S, mode, grouping, recipe, page_key)
            push!(entries, (;
                page_key,
                view_key,
                axes = make_axes(S, mode, recipe, page_key, view_key)
            ))
        end
    end
    return entries
end
