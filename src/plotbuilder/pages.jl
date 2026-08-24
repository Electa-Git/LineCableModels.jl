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

Construct every renderer-independent page for a resolved recipe.
"""
function make_pages(
        ::Type{S}, mode::Val, grouping::Val, recipe::PlotRecipe
) where {S <: AbstractPlotDefinition}
    pages = PageSpec[]
    for page_key in _page_keys(S, mode, grouping, recipe)
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
