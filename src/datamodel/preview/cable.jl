PlotBuilder.dispatch_on(::Type{CablePreviewPlotDefinition}) = CableDesign
function PlotBuilder.input_kwargs(::Type{CablePreviewPlotDefinition})
    (
        :x_offset,
        :y_offset,
        :display_legend,
        :display_id,
        :display_colorbars
    )
end
PlotBuilder.renderer_kwargs(::Type{CablePreviewPlotDefinition}) = (:size,)
function PlotBuilder.input_defaults(::Type{CablePreviewPlotDefinition}, ::CableDesign)
    (;
        x_offset = 0.0,
        y_offset = 0.0,
        display_legend = true,
        display_id = false,
        display_colorbars = true
    )
end
function PlotBuilder.renderer_defaults(::Type{CablePreviewPlotDefinition}, ::CableDesign)
    (; size = (900, 700))
end

function PlotBuilder.resolve(::Type{CablePreviewPlotDefinition}, recipe::PlotBuilder.PlotRecipe)
    recipe.input.x_offset isa Real || throw(ArgumentError("x_offset must be real"))
    recipe.input.y_offset isa Real || throw(ArgumentError("y_offset must be real"))
    all(name -> getproperty(recipe.input, name) isa Bool,
        (:display_legend, :display_id, :display_colorbars)) || throw(
        ArgumentError("display_legend, display_id, and display_colorbars must be Bool"),
    )
    recipe.renderer.size isa Tuple{Int, Int} || throw(
        ArgumentError("size must be a tuple of two integers"),
    )
    return recipe
end

function PlotBuilder.fetch(::Type{CablePreviewPlotDefinition}, recipe::PlotBuilder.PlotRecipe)
    design = recipe.object
    series = _design_series(
        design,
        recipe.input.x_offset,
        recipe.input.y_offset;
        display_legend = recipe.input.display_legend
    )
    title = _cable_title(Val(recipe.input.display_id), design)
    colorbars = _cable_colorbars(Val(recipe.input.display_colorbars), design)
    identity = (; kind = :cable, id = design.cable_id)
    export_name = design.cable_id
    return PlotBuilder.PlotRecipe(
        CablePreviewPlotDefinition,
        design,
        merge(recipe.input, (; series, title, colorbars, identity, export_name)),
        recipe.renderer
    )
end

function PlotBuilder.make_axes(
        ::Type{CablePreviewPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    xaxis, yaxis = _distance_axes()
    return (; xaxis, yaxis, zaxis = nothing)
end

function PlotBuilder.make_series(
        ::Type{CablePreviewPlotDefinition},
        mode::Val,
        grouping::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key,
        axes::NamedTuple
)
    return recipe.input.series
end

_cable_title(::Val{false}, design) = "Cable design preview"
_cable_title(::Val{true}, design) = "Cable design preview: $(design.cable_id)"

_cable_colorbars(::Val{false}, design) = PlotBuilder.ColorbarSpec[]
_cable_colorbars(::Val{true}, design) = _colorbar_specs(_property_ranges(design)...)

function PlotBuilder.default_title(
        ::Type{CablePreviewPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return recipe.input.title
end

function PlotBuilder.view_key(
        ::Type{CablePreviewPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return (; kind = :cable)
end

function PlotBuilder.view_aspect(
        ::Type{CablePreviewPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return :data
end

function PlotBuilder.default_figsize(
        ::Type{CablePreviewPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return recipe.renderer.size
end

function PlotBuilder.layout_spec(
        ::Type{CablePreviewPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return :preview
end

function PlotBuilder.colorbar_specs(
        ::Type{CablePreviewPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return recipe.input.colorbars
end

function PlotBuilder.legend_spec(
        ::Type{CablePreviewPlotDefinition}, mode::Val,
        recipe::PlotBuilder.PlotRecipe, page_key)
    return PlotBuilder.LegendSpec(enabled = recipe.input.display_legend)
end

function PlotBuilder.page_identity(
        ::Type{CablePreviewPlotDefinition}, mode::Val,
        recipe::PlotBuilder.PlotRecipe, page_key)
    return recipe.input.identity
end

function PlotBuilder.export_spec(
        ::Type{CablePreviewPlotDefinition}, mode::Val,
        recipe::PlotBuilder.PlotRecipe, page_key, title::AbstractString)
    return PlotBuilder.ExportSpec(
        theme = recipe.renderer.export_theme,
        name = recipe.input.export_name,
        open_file = recipe.renderer.open_export
    )
end
