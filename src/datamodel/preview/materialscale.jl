PlotBuilder.dispatch_on(::Type{MaterialScalePlotDefinition}) = Nothing
PlotBuilder.renderer_kwargs(::Type{MaterialScalePlotDefinition}) = (:size,)
function PlotBuilder.renderer_defaults(::Type{MaterialScalePlotDefinition}, ::Nothing)
    (; size = (800, 400))
end

function PlotBuilder.resolve(::Type{MaterialScalePlotDefinition}, recipe::PlotBuilder.PlotRecipe)
    recipe.renderer.size isa Tuple{Int, Int} || throw(
        ArgumentError("size must be a tuple of two integers"),
    )
    return recipe
end

function PlotBuilder._composition(
        ::Type{MaterialScalePlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe
)
    return Val(:empty)
end

function PlotBuilder.default_title(
        ::Type{MaterialScalePlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return "Material property colour scale"
end

function PlotBuilder.default_figsize(
        ::Type{MaterialScalePlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return recipe.renderer.size
end

function PlotBuilder.layout_spec(
        ::Type{MaterialScalePlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return :material_scale
end

function PlotBuilder.colorbar_specs(
        ::Type{MaterialScalePlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return _colorbar_specs(
        (_RHO_MIN, _RHO_MAX),
        (1.0, 300.0),
        (1.0, 1000.0)
    )
end

function PlotBuilder.control_spec(
        ::Type{MaterialScalePlotDefinition}, mode::Val,
        recipe::PlotBuilder.PlotRecipe, page_key)
    return PlotBuilder.ControlSpec(reset = false)
end

function PlotBuilder.legend_spec(
        ::Type{MaterialScalePlotDefinition}, mode::Val,
        recipe::PlotBuilder.PlotRecipe, page_key)
    return PlotBuilder.LegendSpec(enabled = false)
end

function PlotBuilder.page_identity(
        ::Type{MaterialScalePlotDefinition}, mode::Val,
        recipe::PlotBuilder.PlotRecipe, page_key)
    return (; kind = :material_scale)
end

function PlotBuilder.export_spec(
        ::Type{MaterialScalePlotDefinition}, mode::Val,
        recipe::PlotBuilder.PlotRecipe, page_key, title::AbstractString)
    return PlotBuilder.ExportSpec(
        theme = recipe.renderer.export_theme,
        name = "material_scale",
        open_file = recipe.renderer.open_export
    )
end
