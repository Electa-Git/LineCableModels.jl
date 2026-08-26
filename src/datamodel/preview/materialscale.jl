PlotBuilder.entitle(::Type{MaterialScalePlotDefinition}, value::Nothing) = value
function PlotBuilder.renderer_defaults(::Type{MaterialScalePlotDefinition}, ::Nothing)
    return (; size = (800, 400))
end

function PlotBuilder.resolve(
        ::Type{MaterialScalePlotDefinition},
        ::Nothing,
        request::NamedTuple
)
    request.renderer.size isa Tuple{Int, Int} || throw(
        ArgumentError("size must be a tuple of two integers"),
    )
    all(>(0), request.renderer.size) || throw(
        ArgumentError("size dimensions must be positive"),
    )
    return request
end

function PlotBuilder.fetch(
        ::Type{MaterialScalePlotDefinition},
        ::Nothing,
        request::NamedTuple
)
    title = "Material property colour scale"
    colorbars = _colorbar_specs(
        (_RHO_MIN, _RHO_MAX),
        (1.0, 300.0),
        (1.0, 1000.0)
    )
    return PlotBuilder.PlotPage[
        PlotBuilder.PlotPage(
        title,
        request.renderer.size,
        (; kind = :material_scale),
        nothing;
        legend = PlotBuilder.LegendDefinition(enabled = false),
        colorbars,
        export_definition = PlotBuilder.ExportDefinition(
            theme = request.renderer.export_theme,
            name = "material_scale",
            open_file = request.renderer.open_export
        )
    ),
    ]
end
