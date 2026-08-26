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

function PlotBuilder.resolve(
        ::Type{CablePreviewPlotDefinition},
        ::CableDesign,
        request::NamedTuple
)
    request.input.x_offset isa Real || throw(ArgumentError("x_offset must be real"))
    request.input.y_offset isa Real || throw(ArgumentError("y_offset must be real"))
    all(name -> getproperty(request.input, name) isa Bool,
        (:display_legend, :display_id, :display_colorbars)) || throw(
        ArgumentError("display_legend, display_id, and display_colorbars must be Bool"),
    )
    request.renderer.size isa Tuple{Int, Int} || throw(
        ArgumentError("size must be a tuple of two integers"),
    )
    return request
end

function PlotBuilder.fetch(
        ::Type{CablePreviewPlotDefinition},
        design::CableDesign,
        request::NamedTuple
)
    polygons = _design_shapes(
        design,
        request.input.x_offset,
        request.input.y_offset;
        display_legend = request.input.display_legend
    )
    title = _cable_title(Val(request.input.display_id), design)
    colorbars = _cable_colorbars(Val(request.input.display_colorbars), design)
    identity = (; kind = :cable, id = design.cable_id)
    payload = PreviewPayload(
        polygons,
        (),
        title,
        nothing,
        colorbars,
        PlotBuilder.LegendDefinition(enabled = request.input.display_legend),
        identity,
        PlotBuilder.ExportDefinition(
            theme = request.renderer.export_theme,
            name = design.cable_id,
            open_file = request.renderer.open_export
        ),
        nothing
    )
    return PlotBuilder.PlotPage[
        PlotBuilder.PlotPage(title, request.renderer.size, identity, payload),
    ]
end

_cable_title(::Val{false}, design) = "Cable design preview"
_cable_title(::Val{true}, design) = "Cable design preview: $(design.cable_id)"

_cable_colorbars(::Val{false}, design) = ()
_cable_colorbars(::Val{true}, design) = _colorbar_specs(_property_ranges(design)...)
