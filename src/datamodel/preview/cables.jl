# Stage 1 — reject an empty collection while the source is still being
# entitled. Later stages may therefore assume that at least one panel exists.
function PlotBuilder.entitle(
        ::Type{CableCollectionPreviewPlotDefinition},
        designs::AbstractVector{<:CableDesign}
)
    isempty(designs) && throw(ArgumentError(
        "a cable collection preview requires at least one design",
    ))
    return designs
end

# Stage 2 — declare the complete semantic keyword surface. `layout = nothing`
# asks the definition to select a near-square grid; `(rows, columns)` fixes it.
function PlotBuilder.input_defaults(
        ::Type{CableCollectionPreviewPlotDefinition},
        ::AbstractVector{<:CableDesign}
)
    return (; layout = nothing, display_colorbars = true)
end

# Stage 3 — renderer choices stay separate from semantic layout. `size` is the
# complete canvas size and remains caller-overridable through the common parser.
function PlotBuilder.renderer_defaults(
        ::Type{CableCollectionPreviewPlotDefinition},
        ::AbstractVector{<:CableDesign}
)
    return (; size = (1200, 900))
end

# PlotBuilder's common `parse` stage consumes the two declarations above,
# rejects unknown keywords, and separates semantic input from renderer input.
# This definition does not replace that fixed grammar step.

# Stage 4 — normalize the optional layout into one fixed `(rows, columns)`
# value. The near-square rule fills rows left-to-right and leaves at most the
# trailing slots empty. Explicit layouts may be larger but never too small.
function PlotBuilder.resolve(
        ::Type{CableCollectionPreviewPlotDefinition},
        designs::AbstractVector{<:CableDesign},
        request::NamedTuple
)
    input = request.input
    input.display_colorbars isa Bool || throw(
        ArgumentError("display_colorbars must be Bool"),
    )
    request.renderer.size isa Tuple{Int, Int} || throw(
        ArgumentError("size must be a tuple of two integers"),
    )
    all(>(0), request.renderer.size) || throw(
        ArgumentError("size dimensions must be positive"),
    )

    if input.layout === nothing
        columns = ceil(Int, sqrt(length(designs)))
        rows = cld(length(designs), columns)
    else
        input.layout isa Tuple && length(input.layout) == 2 &&
            all(value -> value isa Integer && !(value isa Bool), input.layout) || throw(
            ArgumentError("layout must be a tuple of two positive integers or nothing"),
        )
        rows, columns = Int.(input.layout)
        rows > 0 && columns > 0 || throw(
            ArgumentError("layout dimensions must be positive"),
        )
        rows * columns >= length(designs) || throw(DimensionMismatch(
            "layout provides $(rows * columns) slots for $(length(designs)) cable designs",
        ))
    end

    return merge(request, (; input = merge(input, (; layout = (rows, columns)))))
end

# Stage 5 — materialise all renderer-neutral panel data. Each panel reuses the
# scalar `PreviewPayload`; only its position and cable identifier are contingent
# collection information, so they remain NamedTuple fields rather than new
# owned result types. Labels are disabled here because this recipe never owns a
# legend.
function PlotBuilder.fetch(
        ::Type{CableCollectionPreviewPlotDefinition},
        designs::AbstractVector{<:CableDesign},
        request::NamedTuple
)
    rows, columns = request.input.layout
    panels = map(enumerate(designs)) do (index, design)
        position = (cld(index, columns), mod1(index, columns))
        payload = PreviewPayload(
            _design_shapes(design, 0.0, 0.0; display_legend = false),
            (),
            nothing,
            nothing
        )
        return (; position, title = design.cable_id, payload)
    end

    # One scan over the complete source gives the shell three shared material
    # scales. The Makie extension receives no design objects and cannot drift
    # into calculating a different scale for each subplot.
    colorbars = request.input.display_colorbars ?
                _colorbar_specs(_property_ranges(designs)...) : ()
    identity = (;
        kind = :cable_collection,
        ids = Tuple(design.cable_id for design in designs),
        layout = (rows, columns)
    )
    payload = (; panels, layout = (rows, columns))

    # Stage 6 — complete one page. Legend absence and shared colorbars are page
    # declarations; the Makie extension's definition-dispatched placement stage
    # decides where to draw them. Subplot data stays inside the definition-owned
    # payload above.
    return PlotBuilder.PlotPage[
        PlotBuilder.PlotPage(
            "Cable design previews",
            request.renderer.size,
            identity,
            payload;
            legend = PlotBuilder.LegendDefinition(enabled = false),
            colorbars,
            export_definition = PlotBuilder.ExportDefinition(
                theme = request.renderer.export_theme,
                name = "cable_design_previews",
                open_file = request.renderer.open_export
            )
        ),
    ]
end

# Stage 7 — no `finish` specialization is needed. PlotBuilder's common final
# stage checks page cardinality and seals this page into a `PlotRecipe`.
