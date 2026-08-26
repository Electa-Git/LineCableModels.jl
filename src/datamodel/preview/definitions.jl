struct CablePreviewPlotDefinition <: PlotBuilder.AbstractPlotDefinition end
struct SystemPreviewPlotDefinition <: PlotBuilder.AbstractPlotDefinition end
struct MaterialScalePlotDefinition <: PlotBuilder.AbstractPlotDefinition end

"""
$(TYPEDEF)

Store one filled polygon prepared for a cable or system preview.

$(TYPEDFIELDS)
"""
struct PreviewPolygon{G, C, S}
    "Renderer-neutral polygon geometry."
    geometry::G
    "Legend text, or `nothing` when the polygon introduces no entry."
    label::Union{Nothing, String}
    "Identity shared by polygons controlled by one legend entry."
    group::Symbol
    "Final fill color after material-property mapping."
    color::C
    "Polygon outline color."
    stroke::S
    "Polygon outline width in pixels."
    width::Float64
end

"""
$(TYPEDEF)

Store one horizontal reference prepared for a system preview.

$(TYPEDFIELDS)
"""
struct PreviewReferenceLine{V, C}
    "Vertical coordinates in meters."
    values::V
    "Legend-group identity."
    group::Symbol
    "Line color."
    color::C
    "Line width in pixels."
    width::Float64
end

"""
$(TYPEDEF)

Store the detached geometry and shell declarations needed to redraw one preview.

$(TYPEDFIELDS)
"""
struct PreviewPayload{P, R, L, C, K, E, S}
    "Prepared filled polygons in draw order."
    polygons::P
    "Prepared horizontal references in draw order."
    references::R
    "Displayed page and axis title."
    title::String
    "Explicit `(x, y)` limits, or `nothing` for fitted limits."
    limits::L
    "Material color scales shown by the standard shell."
    colorbars::C
    "Legend behavior supplied to the standard shell."
    legend::PlotBuilder.LegendDefinition
    "Semantic page identity."
    key::K
    "SVG export behavior supplied to the standard shell."
    export_definition::E
    "Captured runtime state used for current-state SVG replay."
    runtime::S
end

"""
$(TYPEDEF)

Store the material scales displayed by a colorbar-only page.

$(TYPEDFIELDS)
"""
struct MaterialScalePayload{C, E}
    "Material color scales shown by the standard shell."
    colorbars::C
    "Legend behavior supplied to the standard shell."
    legend::PlotBuilder.LegendDefinition
    "SVG export behavior supplied to the standard shell."
    export_definition::E
end
