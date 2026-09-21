"""
    PlotBuilder

Provide plotting functions and live figure handles for the Makie extension.
Load a Makie backend to draw figures; returned handles expose the Makie objects
for further editing.
"""
module PlotBuilder

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES

export UIPlot, plot, preview, show_material_scale, export_svg
export figurelegend!, panellegend!, figuretitle!, paneltitle!
export figurecolorbars!, axisscale!, resetview!, addwidget!, removewidget!
export plotwindow, materialcolors, materialscale!

include("handle.jl")
include("interfaces.jl")

end # module PlotBuilder
