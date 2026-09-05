"""
    PlotBuilder

Provide the thin, renderer-optional plotting surface used by the Makie
extension.  This module owns live plot handles and addon entry points only; it
does not define renderer-independent plot specifications or scientific data.
"""
module PlotBuilder

using DocStringExtensions: TYPEDEF, TYPEDFIELDS

export UIPlot, plot, preview, show_material_scale, export_svg
export figurelegend!, panellegend!, figuretitle!, paneltitle!
export plotwindow, materialcolors, materialscale!

include("handle.jl")
include("interfaces.jl")

end # module PlotBuilder
