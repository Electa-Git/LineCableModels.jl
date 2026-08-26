"""
    PlotBuilder

Build renderer-independent plotting recipes from domain definitions. Optional Makie
extensions render the completed `PlotRecipe` values.
"""
module PlotBuilder

using DocStringExtensions: SIGNATURES, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES

export AbstractPlotDefinition, PlotPage, PlotRecipe
export LegendDefinition, ColorbarDefinition, ExportDefinition,
       AbstractWidgetDefinition
export UIPlot
export make_render, export_svg, plotwindow, axis!, register!
export backend_available, current_backend_symbol, ensure_backend!, make_screen,
       next_fignum, renderfig, set_backend!, with_backend
export dispatch_on, input_kwargs, renderer_kwargs, input_defaults, renderer_defaults
export entitle, parse, resolve, fetch, finish

include("backends.jl")
include("types.jl")
include("uiplot.jl")
include("interfaces.jl")
include("render.jl")
include("base.jl")

end # module PlotBuilder
