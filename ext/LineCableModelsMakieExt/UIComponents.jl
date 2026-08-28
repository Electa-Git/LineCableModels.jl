"""
    LineCableModelsMakieExt.UIComponents

Draw detached PlotBuilder pages through the standard Makie shell.
"""
module UIComponents

using Makie
using Dates
using LinearAlgebra: diag
using Printf: @sprintf

import LineCableModels
import LineCableModels.PlotBuilder
import LineCableModels.PlotBuilder: validate_export_theme
using LineCableModels: nominal, uncertainty
using LineCableModels.Units: label
using LineCableModels.PlotBuilder:
                                   ColorbarDefinition, ExportDefinition,
                                   LegendDefinition, PlotPage, PlotRecipe, UIPlot

export build
public build_widget!, toolbar_button!, bind_widget_callback!
public place_legend!, place_colorbars!

include("context.jl")
include("uicomponents/theme.jl")
include("uicomponents/axes.jl")
include("uicomponents/drawing.jl")
include("uicomponents/docks.jl")
include("shell.jl")
include("uicomponents/export.jl")
include("native.jl")
include("lineparameters.jl")
include("montecarlo.jl")
include("previews.jl")

end # module UIComponents
