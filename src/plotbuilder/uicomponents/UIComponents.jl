module UIComponents

using Makie

import ..BackendHandler
using ..PlotUIComponents: PlotAssembly

import ..PlotBuilder: AbstractPlotSpec, RenderSpec, PageSpec, ViewSpec, SeriesSpec, AxisSpec

export build,
       UIContext, UILayoutSpec, UIContainerSpec, UISlotSpec,
       UIFigure, UIPanel, PlotAssembly

export display!

include("themes.jl")
include("types.jl")
include("layoutspecs.jl")
include("actions.jl")
include("draw.jl")
include("widgets.jl")
include("pipeline.jl")

end # module UIComponents
