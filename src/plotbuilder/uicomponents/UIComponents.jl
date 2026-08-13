module UIComponents

using Makie

import ..BackendHandler

import ..PlotBuilder: AbstractPlotSpec, RenderSpec, PageSpec, ViewSpec, SeriesSpec, AxisSpec

export build, export_svg!,
	UIContext, UILayoutSpec, UIContainerSpec, UISlotSpec,
	UIFigure, UIPanel, PlotAssembly

export build_context, display!

include("themes.jl")
include("types.jl")
include("layoutspecs.jl")
include("actions.jl")
include("draw.jl")
include("widgets.jl")
include("pipeline.jl")

end # module UIComponents
