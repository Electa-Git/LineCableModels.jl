module PlotBuilder

import ..UnitHandler: Units, QuantityTag

export AbstractPlotSpec, AxisSpec, SeriesSpec, ViewSpec, PageSpec, RenderSpec, UIPlot
export make_render, export_svg

function control_definitions(;
        reset::Bool = true,
        export_svg::Bool = true,
        xlog::Bool = false,
        ylog::Bool = false,
        legend::Bool = true,
        visibility::Bool = true,
        zoom::Bool = true
)
    return (; reset, export_svg, xlog, ylog, legend, visibility, zoom)
end

include("backendhandler/BackendHandler.jl")
using .BackendHandler

include("types.jl")

function make_render(::Type{S}, object; kwargs...) where {S <: AbstractPlotSpec}
    throw(MethodError(make_render, (S, object)))
end

function export_svg end

end # module PlotBuilder
