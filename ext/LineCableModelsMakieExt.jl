module LineCableModelsMakieExt

using LineCableModels
using Makie
using Dates: format, now

const PlotBuilder = LineCableModels.PlotBuilder
const BackendHandler = LineCableModels.PlotBuilder.BackendHandler

function current_backend_symbol()
    name = nameof(Makie.current_backend())
    name === :CairoMakie && return :cairo
    name === :GLMakie && return :gl
    name === :WGLMakie && return :wgl
    return :unknown
end

renderfig(fig) = display(fig)

include(joinpath(
    @__DIR__, "..", "src", "plotbuilder", "plotuicomponents", "PlotUIComponents.jl"))
include(joinpath(@__DIR__, "..", "src", "plotbuilder", "uicomponents", "UIComponents.jl"))

using LineCableModels.PlotBuilder: AbstractPlotSpec
using LineCableModels.PlotBuilder.BackendHandler: next_fignum
using .PlotUIComponents: ControlButtonSpec, ControlReaction, ICON_TTF, MI_REFRESH, MI_SAVE,
                         _make_window, _run_plot_pipeline, with_icon, with_plot_theme

include(joinpath(@__DIR__, "..", "src", "plotbuilder", "plotspecs.jl"))

module DataModelPreview

using Makie
using Colors
using Printf
using Dates
using Statistics
using LineCableModels.DataModel
using LineCableModels.DataModel.BaseParams: calc_circstrands_coords
import LineCableModels.DataModel: AbstractCablePart, preview
using LineCableModels.Utils: is_in_testset, to_nominal
using LineCableModels.PlotBuilder.BackendHandler: current_backend_symbol, ensure_backend!,
                                                  next_fignum, renderfig
using ..PlotUIComponents: ICON_TTF, MI_REFRESH, MI_SAVE, gl_screen, with_icon

include(joinpath(@__DIR__, "..", "src", "datamodel", "preview.jl"))

end # module DataModelPreview

module EnginePlots

using Makie
using Measurements: Measurements
using LineCableModels.Engine
import LineCableModels.Engine: get_description, plot
using LineCableModels.Commons: ModalDomain, PhaseDomain, domain
using LineCableModels.PlotBuilder
using ..PlotUIComponents

include(joinpath(@__DIR__, "..", "src", "engine", "plot.jl"))

end # module EnginePlots

module UQPlots

using Makie
using Printf
using Distributions
using StatsBase
using LineCableModels.UQ: CableDesignMC, LineParametersMC, LineParametersPDF
using LineCableModels.PlotBuilder
using ..EnginePlots
using ..PlotUIComponents

include(joinpath(@__DIR__, "..", "src", "uq", "plot.jl"))

end # module UQPlots

end # module LineCableModelsMakieExt
