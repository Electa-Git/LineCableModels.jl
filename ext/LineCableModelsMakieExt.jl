module LineCableModelsMakieExt

using LineCableModels
using Makie

const PlotBuilder = LineCableModels.PlotBuilder
const BackendHandler = PlotBuilder.BackendHandler

function current_backend_symbol()
    name = nameof(Makie.current_backend())
    name === :CairoMakie && return :cairo
    name === :GLMakie && return :gl
    name === :WGLMakie && return :wgl
    return :unknown
end

renderfig(figure) = display(figure)

include(joinpath(
    @__DIR__,
    "..",
    "src",
    "plotbuilder",
    "uicomponents",
    "UIComponents.jl"
))
using .UIComponents

import LineCableModels.Engine: plot
import LineCableModels.DataModel: preview, show_material_scale

function _scale_symbol(value)
    value isa Symbol && return value
    value === Makie.identity && return :linear
    value === Makie.log10 && return :log10
    throw(ArgumentError("axis scale must be :linear, :log10, Makie.identity, or Makie.log10"))
end

function plot(
        object::LineCableModels.SeriesImpedance,
        frequencies;
        backend = nothing,
        display_plot::Bool = true,
        xscale = :linear,
        yscale = :linear,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        object;
        frequencies,
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        kwargs...
    )
    return UIComponents.build(render_spec; backend, display = display_plot)
end

function Makie.plot(object::LineCableModels.SeriesImpedance, frequencies; kwargs...)
    plot(object, frequencies; kwargs...)
end

function plot(
        object::LineCableModels.ShuntAdmittance,
        frequencies;
        backend = nothing,
        display_plot::Bool = true,
        xscale = :linear,
        yscale = :linear,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        object;
        frequencies,
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        kwargs...
    )
    return UIComponents.build(render_spec; backend, display = display_plot)
end

function Makie.plot(object::LineCableModels.ShuntAdmittance, frequencies; kwargs...)
    plot(object, frequencies; kwargs...)
end

function plot(
        parameters::LineCableModels.LineParameters;
        backend = nothing,
        display_plot::Bool = true,
        xscale = :linear,
        yscale = :linear,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        parameters;
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        kwargs...
    )
    return UIComponents.build(render_spec; backend, display = display_plot)
end

function Makie.plot(parameters::LineCableModels.LineParameters; kwargs...)
    plot(parameters; kwargs...)
end

function _quantity_symbol(expression)
    expression isa Symbol && return expression, nothing
    if expression isa Expr && expression.head === :ref && length(expression.args) == 4
        return Symbol(expression.args[1]), Tuple(Int.(expression.args[2:4]))
    end
    throw(ArgumentError("use a quantity Symbol and optional ijk=(i,j,k)"))
end

function plot(
        result::Union{LineCableModels.CableConstantsMC, LineCableModels.LineParametersMC},
        expression = :R;
        ijk = nothing,
        backend = nothing,
        display_plot::Bool = true,
        kwargs...
)
    quantity, parsed_indices = _quantity_symbol(expression)
    selection = ijk === nothing ? parsed_indices : ijk
    render_spec = PlotBuilder.make_render(
        LineCableModels.UQ.MCDistributionPlotSpec,
        result;
        quantity,
        ijk = selection,
        kwargs...
    )
    return only(UIComponents.build(render_spec; backend, display = display_plot))
end

function Makie.plot(
        result::Union{LineCableModels.CableConstantsMC, LineCableModels.LineParametersMC},
        expression = :R;
        kwargs...
)
    plot(result, expression; kwargs...)
end

function preview(
        design::LineCableModels.CableDesign;
        backend = nothing,
        display_plot::Bool = true,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.DataModel.CablePreviewPlotSpec,
        design;
        kwargs...
    )
    return only(UIComponents.build(render_spec; backend, display = display_plot))
end

function preview(
        system::LineCableModels.LineCableSystem;
        backend = nothing,
        display_plot::Bool = true,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.DataModel.SystemPreviewPlotSpec,
        system;
        kwargs...
    )
    return only(UIComponents.build(render_spec; backend, display = display_plot))
end

function show_material_scale(
        ; backend = nothing,
        display_plot::Bool = true,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.DataModel.MaterialScalePlotSpec,
        nothing;
        kwargs...
    )
    return only(UIComponents.build(render_spec; backend, display = display_plot))
end

end # module LineCableModelsMakieExt
