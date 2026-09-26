"""
    LineCableModelsMakieExt

Add compact high-level LineCableModels plotting methods to native Makie.
"""
module LineCableModelsMakieExt

import LineCableModels
using LineCableModels: EarthLayer, Material, RadialDielectric,
                       SeriesImpedance, ShuntAdmittance, label, nominal,
                       observables, outer_radius
import Makie
using Makie: Auto, Axis, Button, Colorbar, DataAspect, Figure,
             Fixed, GridLayout, Label, Legend, LineElement,
             Mixed,
             Observable, Outside, Rect2f, Relative, Theme, Toggle,
             colgap!, colsize!, content, errorbars!,
             hlines!, hspan!, lift, lines!,
             off, on, onany, poly!, reset_limits!, rowgap!,
             rowsize!, scatter!, stairs!, text!, to_value, translate!, update!, widths,
             with_theme
using Printf: @sprintf
using Colors: HSV, Oklab, RGB, RGBA, blue, green, red
import Dates
import Base: resize!
using Statistics: mean

import LineCableModels.Units
import LineCableModels.Engine
import LineCableModels.DataModel
import LineCableModels.Grammar
import LineCableModels.ImportExport
import LineCableModels.UQ
import Makie.GridLayoutBase
using Makie.GridLayoutBase: nrows, offsets, with_updates_suspended
import LineCableModels.Grammar:
                                request_identity

struct _Omitted end
const _omitted = _Omitted()

function current_backend_symbol()
    backend = Makie.current_backend()
    backend isa Module || return :none
    name = nameof(backend)
    name === :CairoMakie && return :cairo
    name === :GLMakie && return :gl
    name === :WGLMakie && return :wgl
    return :unknown
end

include("recipes/line_data.jl")
include("recipes/comparison_data.jl")
include("material_colors.jl")
include("recipes/preview_types.jl")
include("recipes/preview_data.jl")
include("attributes.jl")
include("series_styles.jl")
include("shell.jl")
include("controls.jl")
include("guides.jl")
include("layout.jl")
include("recipes/line_facets.jl")
include("recipes/preview_render.jl")
include("montecarlo.jl")
include("export_presentation.jl")
include("native_export.jl")

import LineCableModels.PlotBuilder: plot, preview, show_material_scale

include("plot.jl")
include("recipes/formulation_comparisons.jl")

function preview(
        design::DataModel.CableDesign;
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    return _addon_preview(
        design;
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

# The public extension method consumes backend/display choices and forwards the
# remaining preview options unchanged. DataModel retains only detached geometry
# and material attributes; Makie objects and backend state stay here.
function preview(
        designs::AbstractVector{<:DataModel.CableDesign};
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    return _addon_preview(
        designs;
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

function preview(
        system::DataModel.LineCableSystem;
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    return _addon_preview(
        system;
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

function show_material_scale(
        ; backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    return _addon_material_scale(;
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

end # module LineCableModelsMakieExt
