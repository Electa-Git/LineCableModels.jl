"""
    PlotBuilder

Build renderer-independent plotting recipes from domain definitions. Optional Makie
extensions render the completed `PlotRecipe` values.
"""
module PlotBuilder

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES

import ..QuantityUnits: Units, QuantityTag, display_unit, get_label
import ..LineCableModels: nominal, standard_uncertainty, validate
import ..Validation

export AbstractPlotDefinition, PlotRecipe
export UIPlot
export make_render, export_svg
export backend_available, current_backend_symbol, ensure_backend!, make_screen,
       next_fignum, renderfig, set_backend!, with_backend
export dispatch_on, input_kwargs, renderer_kwargs, input_defaults, renderer_defaults
export entitle, parse_kwargs, resolve_input, observe
export geom_axes, axis_quantity, axis_unit, axis_label
export axis_scale, axis_scales, axis_exponent, axis_attributes
export plot_kind, series_data, legend_label, series_group, series_visible,
       series_attributes
export default_title, default_figsize, layout_spec, layout_preset, page_identity
export view_key, view_placement, view_aspect, view_limits, view_attributes
export control_spec, legend_spec, colorbar_specs, status_spec, export_spec
export make_axes, make_series, make_views, make_pages, decorate, finish, validate

include("backends.jl")
include("types.jl")
include("validation.jl")
include("uiplot.jl")
include("interfaces.jl")
include("composition.jl")
include("axes.jl")
include("series.jl")
include("views.jl")
include("pages.jl")
include("render.jl")
include("base.jl")

end # module PlotBuilder
