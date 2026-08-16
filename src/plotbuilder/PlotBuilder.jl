"""
    PlotBuilder

Build backend-neutral plotting specifications from domain recipes. Optional
Makie extensions render the resulting `RenderSpec` values.
"""
module PlotBuilder

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES

import ..UnitHandler: Units, QuantityTag, display_unit, get_label
import ..Utils: to_nominal, uncertainty_value

export AbstractPlotSpec, PlotRecipe
export AbstractTrackSize, FixedTrack, RelativeTrack, ContentTrack
export GridArea, GridSpec, SlotSpec, LayoutSpec, PlacementSpec
export ControlSpec, LegendSpec, ColorbarSpec, StatusSpec, ExportSpec
export AxisSpec, SeriesSpec, ViewSpec, PageSpec, RenderSpec, UIPlot
export make_render, export_svg
export dispatch_on, input_kwargs, renderer_kwargs, input_defaults, renderer_defaults
export parse_kwargs, resolve_input, recipe_mode, grouping_mode
export page_facets, group_facets, geom_axes, axis_quantity, axis_unit, axis_label
export axis_scale, axis_scales, axis_exponent, axis_attributes
export plot_kind, series_data, legend_label, series_group, series_visible,
       series_attributes
export default_title, default_figsize, layout_spec, layout_preset, page_identity
export view_key, view_placement, view_aspect, view_limits, view_attributes
export control_spec, legend_spec, colorbar_specs, status_spec, export_spec
export make_axes, make_series, make_views, make_pages, validate

const EXPORT_THEMES = (:default, :publication)

function _validate_export_theme(value::Symbol)
    value in EXPORT_THEMES || throw(
        ArgumentError("export_theme must be :default or :publication"),
    )
    return value
end

include("backendhandler/BackendHandler.jl")
using .BackendHandler

include("types.jl")
include("grammar.jl")

"""
$(TYPEDSIGNATURES)

Export the current typed state of a `UIPlot` through an explicitly loaded
CairoMakie extension.
"""
function export_svg end

end # module PlotBuilder
