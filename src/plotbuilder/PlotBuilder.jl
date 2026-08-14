module PlotBuilder

import ..UnitHandler: Units, QuantityTag, display_unit, get_label

export AbstractPlotSpec, AxisSpec, SeriesSpec, ViewSpec, PageSpec, RenderSpec, UIPlot
export make_render, export_svg
export dispatch_on, input_kwargs, renderer_kwargs, input_defaults, renderer_defaults
export parse_kwargs, resolve_input, recipe_mode, grouping_mode
export page_facets, group_facets, geom_axes, axis_quantity, axis_unit, axis_label
export axis_scale, plot_kind, series_data, legend_label, series_attributes
export default_title, default_figsize, figure_layout, enable_logscale
export view_key, view_attributes, page_kwargs, make_axes, make_series, make_views,
       make_pages

const EXPORT_THEMES = (:default, :publication)

function _validate_export_theme(value::Symbol)
    value in EXPORT_THEMES || throw(
        ArgumentError("export_theme must be :default or :publication"),
    )
    return value
end

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
include("grammar.jl")

function export_svg end

end # module PlotBuilder
