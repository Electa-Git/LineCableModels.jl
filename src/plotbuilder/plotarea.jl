"""
	build_plot_areas(::Type{S}, nt, axes, series) where {S<:AbstractPlotSpec}

Groups PlotSeries into PlotArea values.

Default behavior:
- all series share the same axes `axes.xaxis`, `axes.yaxis`, `axes.zaxis`,
- one PlotArea is created,
- `title` is taken from `default_title(S, nt)`,
- `key` is the empty NamedTuple `(; )` (no faceting semantics).

Specs that require multiple panels (e.g. grids over indices or frequencies)
should override this method and partition `series` accordingly, setting
a meaningful `key` for each PlotArea.
"""
function build_plot_areas(
	::Type{S},
	nt::NamedTuple,
	axes::NamedTuple,
	series::Vector{PlotSeries},
) where {S <: AbstractPlotSpec}
	title = default_title(S, nt)
	key   = (;)

	xaxis = axes.xaxis
	yaxis = axes.yaxis
	zaxis = axes.zaxis

	plot_area = PlotArea(xaxis, yaxis, zaxis, title, series, key)
	return PlotArea[plot_area]
end