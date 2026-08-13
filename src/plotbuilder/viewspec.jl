"""
	make_views(::Type{S}, nt, axes, series) where {S<:AbstractPlotSpec}

Groups SeriesSpec into ViewSpec values.

Default behavior:
- all series share the same axes `axes.xaxis`, `axes.yaxis`, `axes.zaxis`,
- one ViewSpec is created,
- `title` is taken from `default_title(S, nt)`,
- `key` is the empty NamedTuple `(; )` (no faceting semantics).

Specs that require multiple views (e.g. grids over indices or frequencies)
should override this method and partition `series` accordingly, setting
a meaningful `key` for each ViewSpec.
"""
function make_views(
	::Type{S},
	nt::NamedTuple,
	axes::NamedTuple,
	series::Vector{SeriesSpec},
) where {S <: AbstractPlotSpec}
	title = default_title(S, nt)
	key   = (;)

	xaxis = axes.xaxis
	yaxis = axes.yaxis
	zaxis = axes.zaxis

	view = ViewSpec(xaxis, yaxis, zaxis, title, series, key)
	return ViewSpec[view]
end