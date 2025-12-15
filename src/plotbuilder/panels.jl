"""
	build_panels(::Type{S}, nt, axes, series) where {S<:AbstractPlotSpec}

Groups Dataseries into Panel values.

Default behavior:
- all series share the same axes `axes.xaxis`, `axes.yaxis`, `axes.zaxis`,
- one Panel is created,
- `title` is taken from `default_title(S, nt)`,
- `key` is the empty NamedTuple `(; )` (no faceting semantics).

Specs that require multiple panels (e.g. grids over indices or frequencies)
should override this method and partition `series` accordingly, setting
a meaningful `key` for each Panel.
"""
function build_panels(
	::Type{S},
	nt::NamedTuple,
	axes::NamedTuple,
	series::Vector{Dataseries},
) where {S <: AbstractPlotSpec}
	title = default_title(S, nt)
	key   = (;)

	xaxis = axes.xaxis
	yaxis = axes.yaxis
	zaxis = axes.zaxis

	panel = Panel(xaxis, yaxis, zaxis, title, series, key)
	return Panel[panel]
end