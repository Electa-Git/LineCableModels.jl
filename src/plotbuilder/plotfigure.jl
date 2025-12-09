"""
	build_figures(::Type{S}, nt, plot_areas) where {S<:AbstractPlotSpec}

Packs PlotArea values into PlotFigure payloads.

Default behavior:
- a single PlotFigure is created for this spec call,
- `layout` is chosen based on `grouping_mode(S)` and the number of plot areas:
	:grid   → :grid if more than one area, otherwise :single
	:none   → :single
	:overlay → :single
- figure kwargs include at least `figsize = default_figsize(S)`.

Specs that want multiple OS windows or more complex layout policies may
override this method.
"""
function build_figures(
	::Type{S},
	nt::NamedTuple,
	plot_areas::Vector{PlotArea},
) where {S <: AbstractPlotSpec}
	mode   = grouping_mode(S)
	layout = if mode === :grid && length(plot_areas) > 1
		:grid
	else
		:single
	end

	figsize = default_figsize(S)
	fig_kwargs = (figsize = figsize,)

	figure = PlotFigure("", layout, plot_areas, fig_kwargs)
	return PlotFigure[figure]
end

"""
	build_figures(::Type{S}, nt) where {S<:AbstractPlotSpec}

High-level payload builder: from normalized input `nt` to a vector of
PlotFigure values.

Pipeline:

  axes   = build_axes(S, nt)
  series = build_series(S, nt, axes)
  plot_areas  = build_plot_areas(S, nt, axes, series)
  figs   = build_figures(S, nt, plot_areas)

Specs may override `build_series`, `build_plot_areas`, or the
`build_figures(S, nt, plot_areas)` method to encode custom trace, panel, or
window grouping, while reusing the common grammar and axis machinery.
"""
function build_figures(
	::Type{S},
	nt::NamedTuple,
) where {S <: AbstractPlotSpec}
	axes = build_axes(S, nt)
	series = build_series(S, nt, axes)
	plot_areas = build_plot_areas(S, nt, axes, series)
	return build_figures(S, nt, plot_areas)
end
