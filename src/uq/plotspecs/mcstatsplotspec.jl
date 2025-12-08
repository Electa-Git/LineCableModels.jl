
using ..UnitHandler:
	QuantityTag, default_unit, display_unit, get_label, get_symbol, scale_factor
using ..PlotBuilder: AbstractPlotSpec, PlotBuilder

struct MCStatsPlotSpec <: AbstractPlotSpec end

# --- Traits ---------------------------------------------------------------

PlotBuilder.plot_kind(::Type{MCStatsPlotSpec}) = :line

# no logscales for now.
PlotBuilder.enable_logscale(::Type{MCStatsPlotSpec}) = ()

# Dispatch on MC container
PlotBuilder.dispatch_on(::Type{MCStatsPlotSpec}) = LineParametersMC

PlotBuilder.default_figsize(::Type{MCStatsPlotSpec}) = (800, 400)

# Geometric axes: 2D line plot
PlotBuilder.geom_axes(::Type{MCStatsPlotSpec}) = (:x, :y)

# Index convention: all MC stats are per-(i,j) element over frequency
PlotBuilder.index_keys(::Type{MCStatsPlotSpec})  = (:i, :j, :k)
PlotBuilder.ranged_keys(::Type{MCStatsPlotSpec}) = (:k,)

# Grouping: overlay all on single plot
PlotBuilder.grouping_mode(::Type{MCStatsPlotSpec}) = :overlay

# X is always frequency; Y will depend on user kwarg, so the valid possible quantities are defined below. 
PlotBuilder.axis_quantity(::Type{MCStatsPlotSpec}, ::Val{:x}, ::Val{:f}) =
	QuantityTag{:freq}()
PlotBuilder.axis_quantity(::Type{MCStatsPlotSpec}, ::Val{:y}, ::Val{:R}) =
	QuantityTag{:resistance}()
PlotBuilder.axis_quantity(::Type{MCStatsPlotSpec}, ::Val{:y}, ::Val{:L}) =
	QuantityTag{:inductance}()
PlotBuilder.axis_quantity(::Type{MCStatsPlotSpec}, ::Val{:y}, ::Val{:C}) =
	QuantityTag{:capacitance}()
PlotBuilder.axis_quantity(::Type{MCStatsPlotSpec}, ::Val{:y}, ::Val{:G}) =
	QuantityTag{:conductance}()

PlotBuilder.data_container(::Type{MCStatsPlotSpec}, ::Val{:x}) = nothing         # obj.f
PlotBuilder.data_container(::Type{MCStatsPlotSpec}, ::Val{:y}) = :stats          # obj.stats[Sym]


# Define plot title
function PlotBuilder.default_title(::Type{MCStatsPlotSpec}, nt::NamedTuple)
	qx = nt.x_quantity
	qy = nt.y_quantity
	stat_sym = nt.stat

	y_label = get_label(qy)  # "Series resistance"
	x_label = get_label(qx)  # "Frequency"

	return string(y_label, " ", String(stat_sym), " vs. ", x_label)
end

# Define legend labels
function PlotBuilder.legend_labels(::Type{MCStatsPlotSpec}, nt::NamedTuple)
	qy       = nt.y_quantity
	stat_sym = nt.stat
	i        = nt.i
	j        = nt.j

	y_label = get_symbol(qy)

	entry = string(y_label, "[", i, ",", j, "] ", String(stat_sym))
	return [entry]
end

# Semantic knobs:
#   :stat  → which field from stats NamedTuple
PlotBuilder.input_kwargs(::Type{MCStatsPlotSpec}) =
	(:stat,)

# Backend knobs this spec forwards to Makie
PlotBuilder.backend_kwargs(::Type{MCStatsPlotSpec}) =
	(:color, :linewidth, :linestyle, :marker, :markersize)

# Defaults for semantic knobs, given the dispatched object
PlotBuilder.input_defaults(::Type{MCStatsPlotSpec}, ::LineParametersMC) = (
	xdata = :f,
	ydata = :R,
	stat = :mean,
	# i,j,k are handled via index_keys + normalize_kwargs (default 1 or :)
)

# Defaults for backend knobs
PlotBuilder.backend_defaults(::Type{MCStatsPlotSpec}, ::LineParametersMC) = (
	color      = :blue,
	linewidth  = 2,
	marker     = :circle,
	markersize = 6,
)

function PlotBuilder.build_payloads(::Type{S}, nt::NamedTuple) where {S <: MCStatsPlotSpec}
	obj       = nt.obj
	ydata_src = nt.ydata        # :R, :L, :C, :G
	field_src = nt.stat
	i, j      = nt.i, nt.j
	backend   = nt.backend

	xraw = getproperty(obj, nt.xdata)    # lp.f

	stats_arr = obj.stats[ydata_src]
	stats     = stats_arr[i, j, :]
	yraw      = [s[field_src] for s in stats]

	length(yraw) == length(xraw) ||
		error("stats length ($(length(yraw))) != length(f) ($(length(xraw)))")

	axes  = PlotBuilder.build_axes(S, nt)
	xaxis = axes.xaxis
	yaxis = axes.yaxis
	zaxis = axes.zaxis  # `nothing` here, but grammar-consistent

	sx = scale_factor(xaxis.quantity, xaxis.units)
	sy = scale_factor(yaxis.quantity, yaxis.units)

	xdata  = xraw .* sx
	ydata  = yraw .* sy
	title  = PlotBuilder.default_title(S, nt)
	legend = PlotBuilder.legend_labels(S, nt)

	payload = (
		xdata  = xdata,
		ydata  = ydata,
		zdata  = nothing,
		xaxis  = xaxis,
		yaxis  = yaxis,
		zaxis  = zaxis,
		title  = title,
		legend = legend,
		kwargs = backend,
	)

	return [payload]
end
