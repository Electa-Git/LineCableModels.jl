
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
PlotBuilder.index_keys(::Type{MCStatsPlotSpec}) = (:i, :j, :k)
PlotBuilder.ranged_keys(::Type{MCStatsPlotSpec}) = (:k,)

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
	key = nt.field

	y_label = get_label(qy)  # "Series resistance"
	x_label = get_label(qx)  # "Frequency"

	return string(y_label, " ", String(key), " vs. ", x_label)
end

# Define legend labels
function PlotBuilder.legend_labels(::Type{MCStatsPlotSpec}, nt::NamedTuple)
	qy = nt.y_quantity
	key = nt.field
	i = nt.i
	j = nt.j

	y_label = get_symbol(qy)

	entry = string(y_label, "[", i, ",", j, "] ", String(key))
	return [entry]
end

# Semantic knobs:
#   :field  → which field from stats NamedTuple
PlotBuilder.input_kwargs(::Type{MCStatsPlotSpec}) = (:field,)

# Backend knobs this spec forwards to Makie
PlotBuilder.renderer_kwargs(::Type{MCStatsPlotSpec}) = ()

# Defaults for semantic knobs, given the dispatched object
PlotBuilder.input_defaults(::Type{MCStatsPlotSpec}, ::LineParametersMC) = (
	x = :f,
	y = :R,
	field = :mean,
	# i,j,k are handled via index_keys + parse_kwargs (default 1 or :)
)

PlotBuilder.select_field(::Type{MCStatsPlotSpec}, ::Val{:x}) = nothing
PlotBuilder.select_field(::Type{MCStatsPlotSpec}, ::Val{:y}) = :field # or :mean directly

# Defaults for renderer knobs
PlotBuilder.renderer_defaults(::Type{MCStatsPlotSpec}, ::LineParametersMC) = ()
