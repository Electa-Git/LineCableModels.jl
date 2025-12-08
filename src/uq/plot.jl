using Makie
import Makie: plot

using ..Engine:
	Engine,
	SeriesImpedance,
	ShuntAdmittance,
	UnitSpec,
	LP_FIG_SIZE,
	quantity_scale,
	length_scale,
	normalize_quantity_units,
	composite_unit,
	resolve_quantity_prefix,
	autoscale_axis,
	_axis_label,
	_ICON_FN,
	get_description,
	get_unit_symbol,
	ComponentMetadata

include("../plotbuilder/plothelpers.jl")

struct MCPlotHistSpec
	quantity::Symbol
	symbol::String
	title::String
	xlabel::String
	ylabel::String
	values::Union{Nothing, Vector{<:Real}}
	pdf_obj::Union{Nothing, LineParametersPDF}
	bins::Vector{<:Real}
	x_exp::Int
	fig_size::Union{Nothing, Tuple{Int, Int}}
	normalization::Symbol
	data::Symbol
	mode::Symbol
end



function _mc_quantity_metadata()
	sdesc = get_description(SeriesImpedance(zeros(1, 1, 1)))
	symb = Engine.get_symbol(SeriesImpedance(zeros(1, 1, 1)))
	sunit = get_unit_symbol(SeriesImpedance(zeros(1, 1, 1)))

	adesc = get_description(ShuntAdmittance(zeros(1, 1, 1)))
	asymb = Engine.get_symbol(ShuntAdmittance(zeros(1, 1, 1)))
	aunit = get_unit_symbol(ShuntAdmittance(zeros(1, 1, 1)))

	return Dict(
		:R => ComponentMetadata(
			:resistance,
			:resistance,
			symb.resistance,
			sdesc.resistance,
			symb.resistance,
			UnitSpec(sunit.resistance, true),
		),
		:L => ComponentMetadata(
			:inductance,
			:inductance,
			symb.inductance,
			sdesc.inductance,
			symb.inductance,
			UnitSpec(sunit.inductance, true),
		),
		:C => ComponentMetadata(
			:capacitance,
			:capacitance,
			asymb.capacitance,
			adesc.capacitance,
			asymb.capacitance,
			UnitSpec(aunit.capacitance, true),
		),
		:G => ComponentMetadata(
			:conductance,
			:conductance,
			asymb.conductance,
			adesc.conductance,
			asymb.conductance,
			UnitSpec(aunit.conductance, true),
		),
	)
end

const _MC_META = _mc_quantity_metadata()

_quantity_metadata(q::Symbol) =
	get(_MC_META, q) do
		Base.error("Unsupported quantity $(q); use one of :R, :L, :C, :G")
	end

function _parse_values_ref(values, ijk)
	if ijk === nothing
		if values isa Expr && values.head === :ref && length(values.args) == 4
			q = values.args[1]
			q isa Symbol ||
				Base.error("Expected values symbol as first argument in $(values)")
			i, j, k = values.args[2:4]
			return q, (Int(i), Int(j), Int(k))
		else
			Base.error(
				"Provide values as Expr like :R[1,1,1] or pass indices via `ijk = (i,j,k)`",
			)
		end
	else
		ijk isa NTuple{3, Int} ||
			Base.error("ijk must be NTuple{3,Int}, got $(typeof(ijk))")
		values isa Symbol ||
			Base.error("values must be Symbol when ijk is provided; got $(typeof(values))")
		return values, ijk
	end
end

function _build_hist_spec(
	obj::LineParametersMC,
	values_expr;
	ijk::Union{Nothing, NTuple{3, Int}} = nothing,
	length_unit::Symbol = :kilo,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	per_length::Bool = true,
	quantity_units = nothing,
	nbins::Union{Nothing, Int} = nothing,
	normalization::Symbol = :none,
	data::Symbol = :samples,
	mode::Symbol = :hist,
)
	# Force data=:both if mode is not :hist
	plot_data = (mode == :hist) ? data : :both

	vals = nothing
	pdf_obj = nothing

	values_sym, _ijk = _parse_values_ref(values_expr, ijk)
	meta = _quantity_metadata(values_sym)
	units = normalize_quantity_units(quantity_units)
	q_prefix = resolve_quantity_prefix(meta.quantity, units)
	c_scale = (per_length ? length_scale(length_unit) : 1.0) * quantity_scale(q_prefix)
	i, j, k = _ijk

	# --- Load Data (uses plot_data) ---
	if plot_data == :samples || plot_data == :both
		obj.samples === nothing &&
			Base.error("mode=:$mode requires samples, but none are available.")
		samps = getfield(obj.samples, values_sym)

		max_i, max_j, max_k, _ = size(samps)
		(1 <= i <= max_i && 1 <= j <= max_j && 1 <= k <= max_k) || Base.error(
			"indices (i=$(i), j=$(j), k=$(k)) out of bounds for samples size $(size(samps))",
		)

		raw_vals = @view samps[i, j, k, :]
		vals = collect(raw_vals) .* c_scale
	end
	if plot_data == :pdf || plot_data == :both
		obj.pdf === nothing &&
			Base.error("mode=:$mode requires PDF, but none is available.")
		pdfs = getfield(obj.pdf, values_sym)

		max_i, max_j, max_k = size(pdfs)
		(1 <= i <= max_i && 1 <= j <= max_j && 1 <= k <= max_k) || Base.error(
			"indices (i=$(i), j=$(j), k=$(k)) out of bounds for pdf size $(size(pdfs))",
		)

		raw_pdf = pdfs[i, j, k]
		scaled_edges = raw_pdf.edges .* c_scale
		scaled_dens = raw_pdf.dens ./ c_scale
		pdf_obj = LineParametersPDF(scaled_edges, scaled_dens)
	end

	# --- Binning (Conditional) & Scaling ---
	local bin_edges::Vector{<:Real}
	local current_norm::Symbol
	local x_exp::Int

	if mode == :hist
		if pdf_obj !== nothing
			bin_edges = pdf_obj.edges
			current_norm = :pdf
		elseif vals !== nothing
			current_norm = normalization
			hist_nbins = isnothing(nbins) ? _auto_nbins(vals) : nbins
			h_fit = fit(Histogram, vals; nbins = hist_nbins, closed = :left)
			bin_edges = collect(h_fit.edges[1])
		else
			error("No data (samples or PDF) to plot.")
		end
	else
		bin_edges = Float64[] # Not used
		current_norm = :none  # Not used
	end

	# Autoscale is *always* needed
	if vals !== nothing
		_, x_exp = autoscale_axis(vals)
	else
		_, x_exp = autoscale_axis(pdf_obj.edges)
	end

	# --- Labels and Title (Conditional) ---
	local title::String
	local xlabel::String
	local ylabel::String

	xlabel_unit = composite_unit(q_prefix, meta.unit.symbol, per_length, length_unit)
	base_xlabel = string(meta.axis_label, " [", xlabel_unit, "]")
	freq_str = @sprintf("%.4g", obj.f[k])

	if mode == :hist
		title = string(meta.title, " histogram @ f=", freq_str, " Hz")
		xlabel = base_xlabel
		ylabel =
			current_norm == :none ? "count" :
			(current_norm == :pdf ? "density" : String(current_norm))
	elseif mode == :ecdf
		title = string(meta.title, " CDF @ f=", freq_str, " Hz")
		xlabel = base_xlabel
		ylabel = "cumulative probability"
	elseif mode == :qq
		title = string(meta.title, " Q-Q plot @ f=", freq_str, " Hz")
		xlabel = "sample quantiles"
		ylabel = "model quantiles"
	end

	return MCPlotHistSpec(
		values_sym, meta.symbol, title, xlabel, ylabel,
		vals, pdf_obj, bin_edges, x_exp,
		fig_size, current_norm, plot_data, mode,
	)
end

function _hist_specs(
	obj::LineParametersMC,
	values_expr;
	ijk::Union{Nothing, NTuple{3, Int}} = nothing,
	length_unit::Symbol = :kilo,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	per_length::Bool = true,
	quantity_units = nothing,
	nbins::Union{Nothing, Int} = nothing,
	normalization::Symbol = :none,
	data::Symbol = :samples,
	mode::Symbol = :hist,
)
	spec = _build_hist_spec( # Calls LP-MC builder
		obj,
		values_expr;
		ijk = ijk,
		length_unit = length_unit,
		fig_size = fig_size,
		per_length = per_length,
		quantity_units = quantity_units,
		nbins = nbins,
		normalization = normalization,
		data = data,
		mode = mode,
	)
	return [spec]
end

function _default_export_path_hist(spec::MCPlotHistSpec)
	base_title = strip(spec.title)
	name = _sanitize_filename_plot(base_title)
	timestamp = Dates.format(Dates.now(), EXPORT_TIMESTAMP_FORMAT)
	filename = string(name, "_", timestamp, ".", EXPORT_EXTENSION)
	return joinpath(pwd(), filename)
end

function _save_hist_export(spec::MCPlotHistSpec, axis)
	fig = _build_hist_export(spec)
	trim!(fig.layout)
	path = _default_export_path_hist(spec)
	Makie.save(path, fig)
	return path
end

function _build_hist_export(spec::MCPlotHistSpec)
	backend_ctx = _make_window(
		BackendHandler,
		:cairo;
		icons = _ICON_FN,
		icons_font = ICON_TTF,
		interactive_override = false,
		use_latex_fonts = true,
	)
	pipeline_kwargs =
		spec.fig_size === nothing ?
		(; initial_status = "") :
		(; fig_size = spec.fig_size, initial_status = "")
	assembly = with_plot_theme(backend_ctx; mode = :export) do
		_run_plot_pipeline(
			backend_ctx,
			(fig_ctx, ctx, axis) -> _build_hist_plot!(fig_ctx, ctx, axis, spec);
			pipeline_kwargs...,
		)
	end
	ensure_export_background!(assembly.figure)
	return assembly.figure
end

function _render_hist_spec(
	spec::MCPlotHistSpec; # <-- Unified spec
	backend = nothing,
	display_plot::Bool = true,
)
	n = next_fignum()
	backend_ctx = _make_window(
		BackendHandler,
		backend;
		title = "Fig. $(n) – $(spec.title)",
		icons = _ICON_FN,
		icons_font = ICON_TTF,
	)
	pipeline_kwargs =
		spec.fig_size === nothing ?
		(; initial_status = " ") :
		(; fig_size = spec.fig_size, initial_status = " ")
	assembly = with_plot_theme(backend_ctx) do
		_run_plot_pipeline(
			backend_ctx,
			# Calls the single, unified _build_hist_plot!
			(fig_ctx, ctx, axis) -> _build_hist_plot!(fig_ctx, ctx, axis, spec);
			pipeline_kwargs...,
		)
	end
	if display_plot
		_display!(backend_ctx, assembly.figure; title = spec.title)
	end
	return assembly
end

function _display!(backend_ctx, fig::Makie.Figure; title::AbstractString = "")
	if backend_ctx.interactive && backend_ctx.window !== nothing
		display(backend_ctx.window, fig)
		if !isempty(title) && hasproperty(backend_ctx.window, :title)
			backend_ctx.window.title[] = title
		end
	else
		BackendHandler.renderfig(fig)
	end
	return nothing
end

function _build_hist_plot!(fig_ctx, ctx, axis, spec::MCPlotHistSpec)
	axis.title = spec.title

	x_scale = 10.0 ^ spec.x_exp
	axis.ytickformat[] = vals -> TICKFORMATTER(vals)

	# Conditional axis labels
	if spec.mode == :qq
		# Q-Q plot: both axes get the exponent
		axis.xlabel = _axis_label(spec.xlabel, spec.x_exp)
		axis.ylabel = _axis_label(spec.ylabel, spec.x_exp)
	else
		# Hist/ECDF: only x-axis gets the exponent
		axis.xlabel = _axis_label(spec.xlabel, spec.x_exp)
		axis.ylabel = spec.ylabel
	end

	# --- PLOTTING LOGIC SWITCH ---
	if spec.mode == :hist

		scaled_edges = spec.bins ./ x_scale

		if spec.data == :samples || spec.data == :both
			vals_scaled = spec.values ./ x_scale
			hist!(
				axis, vals_scaled;
				bins = scaled_edges,
				normalization = spec.normalization,
				color = :steelblue, strokecolor = :white, strokewidth = 0.5,
				label = "samples",
			)

		end

		if spec.data == :pdf || spec.data == :both
			pdf = spec.pdf_obj
			dens_scaled = pdf.dens .* x_scale
			y_values = [dens_scaled; dens_scaled[end]]
			stairs!(
				axis, scaled_edges, y_values;
				step = :post, color = :red, linewidth = 2, overdraw = true,
				label = "model PDF",
			)

		end

		Makie.autolimits!(axis)
		ylims!(axis, 0, nothing) # Glue bars to x-axis

	elseif spec.mode == :ecdf
		# --- ECDF PLOT ---
		pdf = spec.pdf_obj
		vals_scaled = spec.values ./ x_scale
		ecdf_func = ecdf(vals_scaled)

		# Define plot range from scaled PDF edges
		xmin = minimum(pdf.edges) / x_scale
		xmax = maximum(pdf.edges) / x_scale
		pad = (xmax - xmin) * 0.05
		xs = range(xmin - pad, xmax + pad, length = 500)


		# 1. Plot Model CDF (Theory)
		# We must feed *physical* values (xs .* x_scale) to the cdf function
		model_cdf_data = cdf.(Ref(pdf), xs .* x_scale)
		lines!(axis, xs, model_cdf_data,
			color = :red, linewidth = 2, label = "model CDF",
		)

		# 2. Plot ECDF (Data)
		lines!(axis, xs, ecdf_func.(xs),
			color = :blue, linestyle = :dash, linewidth = 2, label = "empirical",
		)


		Makie.autolimits!(axis)
		ylims!(axis, 0, nothing) # CDFs are bounded [0, 1]

	elseif spec.mode == :qq
		# --- Q-Q PLOT ---
		vals_scaled = spec.values ./ x_scale
		sample_quantiles = sort(vals_scaled)
		n = length(sample_quantiles)
		probs = ((1:n) .- 0.5) ./ n

		pdf = spec.pdf_obj
		s = sampler(pdf) # Sampler on physically-scaled PDF

		model_quantiles_physical = quantile.(Ref(s), probs)
		model_quantiles_scaled = model_quantiles_physical ./ x_scale

		# 1. Plot the quantiles
		scatter!(axis, sample_quantiles, model_quantiles_scaled,
			color = :steelblue, markersize = 6, label = "quantiles",
		)


		# 2. Plot the y=x line
		diag_min = min(sample_quantiles[1], model_quantiles_scaled[1])
		diag_max = max(sample_quantiles[end], model_quantiles_scaled[end])
		lines!(axis, [diag_min, diag_max], [diag_min, diag_max],
			color = :black, linestyle = :dash, linewidth = 2, label = "perfect fit",
		)

		Makie.autolimits!(axis)
		# No ylims! for Q-Q
	end

	# --- Buttons (identical) ---
	buttons = [
		ControlButtonSpec(
			(_ctx, _btn) -> (Makie.autolimits!(axis); nothing), # Use autolimits
			icon = MI_REFRESH,
			on_success = ControlReaction(status_string = "Axis limits reset"),
		),
		ControlButtonSpec(
			(_ctx, _btn) -> _save_hist_export(spec, axis),
			icon = MI_SAVE,
			on_success = ControlReaction(
				status_string = path -> string("Saved SVG to ", basename(path)),
			),
		),
	]

	legend_builder =
		parent ->
			Makie.Legend(
				parent,
				axis;
				orientation = :vertical,
			)

	# --- Return (identical) ---
	return PlotBuildArtifacts(
		axis = axis,
		legends = legend_builder,
		colorbars = Any[],
		control_buttons = buttons,
		control_toggles = ControlToggleSpec[],
		status_message = nothing,
	)
end

function plot(
	obj::LineParametersMC,
	values_expr;
	ijk::Union{Nothing, NTuple{3, Int}} = nothing,
	length_unit::Symbol = :kilo,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	per_length::Bool = true,
	quantity_units = nothing,
	nbins::Union{Nothing, Int} = nothing,
	normalization::Symbol = :none,
	data::Symbol = :samples,
	mode::Symbol = :hist,
	backend = nothing,
	display_plot::Bool = true,
)
	data in (:samples, :pdf, :both) ||
		Base.error("`data` must be one of :samples, :pdf, or :both")

	mode in (:hist, :ecdf, :qq) ||
		Base.error("`mode` must be one of :hist, :ecdf, or :qq")

	specs = _hist_specs(
		obj,
		values_expr;
		ijk = ijk,
		length_unit = length_unit,
		fig_size = fig_size,
		per_length = per_length,
		quantity_units = quantity_units,
		nbins = nbins,
		normalization = normalization,
		data = data,
		mode = mode,
	)
	spec = first(specs)
	return _render_hist_spec(spec; backend = backend, display_plot = display_plot)
end


# ### Histogram methods for CableDesignMC -> CableDesignHistSpec

const _CABLE_DESIGN_SUPPORTED_QUANTITIES = (:R, :L, :C)

function _parse_cabledesign_quantity(values_expr)
	values_expr isa Symbol && return values_expr
	values_expr isa Expr &&
		values_expr.head === :quote &&
		length(values_expr.args) == 1 &&
		values_expr.args[1] isa Symbol &&
		return values_expr.args[1]
	Base.error(
		"Provide values as Symbol (:R, :L or :C) when plotting CableDesignMC samples",
	)
end

# --- Builder for CableDesign (New, refactored logic) ---
function _build_hist_spec(
	obj::CableDesignMC,
	values_expr;
	length_unit::Symbol = :kilo,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	per_length::Bool = true,
	quantity_units = nothing,
	nbins::Union{Nothing, Int} = nothing,
	normalization::Symbol = :none,
	data::Symbol = :samples,
	mode::Symbol = :hist, # <-- ADDED
)
	# MODIFIED: Force data=:both
	plot_data = (mode == :hist) ? data : :both

	vals = nothing
	pdf_obj = nothing

	values_sym = _parse_cabledesign_quantity(values_expr)
	values_sym in _CABLE_DESIGN_SUPPORTED_QUANTITIES || Base.error(
		"CableDesignMC provides samples only for :R, :L and :C; got $(values_sym)",
	)

	meta = _quantity_metadata(values_sym)
	units = normalize_quantity_units(quantity_units)
	q_prefix = resolve_quantity_prefix(meta.quantity, units)
	c_scale = (per_length ? length_scale(length_unit) : 1.0) * quantity_scale(q_prefix)

	# --- Load Data (uses plot_data) ---
	if plot_data == :samples || plot_data == :both
		obj.samples === nothing &&
			Base.error("mode=:$mode requires samples, but none are available.")
		samps = getfield(obj.samples, values_sym)
		vals = collect(samps) .* c_scale
	end
	if plot_data == :pdf || plot_data == :both
		obj.pdf === nothing &&
			Base.error("mode=:$mode requires PDF, but none is available.")
		raw_pdf = getfield(obj.pdf, values_sym)
		scaled_edges = raw_pdf.edges .* c_scale
		scaled_dens = raw_pdf.dens ./ c_scale
		pdf_obj = LineParametersPDF(scaled_edges, scaled_dens)
	end

	# --- Binning (Conditional) & Scaling ---
	local bin_edges::Vector{<:Real}
	local current_norm::Symbol
	local x_exp::Int

	if mode == :hist
		if pdf_obj !== nothing
			bin_edges = pdf_obj.edges
			current_norm = :pdf
		elseif vals !== nothing
			current_norm = normalization
			hist_nbins = isnothing(nbins) ? _auto_nbins(vals) : nbins
			h_fit = fit(Histogram, vals; nbins = hist_nbins, closed = :left)
			bin_edges = collect(h_fit.edges[1])
		else
			error("No data (samples or PDF) to plot.")
		end
	else
		bin_edges = Float64[]
		current_norm = :none
	end

	if vals !== nothing
		_, x_exp = autoscale_axis(vals)
	else
		_, x_exp = autoscale_axis(pdf_obj.edges)
	end

	# --- Labels and Title (Conditional) ---
	local title::String
	local xlabel::String
	local ylabel::String

	xlabel_unit = composite_unit(q_prefix, meta.unit.symbol, per_length, length_unit)
	base_xlabel = string(meta.axis_label, " [", xlabel_unit, "]")

	if mode == :hist
		title = string(meta.title, " histogram (base values)")
		xlabel = base_xlabel
		ylabel =
			current_norm == :none ? "count" :
			(current_norm == :pdf ? "density" : String(current_norm))
	elseif mode == :ecdf
		title = string(meta.title, " CDF (base values)")
		xlabel = base_xlabel
		ylabel = "cumulative probability"
	elseif mode == :qq
		title = string(meta.title, " Q-Q plot (base values)")
		xlabel = "sampled quantiles"
		ylabel = "model quantiles"
	end

	return MCPlotHistSpec(
		values_sym, meta.symbol, title, xlabel, ylabel,
		vals, pdf_obj, bin_edges, x_exp,
		fig_size, current_norm, plot_data, mode,
	)
end

function _hist_specs(
	obj::CableDesignMC,
	values_expr;
	length_unit::Symbol = :kilo,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	per_length::Bool = true,
	quantity_units = nothing,
	nbins::Union{Nothing, Int} = nothing,
	normalization::Symbol = :none,
	data::Symbol = :samples,
	mode::Symbol = :hist,
)
	spec = _build_hist_spec( # Calls CD-MC builder
		obj,
		values_expr;
		length_unit = length_unit,
		fig_size = fig_size,
		per_length = per_length,
		quantity_units = quantity_units,
		nbins = nbins,
		normalization = normalization,
		data = data,
		mode = mode,
	)
	return [spec]
end

# --- hist for CableDesignMC ---
function plot(
	obj::CableDesignMC,
	values_expr;
	length_unit::Symbol = :kilo,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	quantity_units = nothing,
	nbins::Union{Nothing, Int} = nothing, # <-- Now Union{Nothing, Int}
	normalization::Symbol = :none,
	data::Symbol = :samples, # <-- New arg
	mode::Symbol = :hist,
	backend = nothing,
	display_plot::Bool = true,
)
	data in (:samples, :pdf, :both) ||
		Base.error("`data` must be one of :samples, :pdf, or :both")

	mode in (:hist, :ecdf, :qq) ||
		Base.error("`mode` must be one of :hist, :ecdf, or :qq")

	specs = _hist_specs( # Dispatches to CableDesignMC version
		obj,
		values_expr;
		length_unit = length_unit,
		fig_size = fig_size,
		per_length = true,
		quantity_units = quantity_units,
		nbins = nbins,
		normalization = normalization,
		data = data,
		mode = mode,
	)
	spec = first(specs)
	return _render_hist_spec(spec; backend = backend, display_plot = display_plot)
end

