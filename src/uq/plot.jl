using Makie
import Makie: hist

using Base: basename
using Dates
using Printf: @sprintf

using ..Engine:
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
	get_symbol,
	get_unit_symbol,
	ComponentMetadata

import ..BackendHandler: BackendHandler, next_fignum

using ..PlotUIComponents:
	PlotAssembly,
	PlotBuildArtifacts,
	ControlButtonSpec,
	ControlReaction,
	ControlToggleSpec,
	TICKFORMATTER,
	MI_REFRESH,
	MI_SAVE,
	ICON_TTF,
	_make_window,
	_run_plot_pipeline,
	with_plot_theme,
	ensure_export_background!

const _HIST_EXTENSION = "svg"

struct LineParametersHistSpec
	quantity::Symbol
	symbol::String
	title::String
	xlabel::String
	ylabel::String
	values::Vector{<:Real}
	nbins::Int
	x_exp::Int
	fig_size::Union{Nothing, Tuple{Int, Int}}
	freq::Real
	normalization::Symbol
end

function _mc_quantity_metadata()
	sdesc = get_description(SeriesImpedance(zeros(1, 1, 1)))
	symb = get_symbol(SeriesImpedance(zeros(1, 1, 1)))
	sunit = get_unit_symbol(SeriesImpedance(zeros(1, 1, 1)))

	adesc = get_description(ShuntAdmittance(zeros(1, 1, 1)))
	asymb = get_symbol(ShuntAdmittance(zeros(1, 1, 1)))
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
	obj::LineParametersMCSummary,
	values_expr;
	ijk::Union{Nothing, NTuple{3, Int}} = nothing,
	length_unit::Symbol = :kilo,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	per_length::Bool = true,
	quantity_units = nothing,
	nbins::Int = 15,
	normalization::Symbol = :none,
)
	obj.samples === nothing &&
		Base.error("Sampled values are not available in this LineParametersMCSummary")

	values_sym, _ijk = _parse_values_ref(values_expr, ijk)
	meta = _quantity_metadata(values_sym)
	units = normalize_quantity_units(quantity_units)
	q_prefix = resolve_quantity_prefix(meta.quantity, units)

	l_scale = per_length ? length_scale(length_unit) : 1.0
	y_scale = quantity_scale(q_prefix)

	i, j, k = _ijk
	samps = getfield(obj.samples, values_sym)
	max_i, max_j, max_k, _ = size(samps)
	(1 <= i <= max_i && 1 <= j <= max_j && 1 <= k <= max_k) || Base.error(
		"indices (i=$(i), j=$(j), k=$(k)) out of bounds for samples size $(size(samps))",
	)

	raw_vals = @view samps[i, j, k, :]
	vals = collect(raw_vals) .* y_scale .* l_scale

	_, x_exp = autoscale_axis(vals)

	xlabel_unit = composite_unit(q_prefix, meta.unit.symbol, per_length, length_unit)
	xlabel = string(meta.axis_label, " [", xlabel_unit, "]")
	ylabel =
		normalization == :none ?
		"count" :
		String(normalization)


	freq_val = obj.f[k]
	# freq_str = @sprintf("%.4g", freq_val)
	title = string(meta.title, " histogram @ f=", @sprintf("%.4g", freq_val), " Hz")

	return LineParametersHistSpec(
		values_sym,
		meta.symbol,
		title,
		xlabel,
		ylabel,
		vals,
		nbins,
		x_exp,
		fig_size,
		freq_val,
		normalization,
	)
end

function lineparametermc_hist_specs(
	obj::LineParametersMCSummary,
	values_expr;
	ijk::Union{Nothing, NTuple{3, Int}} = nothing,
	length_unit::Symbol = :kilo,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	per_length::Bool = true,
	quantity_units = nothing,
	nbins::Int = 15,
	normalization::Symbol = :none,
)
	spec = _build_hist_spec(
		obj,
		values_expr;
		ijk = ijk,
		length_unit = length_unit,
		fig_size = fig_size,
		per_length = per_length,
		quantity_units = quantity_units,
		nbins = nbins,
		normalization = normalization,
	)
	return [spec]
end

function _default_export_path_hist(spec::LineParametersHistSpec)
	base_title = strip(spec.title)
	name = strip(replace(lowercase(base_title), r"[^0-9a-z]+" => "_"), '_')
	isempty(name) && (name = "lineparameters_hist")
	timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
	filename = string(name, "_", timestamp, ".", _HIST_EXTENSION)
	return joinpath(pwd(), filename)
end

function _save_hist_export(spec::LineParametersHistSpec, axis)
	fig = build_hist_export(spec)
	trim!(fig.layout)
	path = _default_export_path_hist(spec)
	Makie.save(path, fig)
	return path
end

function build_hist_export(spec::LineParametersHistSpec)
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
	spec::LineParametersHistSpec;
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

function _build_hist_plot!(fig_ctx, ctx, axis, spec::LineParametersHistSpec)
	axis.title  = spec.title
	axis.xlabel = _axis_label(spec.xlabel, spec.x_exp)
	axis.ylabel = spec.ylabel

	x_scale = 10.0 ^ spec.x_exp

	axis.ytickformat[] = vals -> TICKFORMATTER(vals)

	vals_scaled = spec.values ./ x_scale

	hist!(
		axis,
		vals_scaled;
		bins = spec.nbins,
		normalization = spec.normalization,
		color = :steelblue,
		strokecolor = :white,
		strokewidth = 0.5,
	)

	Makie.autolimits!(axis)

	buttons = [
		ControlButtonSpec(
			(_ctx, _btn) -> (Makie.reset_limits!(axis); nothing);
			icon = MI_REFRESH,
			on_success = ControlReaction(status_string = "Axis limits reset"),
		),
		ControlButtonSpec(
			(_ctx, _btn) -> _save_hist_export(spec, axis);
			icon = MI_SAVE,
			on_success = ControlReaction(
				status_string = path -> string("Saved SVG to ", basename(path)),
			),
		),
	]

	return PlotBuildArtifacts(
		axis = axis,
		legends = parent ->
			Makie.Legend(
				parent,
				[
					Makie.PolyElement(
						color = :steelblue,
						strokecolor = :white,
						strokewidth = 0.5,
					),
				],
				["MC samples"];
				orientation = :vertical,
			),
		colorbars = Any[],
		control_buttons = buttons,
		control_toggles = ControlToggleSpec[],
		status_message = nothing,
	)
end

function hist(
	obj::LineParametersMCSummary,
	values_expr;
	ijk::Union{Nothing, NTuple{3, Int}} = nothing,
	length_unit::Symbol = :kilo,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	per_length::Bool = true,
	quantity_units = nothing,
	nbins::Int = 15,
	normalization::Symbol = :none,
	backend = nothing,
	display_plot::Bool = true,
)
	specs = lineparametermc_hist_specs(
		obj,
		values_expr;
		ijk = ijk,
		length_unit = length_unit,
		fig_size = fig_size,
		per_length = per_length,
		quantity_units = quantity_units,
		nbins = nbins,
		normalization = normalization,
	)
	spec = first(specs)
	return _render_hist_spec(spec; backend = backend, display_plot = display_plot)
end


### Histogram methods for CableDesignMCSummary -> CableDesignHistSpec

struct CableDesignHistSpec
	quantity::Symbol
	symbol::String
	title::String
	xlabel::String
	ylabel::String
	values::Vector{<:Real}
	nbins::Int
	x_exp::Int
	fig_size::Union{Nothing, Tuple{Int, Int}}
	normalization::Symbol
end

const _CABLE_DESIGN_SUPPORTED_QUANTITIES = (:R, :L, :C)

function _parse_cabledesign_quantity(values_expr)
	values_expr isa Symbol && return values_expr
	values_expr isa Expr &&
		values_expr.head === :quote &&
		length(values_expr.args) == 1 &&
		values_expr.args[1] isa Symbol &&
		return values_expr.args[1]
	Base.error(
		"Provide values as Symbol (:R, :L or :C) when plotting CableDesignMCSummary samples",
	)
end

function _build_cabledesign_hist_spec(
	obj::CableDesignMCSummary,
	values_expr;
	length_unit::Symbol = :kilo,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	per_length::Bool = true,
	quantity_units = nothing,
	nbins::Int = 15,
	normalization::Symbol = :none,
)
	obj.samples === nothing &&
		Base.error("Sampled values are not available in this CableDesignMCSummary")

	values_sym = _parse_cabledesign_quantity(values_expr)
	values_sym in _CABLE_DESIGN_SUPPORTED_QUANTITIES || Base.error(
		"CableDesignMCSummary provides samples only for :R, :L and :C; got $(values_sym)",
	)

	meta = _quantity_metadata(values_sym)
	units = normalize_quantity_units(quantity_units)
	q_prefix = resolve_quantity_prefix(meta.quantity, units)

	l_scale = per_length ? length_scale(length_unit) : 1.0
	y_scale = quantity_scale(q_prefix)

	samps = getfield(obj.samples, values_sym)
	vals = collect(samps) .* y_scale .* l_scale

	_, x_exp = autoscale_axis(vals)

	xlabel_unit = composite_unit(q_prefix, meta.unit.symbol, per_length, length_unit)
	xlabel = string(meta.axis_label, " [", xlabel_unit, "]")
	ylabel =
		normalization == :none ?
		"count" :
		String(normalization)

	title = string(meta.title, " histogram (base values)")

	return CableDesignHistSpec(
		values_sym,
		meta.symbol,
		title,
		xlabel,
		ylabel,
		vals,
		nbins,
		x_exp,
		fig_size,
		normalization,
	)
end

function cabledesign_hist_specs(
	obj::CableDesignMCSummary,
	values_expr;
	length_unit::Symbol = :kilo,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	per_length::Bool = true,
	quantity_units = nothing,
	nbins::Int = 15,
	normalization::Symbol = :none,
)
	spec = _build_cabledesign_hist_spec(
		obj,
		values_expr;
		length_unit = length_unit,
		fig_size = fig_size,
		per_length = per_length,
		quantity_units = quantity_units,
		nbins = nbins,
		normalization = normalization,
	)
	return [spec]
end

function _default_export_path_hist(spec::CableDesignHistSpec)
	base_title = strip(spec.title)
	name = strip(replace(lowercase(base_title), r"[^0-9a-z]+" => "_"), '_')
	isempty(name) && (name = "cabledesign_hist")
	timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
	filename = string(name, "_", timestamp, ".", _HIST_EXTENSION)
	return joinpath(pwd(), filename)
end

function _save_hist_export(spec::CableDesignHistSpec, axis)
	fig = build_hist_export(spec)
	trim!(fig.layout)
	path = _default_export_path_hist(spec)
	Makie.save(path, fig)
	return path
end

function build_hist_export(spec::CableDesignHistSpec)
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
	spec::CableDesignHistSpec;
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
			(fig_ctx, ctx, axis) -> _build_hist_plot!(fig_ctx, ctx, axis, spec);
			pipeline_kwargs...,
		)
	end
	if display_plot
		_display!(backend_ctx, assembly.figure; title = spec.title)
	end
	return assembly
end

function _build_hist_plot!(fig_ctx, ctx, axis, spec::CableDesignHistSpec)
	axis.title  = spec.title
	axis.xlabel = _axis_label(spec.xlabel, spec.x_exp)
	axis.ylabel = spec.ylabel

	x_scale = 10.0 ^ spec.x_exp

	axis.ytickformat[] = vals -> TICKFORMATTER(vals)

	vals_scaled = spec.values ./ x_scale

	hist!(
		axis,
		vals_scaled;
		bins = spec.nbins,
		normalization = spec.normalization,
		color = :steelblue,
		strokecolor = :white,
		strokewidth = 0.5,
	)

	Makie.autolimits!(axis)

	buttons = [
		ControlButtonSpec(
			(_ctx, _btn) -> (Makie.reset_limits!(axis); nothing);
			icon = MI_REFRESH,
			on_success = ControlReaction(status_string = "Axis limits reset"),
		),
		ControlButtonSpec(
			(_ctx, _btn) -> _save_hist_export(spec, axis);
			icon = MI_SAVE,
			on_success = ControlReaction(
				status_string = path -> string("Saved SVG to ", basename(path)),
			),
		),
	]

	return PlotBuildArtifacts(
		axis = axis,
		legends = parent ->
			Makie.Legend(
				parent,
				[
					Makie.PolyElement(
						color = :steelblue,
						strokecolor = :white,
						strokewidth = 0.5,
					),
				],
				["MC samples"];
				orientation = :vertical,
			),
		colorbars = Any[],
		control_buttons = buttons,
		control_toggles = ControlToggleSpec[],
		status_message = nothing,
	)
end

function hist(
	obj::CableDesignMCSummary,
	values_expr;
	length_unit::Symbol = :kilo,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	# per_length::Bool = true,
	quantity_units = nothing,
	nbins::Int = 15,
	normalization::Symbol = :none,
	backend = nothing,
	display_plot::Bool = true,
)
	specs = cabledesign_hist_specs(
		obj,
		values_expr;
		length_unit = length_unit,
		fig_size = fig_size,
		per_length = true,
		quantity_units = quantity_units,
		nbins = nbins,
		normalization = normalization,
	)
	spec = first(specs)
	return _render_hist_spec(spec; backend = backend, display_plot = display_plot)
end
