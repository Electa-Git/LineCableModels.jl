
using Makie
import Makie: plot

include("../plotbuilder/plothelpers.jl")

using Measurements: Measurements

const _ICON_FN =
	(icon; text = nothing, kwargs...) ->
		with_icon(icon; text = text === nothing ? "" : text, kwargs...)

const LP_FIG_SIZE = (800, 400)

const METRIC_PREFIX_EXPONENT = Dict(
	:yocto => -24,
	:zepto => -21,
	:atto => -18,
	:femto => -15,
	:pico => -12,
	:nano => -9,
	:micro => -6,
	:milli => -3,
	:centi => -2,
	:deci => -1,
	:base => 0,
	:deca => 1,
	:hecto => 2,
	:kilo => 3,
	:mega => 6,
	:giga => 9,
	:tera => 12,
	:peta => 15,
	:exa => 18,
	:zetta => 21,
	:yotta => 24,
)

const METRIC_PREFIX_SYMBOL = Dict(
	:yocto => "y",
	:zepto => "z",
	:atto => "a",
	:femto => "f",
	:pico => "p",
	:nano => "n",
	:micro => "μ",
	:milli => "m",
	:centi => "c",
	:deci => "d",
	:base => "",
	:deca => "da",
	:hecto => "h",
	:kilo => "k",
	:mega => "M",
	:giga => "G",
	:tera => "T",
	:peta => "P",
	:exa => "E",
	:zetta => "Z",
	:yotta => "Y",
)

const DEFAULT_QUANTITY_UNITS = Dict(
	:impedance => :base,
	:admittance => :base,
	:resistance => :base,
	:inductance => :milli,
	:conductance => :base,
	:capacitance => :micro,
	:angle => :base,
)

struct UnitSpec
	symbol::String
	per_length::Bool
end

struct ComponentMetadata
	component::Symbol
	quantity::Symbol
	symbol::String
	title::String
	axis_label::String
	unit::UnitSpec
end

struct LineParametersPlotSpec <: AbstractPlotSpec
	parent_kind::Symbol
	component::Symbol
	symbol::String
	title::String
	xlabel::String
	ylabel::String
	freqs::Vector{<:Real}
	raw_freqs::Vector{<:Real}
	curves::Vector{Vector{<:Real}}
	raw_curves::Vector{Vector{<:Real}}
	labels::Vector{String}
	x_exp::Int
	y_exp::Int
	fig_size::Union{Nothing, Tuple{Int, Int}}
	xscale::Base.RefValue{Function}
	yscale::Base.RefValue{Function}
end

get_description(::SeriesImpedance) = (
	impedance = "Series impedance",
	resistance = "Series resistance",
	inductance = "Series inductance",
)

get_symbol(::SeriesImpedance) = (
	impedance = "Z",
	resistance = "R",
	inductance = "L",
)

get_unit_symbol(::SeriesImpedance) = (
	impedance = "Ω",
	resistance = "Ω",
	inductance = "H",
)

get_description(::ShuntAdmittance) = (
	admittance = "Shunt admittance",
	conductance = "Shunt conductance",
	capacitance = "Shunt capacitance",
)

get_symbol(::ShuntAdmittance) = (
	admittance = "Y",
	conductance = "G",
	capacitance = "C",
)

get_unit_symbol(::ShuntAdmittance) = (
	admittance = "S",
	conductance = "S",
	capacitance = "F",
)

parent_kind(::SeriesImpedance) = :series_impedance
parent_kind(::ShuntAdmittance) = :shunt_admittance

metric_exponent(prefix::Symbol) =
	get(METRIC_PREFIX_EXPONENT, prefix) do
		Base.error("Unsupported metric prefix :$(prefix)")
	end

prefix_symbol(prefix::Symbol) =
	get(METRIC_PREFIX_SYMBOL, prefix) do
		Base.error("Unsupported metric prefix :$(prefix)")
	end


quantity_scale(prefix::Symbol) = 10.0 ^ (-metric_exponent(prefix))
length_scale(prefix::Symbol) = 10.0 ^ (metric_exponent(prefix))
frequency_scale(prefix::Symbol) = quantity_scale(prefix)

function unit_text(quantity_prefix::Symbol, base_unit::String)
	ps = prefix_symbol(quantity_prefix)
	return isempty(ps) ? base_unit : string(ps, base_unit)
end

function length_unit_text(prefix::Symbol)
	ps = prefix_symbol(prefix)
	return isempty(ps) ? "m" : string(ps, "m")
end

function composite_unit(
	quantity_prefix::Symbol,
	base_unit::String,
	per_length::Bool,
	length_prefix::Symbol,
)
	numerator = unit_text(quantity_prefix, base_unit)
	if per_length
		denominator = length_unit_text(length_prefix)
		return string(numerator, "/", denominator)
	else
		return numerator
	end
end

function frequency_axis_label(prefix::Symbol)
	unit = unit_text(prefix, "Hz")
	return string("frequency [", unit, "]")
end

function normalize_quantity_units(units)
	table = Dict(DEFAULT_QUANTITY_UNITS)
	if units isa Symbol
		for key in keys(table)
			table[key] = units
		end
	elseif units isa NamedTuple
		for (key, val) in pairs(units)
			table[key] = val
		end
	elseif units isa AbstractDict
		for (key, val) in units
			table[key] = val
		end
	elseif units === nothing
		return table
	else
		Base.error("Unsupported quantity unit specification $(typeof(units))")
	end
	return table
end

function resolve_quantity_prefix(quantity::Symbol, units::AbstractDict{Symbol, Symbol})
	return get(units, quantity, get(DEFAULT_QUANTITY_UNITS, quantity, :base))
end

function resolve_conductors(data_dims::NTuple{3, Int}, con)
	nrows, ncols, _ = data_dims
	if con === nothing
		return collect(1:nrows), collect(1:ncols)
	elseif con isa Tuple && length(con) == 2
		isel = collect_indices(con[1], nrows)
		jsel = collect_indices(con[2], ncols)
		return isel, jsel
	else
		Base.error("Conductor selector must be a tuple (i_sel, j_sel)")
	end
end

function collect_indices(sel, n)
	if sel === nothing
		return collect(1:n)
	elseif sel isa Integer
		(1 <= sel <= n) ||
			Base.error("Index $(sel) out of bounds for dimension of size $(n)")
		return [sel]
	elseif sel isa AbstractVector
		indices = collect(Int, sel)
		for idx in indices
			(1 <= idx <= n) ||
				Base.error("Index $(idx) out of bounds for dimension of size $(n)")
		end
		return indices
	elseif sel isa AbstractRange
		indices = collect(sel)
		for idx in indices
			(1 <= idx <= n) ||
				Base.error("Index $(idx) out of bounds for dimension of size $(n)")
		end
		return indices
	elseif sel isa Colon
		return collect(1:n)
	else
		Base.error("Unsupported selector $(sel)")
	end
end

function components_for(
	obj::SeriesImpedance,
	mode::Symbol,
	coord::Symbol;
	per_length::Bool = true,
)
	desc = get_description(obj)
	sym = get_symbol(obj)
	units = get_unit_symbol(obj)
	if mode == :ZY
		coord in (:cart, :polar) || Base.error("Unsupported coordinate system $(coord)")
		if coord == :cart
			return ComponentMetadata[
				ComponentMetadata(:real, :impedance, sym.impedance,
					string(desc.impedance, " – real part"),
					string("real(", sym.impedance, ")"),
					UnitSpec(units.impedance, per_length)),
				ComponentMetadata(:imag, :impedance, sym.impedance,
					string(desc.impedance, " – imaginary part"),
					string("imag(", sym.impedance, ")"),
					UnitSpec(units.impedance, per_length)),
			]
		else
			return ComponentMetadata[
				ComponentMetadata(:magnitude, :impedance, sym.impedance,
					string(desc.impedance, " – magnitude"),
					string("|", sym.impedance, "|"),
					UnitSpec(units.impedance, per_length)),
				ComponentMetadata(:angle, :angle, sym.impedance,
					string(desc.impedance, " – angle"),
					string("angle(", sym.impedance, ")"),
					UnitSpec("deg", false)),
			]
		end
	elseif mode == :RLCG
		return ComponentMetadata[
			ComponentMetadata(:resistance, :resistance, sym.resistance,
				desc.resistance,
				sym.resistance,
				UnitSpec(units.resistance, per_length)),
			ComponentMetadata(:inductance, :inductance, sym.inductance,
				desc.inductance,
				sym.inductance,
				UnitSpec(units.inductance, per_length)),
		]
	else
		Base.error("Unsupported mode $(mode)")
	end
end

function components_for(
	obj::ShuntAdmittance,
	mode::Symbol,
	coord::Symbol;
	per_length::Bool = true,
)
	desc = get_description(obj)
	sym = get_symbol(obj)
	units = get_unit_symbol(obj)
	if mode == :ZY
		coord in (:cart, :polar) || Base.error("Unsupported coordinate system $(coord)")
		if coord == :cart
			return ComponentMetadata[
				ComponentMetadata(:real, :admittance, sym.admittance,
					string(desc.admittance, " – real part"),
					string("real(", sym.admittance, ")"),
					UnitSpec(units.admittance, per_length)),
				ComponentMetadata(:imag, :admittance, sym.admittance,
					string(desc.admittance, " – imaginary part"),
					string("imag(", sym.admittance, ")"),
					UnitSpec(units.admittance, per_length)),
			]
		else
			return ComponentMetadata[
				ComponentMetadata(:magnitude, :admittance, sym.admittance,
					string(desc.admittance, " – magnitude"),
					string("|", sym.admittance, "|"),
					UnitSpec(units.admittance, per_length)),
				ComponentMetadata(:angle, :angle, sym.admittance,
					string(desc.admittance, " – angle"),
					string("angle(", sym.admittance, ")"),
					UnitSpec("deg", false)),
			]
		end
	elseif mode == :RLCG
		if (coord == :cart || coord == :polar)
			@warn "Ignoring argument :$(coord) for RLCG parameters"
		end
		return ComponentMetadata[
			ComponentMetadata(:conductance, :conductance, sym.conductance,
				desc.conductance,
				sym.conductance,
				UnitSpec(units.conductance, per_length)),
			ComponentMetadata(:capacitance, :capacitance, sym.capacitance,
				desc.capacitance,
				sym.capacitance,
				UnitSpec(units.capacitance, per_length)),
		]
	else
		Base.error("Unsupported mode $(mode)")
	end
end

function component_values(component::Symbol, slice, freqs::Vector{<:Real})
	data = collect(slice)
	if component === :real
		return (real.(data))
	elseif component === :imag
		return (imag.(data))
	elseif component === :magnitude
		return (abs.(data))
	elseif component === :angle
		return rad2deg.((angle.(data)))
	elseif component === :resistance || component === :conductance
		return (real.(data))
	elseif component === :inductance
		imag_part = (imag.(data))
		return reactance_to_l(imag_part, freqs)
	elseif component === :capacitance
		imag_part = (imag.(data))
		return reactance_to_c(imag_part, freqs)
	else
		Base.error("Unsupported component $(component)")
	end
end

function reactance_to_l(imag_part::Vector{<:Real}, freqs::Vector{<:Real})
	result = similar(freqs, promote_type(eltype(imag_part), eltype(freqs)))
	two_pi = 2π
	for idx in eachindex(freqs)
		f = freqs[idx]
		if iszero(f)
			result[idx] = NaN
		else
			result[idx] = imag_part[idx] / (two_pi * f)
		end
	end
	return result
end

function reactance_to_c(imag_part::Vector{<:Real}, freqs::Vector{<:Real})
	result = similar(freqs, promote_type(eltype(imag_part), eltype(freqs)))
	two_pi = 2π
	for idx in eachindex(freqs)
		f = freqs[idx]
		if iszero(f)
			result[idx] = NaN
		else
			result[idx] = imag_part[idx] / (two_pi * f)
		end
	end
	return result
end

function legend_label(symbol::String, i::Int, j::Int)
	return string(symbol, "(", i, ",", j, ")")
end

function _axis_label(base::AbstractString, exp::Int)
	exp == 0 && return base
	return Makie.rich(
		base,
		Makie.rich("  × 10"; font = :regular, fontsize = AXIS_LABEL_FONT_SIZE),
		Makie.rich(
			superscript(string(exp));
			font = :regular,
			fontsize = AXIS_LABEL_FONT_SIZE - 2,
			# baseline_shift = 0.6,
		),
	)
end

# Return scaled data and the exponent factored out for the axis badge.
function autoscale_axis(values::AbstractVector{<:Real}; _threshold = 1e4)
	isempty(values) && return values, 0
	maxval = 0.0
	has_value = false
	for val in values
		if isnan(val)
			continue
		end
		absval = abs(val)
		if !has_value || absval > maxval
			maxval = absval
			has_value = true
		end
	end
	!has_value && return values, 0
	exp = floor(Int, log10(maxval))
	abs(exp) < 3 && return values, 0
	scale = 10.0 ^ exp
	# return values ./ scale, exp
	return values ./ scale, exp
end

function autoscale_axis_stacked(
	curves::AbstractVector{<:AbstractVector{<:Real}};
	_threshold = 1e4,
)
	isempty(curves) && return curves, 0
	maxval = 0.0
	has_value = false
	for curve in curves
		for val in curve
			if isnan(val)
				continue
			end
			absval = abs(val)
			if !has_value || absval > maxval
				maxval = absval
				has_value = true
			end
		end
	end
	!has_value && return curves, 0
	exp = floor(Int, log10(maxval))
	abs(exp) < 3 && return curves, 0
	scale = 10.0 ^ exp
	scaled_curves = [curve ./ scale for curve in curves]
	return scaled_curves, exp
end

function lineparameter_plot_specs(
	obj::SeriesImpedance,
	freqs::AbstractVector;
	mode::Symbol = :ZY,
	coord::Symbol = :cart,
	freq_unit::Symbol = :base,
	length_unit::Symbol = :base,
	quantity_units = nothing,
	con = nothing,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	xscale::Function = Makie.identity,
	yscale::Function = Makie.identity,
	per_length::Bool = true,
)
	freq_vec = collect(freqs)
	nfreq = length(freq_vec)
	if nfreq <= 1
		@warn "Frequency vector has $(nfreq) sample(s); nothing to plot."
		return LineParametersPlotSpec[]
	end
	size(obj.values, 3) == nfreq ||
		Base.error("Frequency vector length does not match impedance samples")
	comps = components_for(obj, mode, coord; per_length = per_length)
	units = normalize_quantity_units(quantity_units)
	freq_scale = frequency_scale(freq_unit)
	raw_freq_axis = freq_vec .* freq_scale
	freq_axis, freq_exp = autoscale_axis(raw_freq_axis)
	xlabel_base = frequency_axis_label(freq_unit)
	(isel, jsel) = resolve_conductors(size(obj.values), con)
	specs = LineParametersPlotSpec[]
	for meta in comps
		q_prefix = resolve_quantity_prefix(meta.quantity, units)
		y_scale = quantity_scale(q_prefix)
		l_scale = meta.unit.per_length ? length_scale(length_unit) : 1.0
		ylabel_unit =
			composite_unit(q_prefix, meta.unit.symbol, meta.unit.per_length, length_unit)
		ylabel_base = string(meta.axis_label, " [", ylabel_unit, "]")

		# collect raw curves and labels
		raw_curves = Vector{Vector{<:Real}}()
		labels = String[]
		for i in isel, j in jsel
			slice = @view obj.values[i, j, :]
			raw_vals = component_values(meta.component, slice, freq_vec)
			push!(raw_curves, (raw_vals .* y_scale .* l_scale))
			push!(labels, legend_label(meta.symbol, i, j))
		end
		curves, y_exp = autoscale_axis_stacked(raw_curves)
		push!(
			specs,
			LineParametersPlotSpec(
				parent_kind(obj),
				meta.component,
				meta.symbol,
				meta.title,
				xlabel_base,
				ylabel_base,
				freq_axis,
				raw_freq_axis,
				curves,
				raw_curves,
				labels,
				freq_exp,
				y_exp,
				fig_size,
				Ref{Function}(xscale),
				Ref{Function}(yscale),
			),
		)
	end
	return specs
end

function lineparameter_plot_specs(
	obj::ShuntAdmittance,
	freqs::AbstractVector;
	mode::Symbol = :ZY,
	coord::Symbol = :cart,
	freq_unit::Symbol = :base,
	length_unit::Symbol = :base,
	quantity_units = nothing,
	con = nothing,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	xscale::Function = Makie.identity,
	yscale::Function = Makie.identity,
	per_length::Bool = true,
)
	freq_vec = collect(freqs)
	nfreq = length(freq_vec)
	if nfreq <= 1
		@warn "Frequency vector has $(nfreq) sample(s); nothing to plot."
		return LineParametersPlotSpec[]
	end
	size(obj.values, 3) == nfreq ||
		Base.error("Frequency vector length does not match admittance samples")
	comps = components_for(obj, mode, coord; per_length = per_length)
	units = normalize_quantity_units(quantity_units)
	freq_scale = frequency_scale(freq_unit)
	raw_freq_axis = freq_vec .* freq_scale
	freq_axis, freq_exp = autoscale_axis(raw_freq_axis)
	xlabel_base = frequency_axis_label(freq_unit)
	(isel, jsel) = resolve_conductors(size(obj.values), con)
	specs = LineParametersPlotSpec[]
	for meta in comps
		q_prefix = resolve_quantity_prefix(meta.quantity, units)
		y_scale = quantity_scale(q_prefix)
		l_scale = meta.unit.per_length ? length_scale(length_unit) : 1.0
		ylabel_unit =
			composite_unit(q_prefix, meta.unit.symbol, meta.unit.per_length, length_unit)
		ylabel_base = string(meta.axis_label, " [", ylabel_unit, "]")

		raw_curves = Vector{Vector{<:Real}}()
		labels = String[]
		for i in isel, j in jsel
			slice = @view obj.values[i, j, :]
			raw_vals = component_values(meta.component, slice, freq_vec)
			push!(raw_curves, (raw_vals .* y_scale .* l_scale))
			push!(labels, legend_label(meta.symbol, i, j))
		end

		curves, y_exp = autoscale_axis_stacked(raw_curves)
		push!(
			specs,
			LineParametersPlotSpec(
				parent_kind(obj),
				meta.component,
				meta.symbol,
				meta.title,
				xlabel_base,
				ylabel_base,
				freq_axis,
				raw_freq_axis,
				curves,
				raw_curves,
				labels,
				freq_exp,
				y_exp,
				fig_size,
				Ref{Function}(xscale),
				Ref{Function}(yscale),
			),
		)
	end
	return specs
end

function lineparameter_plot_specs(
	lp::LineParameters;
	mode::Symbol = :ZY,
	coord::Symbol = :cart,
	freq_unit::Symbol = :base,
	length_unit::Symbol = :base,
	quantity_units = nothing,
	con = nothing,
	fig_size::Union{Nothing, Tuple{Int, Int}} = LP_FIG_SIZE,
	xscale::Function = Makie.identity,
	yscale::Function = Makie.identity,
	per_length::Bool = true,
)
	specs = LineParametersPlotSpec[]
	append!(
		specs,
		lineparameter_plot_specs(lp.Z, lp.f;
			mode = mode,
			coord = coord,
			freq_unit = freq_unit,
			length_unit = length_unit,
			quantity_units = quantity_units,
			con = con,
			fig_size = fig_size,
			xscale = xscale,
			yscale = yscale,
			per_length = per_length,
		),
	)
	append!(
		specs,
		lineparameter_plot_specs(lp.Y, lp.f;
			mode = mode,
			coord = coord,
			freq_unit = freq_unit,
			length_unit = length_unit,
			quantity_units = quantity_units,
			con = con,
			fig_size = fig_size,
			xscale = xscale,
			yscale = yscale,
			per_length = per_length,
		),
	)
	return specs
end

function render_plot_specs(
	specs::Vector{LineParametersPlotSpec};
	backend = nothing,
	display_plot::Bool = true,
)
	assemblies = Dict{Tuple{Symbol, Symbol}, PlotAssembly}()
	for spec in specs
		assembly = _render_spec(spec; backend = backend, display_plot = display_plot)
		assemblies[(spec.parent_kind, spec.component)] = assembly
	end
	return assemblies
end

function plot(
	obj::SeriesImpedance,
	freqs::AbstractVector;
	backend = nothing,
	display_plot::Bool = true,
	per_length::Bool = true,
	kwargs...,
)
	specs = lineparameter_plot_specs(obj, freqs; per_length = per_length, kwargs...)
	return render_plot_specs(specs; backend = backend, display_plot = display_plot)
end

function plot(
	obj::ShuntAdmittance,
	freqs::AbstractVector;
	backend = nothing,
	display_plot::Bool = true,
	per_length::Bool = true,
	kwargs...,
)
	specs = lineparameter_plot_specs(obj, freqs; per_length = per_length, kwargs...)
	return render_plot_specs(specs; backend = backend, display_plot = display_plot)
end

function plot(
	lp::LineParameters;
	backend = nothing,
	display_plot::Bool = true,
	per_length::Bool = true,
	kwargs...,
)
	specs = lineparameter_plot_specs(lp; per_length = per_length, kwargs...)
	return render_plot_specs(specs; backend = backend, display_plot = display_plot)
end

function build_export_figure(spec::LineParametersPlotSpec)
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
			(fig_ctx, ctx, axis) -> _build_plot!(fig_ctx, ctx, axis, spec);
			pipeline_kwargs...,
		)
	end
	ensure_export_background!(assembly.figure)
	return assembly.figure
end

function build_export_figure(
	obj,
	key::Tuple{Symbol, Symbol};
	kwargs...,
)
	specs =
		obj isa LineParametersPlotSpec ? [obj] : lineparameter_plot_specs(obj; kwargs...)
	idx = findfirst(s -> (s.parent_kind, s.component) == key, specs)
	idx === nothing && Base.error("No plot specification found for key $(key)")
	return build_export_figure(specs[idx])
end

function _render_spec(
	spec::LineParametersPlotSpec;
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
			(fig_ctx, ctx, axis) -> _build_plot!(fig_ctx, ctx, axis, spec);
			pipeline_kwargs...,
		)
	end
	if display_plot
		_display!(backend_ctx, assembly.figure; title = spec.title)
	end
	return assembly
end

function _get_axis_data(
	raw_data::Vector{<:Real},
	scaled_data::Vector{<:Real},
	scale_func::Function,
)
	data = scale_func == Makie.log10 ? raw_data : scaled_data
	values = float(Measurements.value.(data))
	errors = if eltype(data) <: Measurements.Measurement
		float(Measurements.uncertainty.(data))
	else
		nothing
	end
	return (; values, errors)
end

function _get_axis_label(base_label::String, exponent::Int, scale_func::Function)
	if scale_func == Makie.log10
		return base_label
	else
		return _axis_label(base_label, exponent)
	end
end

function _build_plot!(fig_ctx, ctx, axis, spec::LineParametersPlotSpec)
	# ---- Axis title & initial labels ----------------------------------------
	axis.title  = spec.title
	axis.xlabel = _get_axis_label(spec.xlabel, spec.x_exp, spec.xscale[])
	axis.ylabel = _get_axis_label(spec.ylabel, spec.y_exp, spec.yscale[])

	# ---- Override global tick formatter for this specialized plot ----------
	axis.xtickformat[] = Makie.automatic
	axis.ytickformat[] = Makie.automatic

	# ---- Helpers ------------------------------------------------------------
	sanitize_log!(v::AbstractVector, is_log::Bool) =
		(is_log && !isempty(v)) ? (v[v .<= 0] .= NaN; v) : v

	_x_data_for(scale) = begin
		xd = _get_axis_data(spec.raw_freqs, spec.freqs, scale)
		sanitize_log!(xd.values, scale == Makie.log10)
		xd
	end

	_y_data_for(i::Int, scale) = begin
		yd = _get_axis_data(spec.raw_curves[i], spec.curves[i], scale)
		sanitize_log!(yd.values, scale == Makie.log10)
		yd
	end

	function _link_visibility!(plot_obj, controller)
		# plot_obj is the Errorbars plot object.
		# controller is the master Lines plot.
		# React to the controller's visibility changes.
		on(controller.visible) do is_visible
			# A. Manually control the visibility of the stem plot directly.
			plot_obj.visible = is_visible

			# B. Manually control the special attribute for the whiskers.
			plot_obj.whisker_visible[] = is_visible
		end
		nothing
	end

	# safe max(abs(.)) ignoring non-finite
	_finite_max_abs(v) = begin
		buf = (x -> abs(x)).(value.(v))
		any(isfinite, buf) ? maximum(x for x in buf if isfinite(x)) : 0.0
	end

	# ---- Select active (non-noise) curves by EPS -------------------------------
	ncurves    = length(spec.curves)
	active_idx = Int[]

	@inbounds for i in 1:ncurves
		# max magnitude of raw curve; works for Real, Complex, and Measurement types
		maxmag = maximum(value.(abs.(spec.raw_curves[i])))
		if maxmag > eps(Float64)          # keep only if anything rises above machine eps
			push!(active_idx, i)
		end
	end

	any_real_curve = !isempty(active_idx)

	# ---- Initial data (x) ---------------------------------------------------
	x_init = _x_data_for(spec.xscale[])
	x_vals_obs = Observable(copy(x_init.values))
	x_errs_obs = x_init.errors === nothing ? nothing : Observable(copy(x_init.errors))

	# ---- Per-curve allocs only for active curves ---------------------------
	palette = Makie.wong_colors()
	ncolors = length(palette)
	nact    = length(active_idx)

	y_vals_obs = Vector{Observable}(undef, nact)
	y_errs_obs = Vector{Union{Nothing, Observable}}(undef, nact)
	line_plots = Vector{Any}(undef, nact)
	yerr_plots = Vector{Any}(undef, nact)
	xerr_plots = Vector{Any}(undef, nact)

	# ---- Draw active curves -------------------------------------------------
	for k in 1:nact
		i = active_idx[k]
		color = palette[mod1(k, ncolors)]   # color by active order
		label = spec.labels[i]

		yd = _y_data_for(i, spec.yscale[])

		y_vals_obs[k] = Observable(copy(yd.values))
		y_errs_obs[k] = yd.errors === nothing ? nothing : Observable(copy(yd.errors))

		# line
		ln = lines!(
			axis,
			x_vals_obs,
			y_vals_obs[k];
			color = color,
			label = label,
			linewidth = 2,
		)
		line_plots[k] = ln

		# Y errorbars: stems + caps; fully follow the line’s visibility
		if y_errs_obs[k] !== nothing
			eb = errorbars!(
				axis, x_vals_obs, y_vals_obs[k], y_errs_obs[k];
				color = :black, direction = :y, whiskerwidth = 3, linewidth = 1,
			)
			_link_visibility!(eb, ln)
			yerr_plots[k] = eb
		else
			yerr_plots[k] = nothing
		end

		# X errorbars: stems + caps; fully follow the line’s visibility
		if x_errs_obs !== nothing
			ebx = errorbars!(
				axis, x_vals_obs, y_vals_obs[k], x_errs_obs;
				color = :black, direction = :x, whiskerwidth = 3, linewidth = 1,
			)
			_link_visibility!(ebx, ln)
			xerr_plots[k] = ebx
		else
			xerr_plots[k] = nothing
		end

	end

	# If nothing to draw, add transparent dummy without legend entry
	if !any_real_curve
		lines!(axis, x_vals_obs, [0]; color = :transparent, label = "No data")
	end

	# ---- Apply initial scales safely ---------------------------------------
	try
		axis.xscale[] = spec.xscale[]
		axis.yscale[] = spec.yscale[]
	catch
		axis.xscale[] = Makie.identity
		axis.yscale[] = Makie.identity
		@warn "Failed to set axis scale; reverted to linear scale."
	end

	# Enforce reasonable limits (avoid microscopic ranges when curves are flat)
	# Helper to compute finite extents
	_finite_extents(v::AbstractVector) = begin
		fv = filter(isfinite, v)
		isempty(fv) && return (NaN, NaN, false)
		return (minimum(fv), maximum(fv), true)
	end

	function _apply_limits!()
		# Helper: smallest positive finite value in a vector
		_min_positive(v::AbstractVector) = begin
			m = Inf
			@inbounds for a in v
				if isfinite(a) && a > 0 && a < m
					m = a
				end
			end
			return m
		end

		# X limits
		x = x_vals_obs[]
		xmin, xmax, okx = _finite_extents(x)
		if okx
			Δx = xmax - xmin
			if Δx <= 0
				xc = (xmax + xmin) / 2
				# minimal span based on magnitude
				Δx = max(1e-12, 1e-3 * max(abs(xc), abs(xmax), abs(xmin), 1.0))
				xmin = xc - Δx / 2
				xmax = xc + Δx / 2
			else
				pad = 0.05 * Δx
				xmin -= pad;
				xmax += pad
			end
			# Guard for log x-axis: lower bound must stay > 0
			if axis.xscale[] == Makie.log10
				posmin = _min_positive(x)
				floor_pos = isfinite(posmin) ? 0.9 * posmin : nextfloat(0.0)
				xmin = max(xmin, floor_pos)
				xmin <= 0 && (xmin = nextfloat(0.0))  # absolute safety
			end
			Makie.xlims!(axis, xmin, xmax)
		end

		# Y limits (consider error bars too)
		ymins = Float64[]
		ymaxs = Float64[]
		@inbounds for k in 1:nact
			y = y_vals_obs[k][]
			ymin, ymax, ok = _finite_extents(y)
			if ok
				if y_errs_obs[k] !== nothing
					e = y_errs_obs[k][]
					eymin, _, okm = _finite_extents(y .- e)
					_, eymax, okp = _finite_extents(y .+ e)
					okm && (ymin = min(ymin, eymin))
					okp && (ymax = max(ymax, eymax))
				end
				push!(ymins, ymin);
				push!(ymaxs, ymax)
			end
		end

		if !isempty(ymins)
			ymin = minimum(ymins)
			ymax = maximum(ymaxs)
			Δy = ymax - ymin
			yc = (ymax + ymin) / 2

			# Minimal span to avoid "micro-zoom" when the curve is essentially flat.
			#  - relative floor: 0.1% of magnitude (>= 1.0 to avoid collapsing near zero)
			#  - absolute floor: 1e-12
			min_span = max(1e-12, 1e-3 * max(abs(yc), abs(ymax), abs(ymin), 1.0))

			if !(Δy > min_span)
				Δy = min_span
				ymin = yc - Δy / 2
				ymax = yc + Δy / 2
			else
				pad = 0.05 * Δy
				ymin -= pad;
				ymax += pad
			end

			# Guard for log y-axis: lower bound must stay > 0
			if axis.yscale[] == Makie.log10
				# find smallest positive among all active curves (and their lower error bars)
				posmin = Inf
				@inbounds for k in 1:nact
					y = y_vals_obs[k][]
					m = _min_positive(y)
					if isfinite(m) && m < posmin
						posmin = m
					end
					if y_errs_obs[k] !== nothing
						e = y_errs_obs[k][]
						# consider lower whiskers
						@inbounds for (yy, ee) in zip(y, e)
							l = yy - ee
							if isfinite(l) && l > 0 && l < posmin
								posmin = l
							end
						end
					end
				end
				floor_pos = isfinite(posmin) ? 0.9 * posmin : nextfloat(0.0)
				ymin = max(ymin, floor_pos)
				ymin <= 0 && (ymin = nextfloat(0.0))  # absolute safety
			end

			Makie.ylims!(axis, ymin, ymax)
		end
		return nothing
	end
	Makie.autolimits!(axis)
	_apply_limits!()

	# ---- Refreshers (update Observables only) ------------------------------
	function _refresh_x!(scale)
		Makie.autolimits!(axis)
		spec.xscale[] = scale
		axis.xscale[] = scale
		axis.xlabel   = _get_axis_label(spec.xlabel, spec.x_exp, scale)

		xd = _x_data_for(scale)
		x_vals_obs[] = xd.values
		if x_errs_obs !== nothing
			x_errs_obs[] = xd.errors
		end

		_apply_limits!()
		nothing
	end

	function _refresh_y!(scale)
		Makie.autolimits!(axis)
		spec.yscale[] = scale
		axis.yscale[] = scale
		axis.ylabel   = _get_axis_label(spec.ylabel, spec.y_exp, scale)

		@inbounds for k in 1:nact
			i = active_idx[k]
			yd = _y_data_for(i, scale)
			y_vals_obs[k][] = yd.values
			if y_errs_obs[k] !== nothing
				y_errs_obs[k][] = yd.errors
			end
		end
		_apply_limits!()
		nothing
	end

	# ---- Buttons ------------------------------------------------------------
	buttons =
		any_real_curve ?
		[
			ControlButtonSpec(
				(_ctx, _btn) -> (Makie.reset_limits!(axis); nothing);
				icon = MI_REFRESH,
				on_success = ControlReaction(status_string = "Axis limits reset"),
			),
			ControlButtonSpec(
				(_ctx, _btn) -> _save_plot_export(spec, axis);
				icon = MI_SAVE,
				on_success = ControlReaction(
					status_string = path -> string("Saved SVG to ", basename(path)),
				),
			),
		] : Any[]

	# ---- Toggles ------------------------------------------------------------
	toggles =
		any_real_curve ?
		[
			ControlToggleSpec(
				(_ctx, _t) -> _refresh_x!(Makie.log10),
				(_ctx, _t) -> _refresh_x!(Makie.identity);
				label = "log x-axis",
				start_active = spec.xscale[] == Makie.log10,
				on_success_on = ControlReaction(status_string = "x-axis scale set to log"),
				on_success_off = ControlReaction(
					status_string = "x-axis scale set to linear",
				),
				on_failure = ControlReaction(status_string = err -> err),
			),
			ControlToggleSpec(
				(_ctx, _t) -> _refresh_y!(Makie.log10),
				(_ctx, _t) -> _refresh_y!(Makie.identity);
				label = "log y-axis",
				start_active = spec.yscale[] == Makie.log10,
				on_success_on = ControlReaction(status_string = "y-axis scale set to log"),
				on_success_off = ControlReaction(
					status_string = "y-axis scale set to linear",
				),
				on_failure = ControlReaction(status_string = err -> err),
			),
		] : Any[]

	# ---- Legend -------------------------------------------------------------
	legend_builder =
		parent ->
			Makie.Legend(
				parent,
				axis;
				orientation = :vertical,
			)

	return PlotBuildArtifacts(
		axis            = axis,
		legends         = legend_builder,
		colorbars       = Any[],
		control_buttons = buttons,
		control_toggles = toggles,
		status_message  = nothing,
	)
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
