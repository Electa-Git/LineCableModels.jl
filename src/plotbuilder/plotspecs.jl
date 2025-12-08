function _sanitize_filename_plot(str::AbstractString)
	sanitized = lowercase(strip(str))
	sanitized = replace(sanitized, r"[^0-9a-z]+" => "_")
	sanitized = strip(sanitized, '_')
	return isempty(sanitized) ? "linecablemodels_plot" : sanitized
end

function _default_export_path(
	spec::AbstractPlotSpec;
	extension::AbstractString = EXPORT_EXTENSION,
)
	base_title = strip(spec.title)
	base = isempty(base_title) ? string(spec.parent_kind, "_", spec.component) : base_title
	name = _sanitize_filename_plot(base)
	timestamp = format(now(), EXPORT_TIMESTAMP_FORMAT)
	filename = string(name, "_", timestamp, ".", extension)
	return joinpath(pwd(), filename)
end

function _save_plot_export(spec::AbstractPlotSpec, axis)
	# Capture current axis scales before building the export figure
	spec.xscale[] = axis.xscale[]
	spec.yscale[] = axis.yscale[]
	fig = build_export_figure(spec)
	trim!(fig.layout)
	path = _default_export_path(spec)
	Makie.save(path, fig)
	return path
end

"""
Parses a values expression like :X[1,1,:] or :X[1,1,1:5].
"""
function _parse_values_expr(values_expr, ijk)

	# --- input is :R[1,1,:] ---
	if ijk === nothing
		if values_expr isa Expr && values_expr.head === :ref &&
		   length(values_expr.args) == 4
			q = values_expr.args[1]
			q isa Symbol ||
				Base.error("Expected values symbol as first argument, got $(q)")

			local i::Int
			local j::Int
			local k::Union{Int, Colon, AbstractRange}

			# Eval the indices to resolve them from Expr.
			# It will correctly resolve 1, :, and 1:5.
			try
				i = eval(values_expr.args[2])
				j = eval(values_expr.args[3])
				k = eval(values_expr.args[4])
			catch e
				Base.error(
					"Failed to parse indices from $(values_expr). Ensure they are valid literals (1, :, 1:5, etc.). Error: $e",
				)
			end

			# Type checking after eval
			i isa Int || Base.error("i index '$i' is not an Int")
			j isa Int || Base.error("j index '$j' is not an Int")
			(k isa Int || k == (:) || k isa AbstractRange) ||
				Base.error("k index '$k' must be Int, ':', or AbstractRange")

			return q, (i, j, k)
		else
			Base.error(
				"Provide values as Expr like :X[1,1,:] or :X[1,1,1:5], or pass symbol and `ijk`",
			)
		end

		# --- input is :X, ijk=(1,1,:) ---
	else
		# This branch already supports non-Int types, it just needs a type assertion.
		ijk isa NTuple{3, Any} ||
			Base.error("ijk must be NTuple{3,Any}, got $(typeof(ijk))")
		values_expr isa Symbol ||
			Base.error(
				"values must be Symbol when ijk is provided; got $(typeof(values_expr))",
			)

		i, j, k = ijk

		# Check types
		i isa Int || Base.error("i in ijk must be Int")
		j isa Int || Base.error("j in ijk must be Int")
		(k isa Int || k == (:) || k isa AbstractRange) ||
			Base.error("k in ijk must be Int, ':', or AbstractRange")

		return values_expr, (i, j, k)
	end
end

function _autoscale_axis(values::AbstractVector{<:Real}; _threshold = 1e4)
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
	_threshold_exp = floor(Int, log10(_threshold))
	abs(exp) < abs(_threshold_exp) && return values, 0
	scale = 10.0 ^ exp
	return values ./ scale, exp
end

function _render_plot_specs(
	spec::AbstractPlotSpec;
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
			# This closure dispatches to the correct _build_plot! method
			# based on the *concrete type* of 'spec'.
			(fig_ctx, ctx, axis) -> _build_plot!(fig_ctx, ctx, axis, spec);
			pipeline_kwargs...,
		)
	end
	if display_plot
		_display!(backend_ctx, assembly.figure; title = spec.title)
	end
	return assembly
end

# --- De-duplicated Button Logic ---
function _build_common_plot_controls(spec::AbstractPlotSpec, axis)
	buttons = [
		ControlButtonSpec(
			(_ctx, _btn) -> (Makie.autolimits!(axis); nothing),
			icon = MI_REFRESH,
			on_success = ControlReaction(status_string = "Axis limits reset"),
		),
		ControlButtonSpec(
			(_ctx, _btn) -> _save_plot_export(spec, axis), # Generic save
			icon = MI_SAVE,
			on_success = ControlReaction(
				status_string = path -> string("Saved SVG to ", basename(path)),
			),
		),
	]
	return buttons
end
