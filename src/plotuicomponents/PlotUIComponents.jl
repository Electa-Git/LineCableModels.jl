module PlotUIComponents

using Makie
import ..BackendHandler: current_backend_symbol, _pkgid

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

const FIG_SIZE = (800, 600)
const FIG_PADDING = (80, 60, 40, 40) # left, right, bottom, top
const CTLBAR_HEIGHT = 36
const STATUSBAR_HEIGHT = 20
const GRID_ROW_GAP = 6
const GRID_COL_GAP = 6
const LEGEND_GAP = 4
const LEGEND_WIDTH = 140
const COLORBAR_GAP = 4
const CTLBAR_GAP = 2
const BUTTON_MIN_WIDTH = 32
const BUTTON_ICON_SIZE = 18
const BUTTON_TEXT_FONT_SIZE = 15
const AXIS_TITLE_FONT_SIZE = 15
const AXIS_LABEL_FONT_SIZE = 14
const AXIS_TICK_FONT_SIZE = 14
const STATUS_FONT_SIZE = 10
const BG_COLOR_INTERACTIVE = :grey90
const BG_COLOR_EXPORT = :white
const ICON_COLOR_ACTIVE = Makie.RGBAf(0.15, 0.15, 0.15, 1.0)
const ICON_COLOR_DISABLED = Makie.RGBAf(0.55, 0.55, 0.55, 1.0)


# -----------------------------------------------------------------------------
# Material UI icons
# -----------------------------------------------------------------------------
const MI_REFRESH = "\uE5D5"  # Material Icons: 'refresh'
const MI_SAVE    = "\uE161"  # Material Icons: 'save'
const ICON_TTF   = joinpath(@__DIR__, "..", "..", "assets", "fonts", "material-icons", "MaterialIcons-Regular.ttf")

# -----------------------------------------------------------------------------
# Data structures
# -----------------------------------------------------------------------------

mutable struct PlotBackendContext
	backend::Symbol
	interactive::Bool
	window::Union{Nothing, Any}
	screen::Union{Nothing, Any}
	use_latex_fonts::Bool
	icons::Function
	icons_font::Union{Nothing, String}
	statusbar::Union{Nothing, Makie.Observable{String}}
end

struct PlotFigureContext
	figure::Makie.Figure
	canvas_node::Any
	legend_grid::Makie.GridLayout
	legend_slot::Any
	colorbar_slot::Any
	ctlbar_node::Makie.GridLayout
	placeholder_node::Makie.GridLayout
	statusbar_node::Makie.GridLayout
end

struct ControlReaction
	status_string::Union{Nothing, String, Function}
	button_color::Union{Nothing, Any}
	button_label::Union{Nothing, String}
	timeout::Union{Nothing, AbstractFloat}
end

ControlReaction(;
	status_string = nothing,
	button_color = nothing,
	button_label = nothing,
	timeout = 1.5,
) =
	ControlReaction(status_string, button_color, button_label, timeout)

struct ControlButtonSpec
	label::Union{Nothing, String}
	icon::Union{Nothing, String}
	action::Function
	on_success::Union{Nothing, ControlReaction}
	on_failure::Union{Nothing, ControlReaction}
end

ControlButtonSpec(
	action::Function;
	label::Union{Nothing, String} = nothing,
	icon::Union{Nothing, String} = nothing,
	on_success::Union{Nothing, ControlReaction} = nothing,
	on_failure::Union{Nothing, ControlReaction} = nothing,
) = ControlButtonSpec(label, icon, action, on_success, on_failure)

struct ControlToggleSpec
	label::Union{Nothing, String}
	action_on::Function
	action_off::Function
	on_success_on::Union{Nothing, ControlReaction}
	on_success_off::Union{Nothing, ControlReaction}
	on_failure::Union{Nothing, ControlReaction}
	start_active::Bool
end

ControlToggleSpec(
	action_on::Function,
	action_off::Function;
	label::Union{Nothing, String} = nothing,
	on_success_on::Union{Nothing, ControlReaction} = nothing,
	on_success_off::Union{Nothing, ControlReaction} = nothing,
	on_failure::Union{Nothing, ControlReaction} = nothing,
	start_active::Bool = false,
) = ControlToggleSpec(
	label,
	action_on,
	action_off,
	on_success_on,
	on_success_off,
	on_failure,
	start_active,
)

struct PlotBuildArtifacts
	axis::Union{Nothing, Makie.Axis}
	legends::Union{Nothing, Any}
	colorbars::Union{Nothing, Vector{Any}}
	control_buttons::Vector{ControlButtonSpec}
	control_toggles::Vector{ControlToggleSpec}
	status_message::Union{Nothing, String}
end

PlotBuildArtifacts(; axis = nothing, legends = nothing, colorbars = nothing,
	control_buttons = ControlButtonSpec[], control_toggles = ControlToggleSpec[],
	status_message = nothing) =
	PlotBuildArtifacts(
		axis,
		legends,
		colorbars,
		control_buttons,
		control_toggles,
		status_message,
	)

struct PlotAssembly
	backend_ctx::PlotBackendContext
	figure_ctx::PlotFigureContext
	figure::Makie.Figure
	axis::Any
	buttons::Vector{Makie.Button}
	legend::Any
	colorbars::Vector{Any}
	status_label::Any
	artifacts::PlotBuildArtifacts
end

# -----------------------------------------------------------------------------
# Backend helpers
# -----------------------------------------------------------------------------

"""Create a GLMakie screen if GL backend is active; otherwise return nothing."""
function gl_screen(title::AbstractString)
	if current_backend_symbol() == :gl
		mod = Base.require(_pkgid(:gl))
		ctor = getproperty(mod, :Screen)
		return Base.invokelatest(ctor; title = String(title))
	end
	return nothing
end

# tiny helper to build "icon + text" labels ergonomically ---
"""
with_icon(icon; text="", isize=14, tsize=12, color=:black, gap=4,
		  dy_icon=-0.18, dy_text=0.0)

- `dy_icon`, `dy_text`: vertical tweaks in *em* units (fraction of that part's fontsize).
  Negative moves down, positive moves up.
"""
with_icon(icon::AbstractString; text::AbstractString = "",
	isize::Int = BUTTON_ICON_SIZE, tsize::Int = BUTTON_TEXT_FONT_SIZE, color = :black,
	gap::Int = 2,
	dy_icon::Float64 = -0.18, dy_text::Float64 = 0.0) =
	text == "" ?
	rich(icon; font = :icons, fontsize = isize, color = color, offset = (0, dy_icon)) :
	rich(
		rich(icon; font = :icons, fontsize = isize, color = color, offset = (0, dy_icon)),
		rich(" "^gap; font = :regular, fontsize = tsize, color = color),
		rich(text; font = :regular, fontsize = tsize, color = color, offset = (0, dy_text)),
	)

function build_backend_context(
	backend::Symbol;
	interactive::Union{Nothing, Bool} = nothing,
	window = nothing,
	screen = nothing,
	icons::Function = (icon; text = nothing, kwargs...) ->
		(text === nothing ? string(icon) : string(text)),
	use_latex_fonts::Bool = false,
	icons_font::Union{Nothing, String} = nothing,
	statusbar::Union{Nothing, Makie.Observable{String}} = nothing,
)
	is_interactive =
		interactive === nothing ? backend in (:gl, :wgl, :wglmakie) : interactive
	chan = statusbar
	if chan === nothing && is_interactive
		chan = Makie.Observable("")
	end
	return PlotBackendContext(
		backend,
		is_interactive,
		window,
		screen,
		use_latex_fonts,
		icons,
		icons_font,
		chan,
	)
end

function attach_window!(ctx::PlotBackendContext; window = nothing, screen = nothing)
	ctx.window = window
	ctx.screen = screen
	return ctx
end

function _make_window(
	backend_handler::Module,
	backend::Union{Nothing, Symbol} = nothing;
	title::AbstractString = "LineCableModels Plot",
	icons::Function = (icon; text = nothing, kwargs...) ->
		(text === nothing ? string(icon) : string(text)),
	use_latex_fonts::Bool = false,
	icons_font::Union{Nothing, String} = nothing,
	statusbar::Union{Nothing, Makie.Observable{String}} = nothing,
	interactive_override::Union{Nothing, Bool} = nothing,
)
	actual_backend = backend_handler.ensure_backend!(backend)
	is_interactive =
		interactive_override === nothing ? actual_backend in (:gl, :wgl) :
		interactive_override
	ctx = build_backend_context(
		actual_backend;
		interactive = is_interactive,
		icons = icons,
		use_latex_fonts = use_latex_fonts,
		icons_font = icons_font,
		statusbar = statusbar,
	)
	if is_interactive && actual_backend == :gl
		scr = gl_screen(title)
		if scr !== nothing
			attach_window!(ctx; window = scr, screen = scr)
		end
	end
	return ctx
end

function theme_for(
	ctx::PlotBackendContext;
	mode::Symbol = ctx.interactive ? :interactive : :export,
)
	background = mode === :interactive ? BG_COLOR_INTERACTIVE : BG_COLOR_EXPORT
	base = Makie.Theme()
	if ctx.use_latex_fonts && mode == :export
		base = merge(base, Makie.theme_latexfonts())
	end
	icon_font = ctx.icons_font
	custom =
		icon_font === nothing ?
		Makie.Theme(
			backgroundcolor = background,
			Axis = (
				titlesize = AXIS_TITLE_FONT_SIZE,
				xlabelsize = AXIS_LABEL_FONT_SIZE,
				ylabelsize = AXIS_LABEL_FONT_SIZE,
				xticklabelsize = AXIS_TICK_FONT_SIZE,
				yticklabelsize = AXIS_TICK_FONT_SIZE,
			),
			Legend = (
				fontsize = AXIS_LABEL_FONT_SIZE,
				labelsize = AXIS_LABEL_FONT_SIZE,
			),
			Colorbar = (
				labelsize = AXIS_LABEL_FONT_SIZE,
				ticklabelsize = AXIS_TICK_FONT_SIZE,
			),
		) :
		Makie.Theme(
			backgroundcolor = background,
			fonts = (; icons = icon_font),
			Axis = (
				titlesize = AXIS_TITLE_FONT_SIZE,
				xlabelsize = AXIS_LABEL_FONT_SIZE,
				ylabelsize = AXIS_LABEL_FONT_SIZE,
				xticklabelsize = AXIS_TICK_FONT_SIZE,
				yticklabelsize = AXIS_TICK_FONT_SIZE,
			),
			Legend = (
				fontsize = AXIS_LABEL_FONT_SIZE,
				labelsize = AXIS_LABEL_FONT_SIZE,
			),
			Colorbar = (
				labelsize = AXIS_LABEL_FONT_SIZE,
				ticklabelsize = AXIS_TICK_FONT_SIZE,
			),
		)
	return merge(base, custom)
end

_configure_theme!(
	ctx::PlotBackendContext;
	mode::Symbol = ctx.interactive ? :interactive : :export,
) =
	theme_for(ctx; mode = mode)

function with_plot_theme(
	f::Function,
	ctx::PlotBackendContext;
	mode::Union{Nothing, Symbol} = nothing,
)
	chosen_mode = mode === nothing ? (ctx.interactive ? :interactive : :export) : mode
	theme = _configure_theme!(ctx; mode = chosen_mode)
	return Makie.with_theme(theme) do
		f()
	end
end

# -----------------------------------------------------------------------------
# Figure helpers
# -----------------------------------------------------------------------------

function _make_figure(
	ctx::PlotBackendContext;
	fig_size::Tuple{Int, Int} = FIG_SIZE,
	figure_padding::NTuple{4, Int} = FIG_PADDING,
	legend_panel_width::Int = LEGEND_WIDTH,
)
	fig = Makie.Figure(; size = fig_size, figure_padding = figure_padding)

	ctlbar_node = fig[1, 1:2] = Makie.GridLayout()
	ctlbar_node.halign = :left
	ctlbar_node.valign = :bottom
	placeholder_node = fig[2, 1:2] = Makie.GridLayout()
	canvas_node = fig[3, 1]
	legend_grid = fig[3, 2] = Makie.GridLayout()
	statusbar_node = fig[4, 1:2] = Makie.GridLayout()
	statusbar_node.halign = :left

	legend_slot = legend_grid[1, 1]
	legend_slot[] = Makie.GridLayout()
	colorbar_slot = legend_grid[2, 1]
	colorbar_slot[] = Makie.GridLayout()

	fig_ctx = PlotFigureContext(
		fig,
		canvas_node,
		legend_grid,
		legend_slot,
		colorbar_slot,
		ctlbar_node,
		placeholder_node,
		statusbar_node,
	)

	_configure_layout!(
		fig_ctx;
		interactive = ctx.interactive,
		legend_panel_width = legend_panel_width,
	)
	return fig_ctx
end

function _configure_layout!(
	fig_ctx::PlotFigureContext;
	interactive::Bool = true,
	legend_panel_width::Int = LEGEND_WIDTH,
)
	layout = fig_ctx.figure.layout

	Makie.rowgap!(layout, GRID_ROW_GAP)
	Makie.colgap!(layout, GRID_COL_GAP)

	Makie.rowsize!(layout, 1, Makie.Fixed(interactive ? CTLBAR_HEIGHT : 0))
	Makie.rowsize!(layout, 2, Makie.Fixed(0))
	Makie.rowsize!(layout, 3, Makie.Relative(1.0))
	Makie.rowsize!(layout, 4, Makie.Fixed(interactive ? STATUSBAR_HEIGHT : 0))

	Makie.colsize!(layout, 1, Makie.Relative(1.0))
	Makie.colsize!(layout, 2, Makie.Fixed(legend_panel_width))

	Makie.rowgap!(fig_ctx.legend_grid, LEGEND_GAP)
	Makie.colgap!(fig_ctx.legend_grid, 0)


	Makie.rowsize!(fig_ctx.legend_grid, 1, Makie.Auto())
	Makie.rowsize!(fig_ctx.legend_grid, 2, Makie.Auto())

	return fig_ctx
end

function _make_canvas!(
	fig_ctx::PlotFigureContext;
	axis_ctor = Makie.Axis,
	axis_options::NamedTuple = NamedTuple(),
)
	axis = axis_ctor(fig_ctx.canvas_node; axis_options...)
	return axis
end

# -----------------------------------------------------------------------------
# Control bar helpers
# -----------------------------------------------------------------------------

function _make_ctlbar!(
	fig_ctx::PlotFigureContext,
	ctx::PlotBackendContext,
	button_specs::AbstractVector{ControlButtonSpec},
	toggle_specs::AbstractVector{ControlToggleSpec};
	button_height::Int = max(CTLBAR_HEIGHT - 12, 32),
	button_gap::Int = CTLBAR_GAP,
)
	if !ctx.interactive || (isempty(button_specs) && isempty(toggle_specs))
		Makie.rowsize!(fig_ctx.figure.layout, 1, Makie.Fixed(0))
		return [], []
	end

	layout = fig_ctx.ctlbar_node
	Makie.rowgap!(layout, 0)
	Makie.colgap!(layout, button_gap)
	Makie.rowsize!(layout, 1, Makie.Fixed(button_height))

	buttons = Makie.Button[]
	toggles = Makie.Toggle[]
	col_idx = 1

	for spec in button_specs
		label = _build_button_label(ctx, spec)
		button_kwargs = (
			; label = label,
			fontsize = BUTTON_TEXT_FONT_SIZE,
			height = button_height,
			halign = :left,
		)
		if spec.icon !== nothing
			width = _preferred_button_width(spec)
			if width !== nothing
				button_kwargs = (; button_kwargs..., width = width)
			end
		end
		button = Makie.Button(layout[1, col_idx]; button_kwargs...)
		push!(buttons, button)
		_wire_button_callback!(button, spec, ctx)
		col_idx += 1
	end

	for spec in toggle_specs
		gl = layout[1, col_idx] = Makie.GridLayout()
		gl.halign = :left
		gl.valign = :center

		toggle = Makie.Toggle(gl[1, 2]; active = spec.start_active)
		if spec.label !== nothing
			Makie.Label(gl[1, 1], spec.label, halign = :right)
		end

		push!(toggles, toggle)
		_wire_toggle_callback!(toggle, spec, ctx)
		col_idx += 1
	end

	return buttons, toggles
end

function _build_button_label(ctx::PlotBackendContext, spec::ControlButtonSpec)
	icon_fn = ctx.icons
	label_text = spec.label === nothing ? "" : spec.label
	if spec.icon === nothing
		return label_text
	end
	try
		return icon_fn(spec.icon; text = label_text, gap = 6)
	catch err
		if err isa MethodError
			return isempty(label_text) ? string(spec.icon) : label_text
		else
			rethrow(err)
		end
	end
end

function _preferred_button_width(spec::ControlButtonSpec)
	if spec.icon !== nothing && spec.label === nothing
		return BUTTON_MIN_WIDTH
	end
	return nothing
end

function _wire_button_callback!(button, spec::ControlButtonSpec, ctx::PlotBackendContext)
	ensure_statusbar!(ctx)

	Makie.on(button.clicks) do _
		Base.@async begin
			try
				result = _invoke_button_action(spec.action, ctx, button)
				_apply_reaction!(ctx, button, spec.on_success, result)
			catch err
				_apply_reaction!(ctx, button, spec.on_failure, sprint(showerror, err))
			end
		end
	end
	return button
end

function _wire_toggle_callback!(toggle, spec::ControlToggleSpec, ctx::PlotBackendContext)
	ensure_statusbar!(ctx)

	Makie.on(toggle.active) do is_active
		Base.@async begin
			original_state = !is_active
			try
				if is_active
					result = _invoke_button_action(spec.action_on, ctx, toggle)
					_apply_reaction!(ctx, toggle, spec.on_success_on, result)
				else
					result = _invoke_button_action(spec.action_off, ctx, toggle)
					_apply_reaction!(ctx, toggle, spec.on_success_off, result)
				end
			catch err
				_apply_reaction!(ctx, toggle, spec.on_failure, sprint(showerror, err))
			end
		end
	end
	return toggle
end

function _invoke_button_action(action::Function, ctx::PlotBackendContext, button)
	try
		return Base.invokelatest(action, ctx, button)
	catch err
		if err isa MethodError && err.f === action
			try
				return Base.invokelatest(action, ctx)
			catch err2
				if err2 isa MethodError && err2.f === action
					return Base.invokelatest(action)
				else
					throw(err2)
				end
			end
		else
			throw(err)
		end
	end
end


function _apply_reaction!(
	ctx::PlotBackendContext,
	button,
	reaction::Union{Nothing, ControlReaction},
	result,
)
	has_color = hasproperty(button, :buttoncolor)
	original_color = has_color ? button.buttoncolor[] : nothing
	has_label = hasproperty(button, :label)
	original_label = has_label ? button.label[] : nothing

	# Determine the status message
	status_msg = nothing
	if reaction !== nothing && reaction.status_string !== nothing
		if reaction.status_string isa Function
			status_msg = reaction.status_string(result)
		else
			status_msg = reaction.status_string
		end
	elseif result isa AbstractString && !isempty(result)
		status_msg = result
	end

	# Apply reaction and status update
	if status_msg !== nothing
		update_status!(ctx, status_msg)
	end

	if reaction !== nothing
		if reaction.button_color !== nothing && has_color
			button.buttoncolor[] = Makie.to_color(reaction.button_color)
		end
		if reaction.button_label !== nothing
			button.label[] = reaction.button_label
		end
	end

	# Handle timeout and UI restoration
	timeout = reaction !== nothing ? reaction.timeout : 1.6
	if timeout !== nothing && isfinite(timeout)
		sleep(timeout)
		if status_msg !== nothing
			clear_status!(ctx)
		end
		if reaction !== nothing
			if reaction.button_color !== nothing && has_color
				button.buttoncolor[] = original_color
			end
			if reaction.button_label !== nothing && has_label
				button.label[] = original_label
			end
		end
	end

	return nothing
end

clear_status!(ctx) = begin
	# non-breaking space keeps the row height while looking empty
	update_status!(ctx, "\u00A0")

end

# -----------------------------------------------------------------------------
# Legend & colorbar helpers
# -----------------------------------------------------------------------------

"""Populate the legend area. Accepts `nothing`, a builder function, or a Makie plot object."""
function _make_legend!(fig_ctx::PlotFigureContext, content; kwargs...)
	slot = fig_ctx.legend_slot
	if content === nothing
		slot[] = Makie.GridLayout()
		Makie.rowsize!(fig_ctx.legend_grid, 1, Makie.Fixed(0))
		return nothing
	end

	Makie.rowsize!(fig_ctx.legend_grid, 1, Makie.Auto())
	container = Makie.GridLayout()
	slot[] = container
	built = _materialize_component!(container[1, 1], content; kwargs...)
	if built !== nothing
		if hasproperty(built, :valign)
			built.valign[] = :top
		end
		if hasproperty(built, :halign)
			built.halign[] = :left
		end
	end
	return built
end

"""Populate the colorbar stack with zero or more builder specs."""
function _make_colorbars!(
	fig_ctx::PlotFigureContext,
	specs::Union{Nothing, AbstractVector};
	kwargs...,
)
	slot = fig_ctx.colorbar_slot
	if specs === nothing || isempty(specs)
		slot[] = Makie.GridLayout()
		Makie.rowsize!(fig_ctx.legend_grid, 2, Makie.Fixed(0))
		return Any[]
	end

	Makie.rowsize!(fig_ctx.legend_grid, 2, Makie.Auto())
	container = Makie.GridLayout()
	slot[] = container
	Makie.rowgap!(container, COLORBAR_GAP)

	built = Any[]
	row = 1
	for spec in specs
		spec === nothing && continue
		node = container[row, 1]
		push!(built, _materialize_component!(node, spec; kwargs...))
		row += 1
	end
	return built
end

function _materialize_component!(parent, spec; kwargs...)
	if spec isa Function
		return spec(parent; kwargs...)
	elseif Makie.isplot(spec)
		parent[] = spec
		return spec
	else
		try
			parent[] = spec
			return spec
		catch err
			if err isa MethodError
				error("Unsupported component specification $(typeof(spec))")
			else
				rethrow(err)
			end
		end
	end
end


# -----------------------------------------------------------------------------
# Status helpers
# -----------------------------------------------------------------------------

function _make_statusbar!(
	fig_ctx::PlotFigureContext,
	ctx::PlotBackendContext;
	initial_message::AbstractString = "",
)
	if !ctx.interactive
		Makie.rowsize!(fig_ctx.figure.layout, 4, Makie.Fixed(0))
		return nothing
	end

	status_obs = ensure_statusbar!(ctx)
	if !isempty(initial_message)
		status_obs[] = String(initial_message)
	end

	label = Makie.Label(fig_ctx.statusbar_node[1, 1];
		text = status_obs,
		fontsize = STATUS_FONT_SIZE,
		halign = :left,
		tellwidth = false,
		tellheight = false,
	)
	return label
end

function ensure_statusbar!(ctx::PlotBackendContext)
	if ctx.statusbar === nothing
		ctx.statusbar = Makie.Observable("")
	end
	return ctx.statusbar
end

function update_status!(ctx::PlotBackendContext, message::AbstractString)
	chan = ensure_statusbar!(ctx)
	chan[] = String(message)
	return chan
end

# -----------------------------------------------------------------------------
# Orchestration helpers
# -----------------------------------------------------------------------------

function _run_plot_pipeline(
	backend_ctx::PlotBackendContext,
	plot_fn::Function;
	fig_size::Tuple{Int, Int} = FIG_SIZE,
	figure_padding::NTuple{4, Int} = FIG_PADDING,
	legend_panel_width::Int = LEGEND_WIDTH,
	axis_ctor = Makie.Axis,
	axis_kwargs::NamedTuple = NamedTuple(),
	extra_buttons::AbstractVector{ControlButtonSpec} = ControlButtonSpec[],
	initial_status::Union{Nothing, String} = nothing,
)
	fig_ctx = _make_figure(
		backend_ctx;
		fig_size = fig_size,
		figure_padding = figure_padding,
		legend_panel_width = legend_panel_width,
	)

	axis =
		isempty(axis_kwargs) ?
		_make_canvas!(fig_ctx; axis_ctor = axis_ctor) :
		_make_canvas!(fig_ctx; axis_ctor = axis_ctor, axis_options = axis_kwargs)

	artifacts = plot_fn(fig_ctx, backend_ctx, axis)
	artifacts = artifacts === nothing ? PlotBuildArtifacts(axis = axis) : artifacts

	axis = artifacts.axis === nothing ? axis : artifacts.axis

	button_specs = ControlButtonSpec[]
	isempty(extra_buttons) || append!(button_specs, extra_buttons)
	isempty(artifacts.control_buttons) || append!(button_specs, artifacts.control_buttons)

	buttons, toggles =
		_make_ctlbar!(fig_ctx, backend_ctx, button_specs, artifacts.control_toggles)

	legend_obj = _make_legend!(fig_ctx, artifacts.legends)
	colorbar_objs = _make_colorbars!(fig_ctx, artifacts.colorbars)

	status_message = artifacts.status_message
	if status_message === nothing
		status_message = initial_status
	end
	status_message = status_message === nothing ? "" : status_message

	status_label = _make_statusbar!(fig_ctx, backend_ctx; initial_message = status_message)
	if !isempty(status_message)
		update_status!(backend_ctx, status_message)
	end

	return PlotAssembly(
		backend_ctx,
		fig_ctx,
		fig_ctx.figure,
		axis,
		buttons,
		legend_obj,
		colorbar_objs,
		status_label,
		artifacts,
	)
end

make_window_context(args...; kwargs...) = _make_window(args...; kwargs...)
make_standard_figure(args...; kwargs...) = _make_figure(args...; kwargs...)
configure_layout!(args...; kwargs...) = _configure_layout!(args...; kwargs...)
make_canvas!(args...; kwargs...) = _make_canvas!(args...; kwargs...)
make_ctlbar!(args...; kwargs...) = _make_ctlbar!(args...; kwargs...)
make_legend!(args...; kwargs...) = _make_legend!(args...; kwargs...)
make_colorbars!(args...; kwargs...) = _make_colorbars!(args...; kwargs...)
make_statusbar!(args...; kwargs...) = _make_statusbar!(args...; kwargs...)
run_plot_pipeline(args...; kwargs...) = _run_plot_pipeline(args...; kwargs...)

function ensure_export_background!(fig)
	if fig !== nothing && hasproperty(fig, :scene)
		fig.scene.backgroundcolor[] = Makie.to_color(BG_COLOR_EXPORT)
	end
	return fig
end

end # module PlotUIComponents
