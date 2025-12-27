# -------------------------
# Context + display routing
# -------------------------

function build_context(;
	backend::Union{Nothing, Symbol} = nothing,
	display::Bool = true,
	interactive_override::Union{Nothing, Bool} = nothing,
	title::AbstractString = "LineCableModels Plot",
	status::Union{Nothing, Makie.Observable{String}} = nothing,
	window = nothing,
	screen = nothing,
	theme::Union{Nothing, Makie.Theme} = nothing,
	kwargs...,
)
	actual_backend = BackendHandler.ensure_backend!(backend)

	interactive =
		interactive_override === nothing ?
		(display && actual_backend in (:gl, :wgl)) :
		interactive_override

	chan = status
	if chan === nothing && interactive
		chan = Makie.Observable("")
	end

	win = window
	scr = screen

	if interactive && actual_backend == :gl && win === nothing && scr === nothing
		scr = BackendHandler.make_screen(title; backend = :gl)
		if scr !== nothing
			win = scr
		end
	end

	return UIContext(actual_backend, interactive, win, scr, chan, theme)
end

function display!(ctx::UIContext, fig::Makie.Figure; title::AbstractString = "")
	if ctx.interactive && ctx.window !== nothing
		display(ctx.window, fig)
		if !isempty(title) && hasproperty(ctx.window, :title)
			try
				ctx.window.title[] = String(title)
			catch
			end
		end
	else
		BackendHandler.renderfig(fig)
	end
	return nothing
end

function display!(ctx::UIContext, assem::PlotAssembly)
	ttl = assem.pbfig.title
	return display!(ctx, assem.uifig.figure; title = ttl)
end

# -------------------------
# Figure / panel build
# -------------------------

function render(r::RenderSpec; backend = nothing, display::Bool = true, kwargs...)
	ctx = build_context(; backend = backend, display = display, kwargs...)
	assemblies = PlotAssembly[]
	for pbfig in r.figures
		assem = render_figure(ctx, r.spec, pbfig; display = display, kwargs...)
		push!(assemblies, assem)
	end
	return assemblies
end

function render_figure(
	ctx::UIContext,
	::Type{S},
	pbfig::PageSpec;
	display::Bool = true,
	kwargs...,
) where {S}
	ls = layout_spec(S, Val(pbfig.layout))
	uifig = build_figure(ls, ctx, pbfig)
	panels = UIPanel[]
	for v in pbfig.views
		push!(panels, build_panel!(uifig, ctx, pbfig, v))
	end
	widgets = build_widgets!(uifig, ctx, S, pbfig, panels; kwargs...)
	assem = PlotAssembly(S, ctx, pbfig, uifig, panels, widgets, (;))
	display && display!(ctx, assem)
	return assem
end

function build_figure(ls::UILayoutSpec, ctx::UIContext, pbfig::PageSpec)
	# Clean up kwargs
	kw = pbfig.kwargs
	kw2 = (; (k => v for (k, v) in pairs(kw) if k != :size && k != :resolution)...)

	fig = Makie.Figure(; size = pbfig.size, kw2...)

	containers = Dict{Symbol, Makie.GridLayout}()
	slots = Dict{Symbol, Any}()

	# Root is always initialized
	root = fig.layout
	containers[:root] = root

	# --- PHASE 1: Materialize Slots as Nested Grids ---
	for s in ls.slots
		parent_gl = containers[s.parent]

		# 1. Get the position
		gp = parent_gl[s.at...]

		# 2. IMMEDIATELY initialize a sub-grid in this slot.
		# This forces the parent_gl to expand to include this row/col.
		subgl = Makie.GridLayout()
		gp[] = subgl

		# 3. Store the sub-grid, not the position.
		slots[s.name] = subgl

		# 4. Now safe to apply sizes to the PARENT because it has expanded.
		if haskey(s.props, :rowsize) && s.at[1] isa Int
			Makie.rowsize!(parent_gl, s.at[1], s.props.rowsize)
		end
		if haskey(s.props, :colsize) && s.at[2] isa Int
			Makie.colsize!(parent_gl, s.at[2], s.props.colsize)
		end
	end

	# Determine panelgrid shape for multi-view scenarios
	n = length(pbfig.views)
	panel_shape = if n > 1
		nr = ceil(Int, sqrt(n))
		nc = ceil(Int, n / nr)
		(nr, nc)
	else
		(1, 1)
	end

	return UIFigure(fig, ls, containers, slots, Ref(0), panel_shape)
end

function build_panel!(uifig::UIFigure, ctx::UIContext, pbfig::PageSpec, view::ViewSpec)
	# Determine the semantic slot
	slot_key =
		(length(pbfig.views) > 1 || !haskey(uifig.slots, :canvas)) ? :panelgrid : :canvas

	# Retrieve the target GridLayout (guaranteed to be a GL now)
	target_gl = uifig.slots[slot_key]

	# Calculate placement INSIDE the slot's grid
	# If it's a panelgrid, we calculate (i, j). If single view, it's just (1, 1).
	ax_pos = if slot_key == :panelgrid && length(pbfig.views) > 1
		uifig.cursor[] += 1
		k = uifig.cursor[]
		nr, nc = uifig.panelshape
		row = (k - 1) ÷ nc + 1
		col = (k - 1) % nc + 1
		target_gl[row, col]
	else
		target_gl[1, 1]
	end

	# Axis creation
	ax = Makie.Axis(ax_pos;
		xlabel = something(view.xaxis.label, ""),
		ylabel = something(view.yaxis.label, ""),
		title  = view.title,
		xscale = (view.xaxis.scale == :log10) ? Makie.log10 : Makie.identity,
		yscale = (view.yaxis.scale == :log10) ? Makie.log10 : Makie.identity,
	)

	plots = Any[]
	for s in view.series
		append!(plots, draw!(ax, s))
	end

	if any(s -> s.label !== nothing, view.series)
		Makie.axislegend(ax)
	end

	return UIPanel(view, ax, plots, (;))
end

function export_svg!(ctx::UIContext, assem::PlotAssembly, path::AbstractString; kwargs...)
	return action_export_svg!(ctx, assem, path; kwargs...)
end
