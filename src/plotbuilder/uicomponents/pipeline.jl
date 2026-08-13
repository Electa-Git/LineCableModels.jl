# -------------------------
# The Architect: make_context
# -------------------------

function make_context(;
	backend::Union{Nothing, Symbol} = nothing,
	display::Bool = true,
	title::AbstractString = "LineCableModels Plot",
	theme::Union{Nothing, Makie.Theme} = nothing,
	use_latex_fonts::Bool = false,
	kwargs...,
)
	active_backend = BackendHandler.ensure_backend!(backend)
	interactive = (display && active_backend in (:gl, :wgl))

	stat = interactive ? Makie.Observable("Ready.") : nothing

	win = nothing
	scr = nothing
	if interactive && active_backend == :gl
		scr = BackendHandler.make_screen(title; backend = :gl)
		win = scr
	end

	# Default theme
	default_theme = make_theme(; interactive, use_latex_fonts)

	# 2. Build official theme (uses ctx.interactive by default)
	built_theme = theme === nothing ? default_theme : merge(default_theme, theme)

	return UIContext(
		active_backend,
		interactive,
		use_latex_fonts,
		win,
		scr,
		stat,
		built_theme,
	)
end

# -------------------------
# The Boss: render
# -------------------------

function build(
	r::RenderSpec{S};
	backend = nothing,
	display::Bool = true,
	kwargs...,
) where {S}
	ctx = make_context(; backend = backend, display = display, kwargs...)
	assemblies = UIPlot[]

	Makie.with_theme(ctx.theme) do
		for page in r.figures # page is PageSpec

			# A. Architect
			layout_s = make_layout(Val(page.layout))

			# B. Constructor (Shell)
			uifig = build_figure(ctx, page, layout_s)

			# C. Constructor (Panels)
			panels = UIPanel[]
			for view in page.views # view is ViewSpec
				push!(panels, build_panel!(uifig, view))
			end

			# D. Constructor (Decorations)
			widgets_s = make_widgets(S, ctx, page, panels)
			widgets_dict = build_toolbar!(uifig, widgets_s, ctx)

			build_statusbar!(uifig, Val(:status), ctx)
			build_legend!(uifig, Val(:legend), panels)

			# E. Assembly
			assem = UIPlot(S, ctx, page, uifig, panels, widgets_dict)
			push!(assemblies, assem)

			if display
				display!(ctx, assem)
			end
		end
	end
	return assemblies
end

# -------------------------
# The Constructors: build_ / build_!
# -------------------------

# Helper to find spec for a slot name (to retrieve attrs)
function get_slot_spec(uifig::UIFigure, name::Symbol)
	idx = findfirst(s -> s.name == name, uifig.layoutspec.slots)
	return idx === nothing ? nothing : uifig.layoutspec.slots[idx]
end

function build_figure(ctx::UIContext, page::PageSpec, ls::UILayoutSpec)
	kw = page.kwargs
	safe_kw = (; (k=>v for (k, v) in pairs(kw) if k != :size && k != :resolution)...)

	fig = Makie.Figure(; size = page.size, safe_kw...)

	containers = Dict{Symbol, Makie.GridLayout}()
	slots = Dict{Symbol, Any}()

	containers[:root] = fig.layout

	# --- PHASE 1: Materialize Slots/Containers ---
	# Uses `s.layout` for Grid properties

	# 1a. Intermediate Containers (and Root configuration)
	for c in ls.containers
		if c.name == :root
			# SPECIAL CASE: Configuration for the main Figure layout
			# Apply gaps, alignmode (padding), etc.
			gl = containers[:root]
			for (k, v) in pairs(c.layout)
				# We use setproperty! or specific Makie functions for gaps
				if k == :rowgap
					Makie.rowgap!(gl, v)
				elseif k == :colgap
					Makie.colgap!(gl, v)
				else
					# alignmode, etc.
					setproperty!(gl, k, v)
				end
			end
		else
			# Standard nested container creation
			parent_gl = containers[c.parent]
			subgl = Makie.GridLayout(; c.layout...)
			parent_gl[c.at...] = subgl
			containers[c.name] = subgl
		end
	end

	# 1b. Slots (Terminals)
	for s in ls.slots
		parent_gl = containers[s.parent]
		subgl = Makie.GridLayout(; s.layout...)
		parent_gl[s.at...] = subgl
		slots[s.name] = subgl

		# --- DEBUG: VISUALIZE SLOTS ---
		# Makie.Box(parent_gl[s.at...], color = (:red, 0.2), strokewidth = 0)
	end

	# --- PHASE 2: Apply Sizing ---

	for (name, sizes) in ls.rowsizes
		gl = get(containers, name, nothing)
		gl === nothing && continue
		for (i, s) in enumerate(sizes)
			Makie.rowsize!(gl, i, s)
		end
	end

	for (name, sizes) in ls.colsizes
		gl = get(containers, name, nothing)
		gl === nothing && continue
		for (i, s) in enumerate(sizes)
			Makie.colsize!(gl, i, s)
		end
	end

	# Grid Shape Logic
	n = length(page.views)
	panel_shape = (1, 1)
	if n > 1
		nr = ceil(Int, sqrt(n))
		nc = ceil(Int, n / nr)
		panel_shape = (nr, nc)
	end

	return UIFigure(fig, ls, containers, slots, Ref(0), panel_shape)
end

function build_panel!(uifig::UIFigure, view::ViewSpec)
	target_gl = uifig.slots[:canvas]
	nr, nc = uifig.panelshape

	ax_pos = if nr > 1 || nc > 1
		uifig.cursor[] += 1
		k = uifig.cursor[]
		row = (k - 1) ÷ nc + 1
		col = (k - 1) % nc + 1
		target_gl[row, col]
	else
		target_gl[1, 1]
	end

	# Retrieve 'attrs' for content
	slot_spec = get_slot_spec(uifig, :canvas)
	slot_attrs = slot_spec !== nothing ? slot_spec.attrs : (;)

	ax = Makie.Axis(ax_pos;
		xlabel = something(view.xaxis.label, ""),
		ylabel = something(view.yaxis.label, ""),
		title  = view.title,
		xscale = (view.xaxis.scale == :log10) ? Makie.log10 : Makie.identity,
		yscale = (view.yaxis.scale == :log10) ? Makie.log10 : Makie.identity,
		slot_attrs...,
	)

	plots = Any[]
	for s in view.series # s is SeriesSpec
		append!(plots, draw!(ax, s))
	end

	return UIPanel(view, ax, plots)
end

function build_toolbar!(uifig::UIFigure, specs::Vector{UIWidgetSpec}, ctx::UIContext)
	dict = Dict{Symbol, Any}()
	haskey(uifig.slots, :toolbar) || return dict

	gl = uifig.slots[:toolbar]
	gl.halign = :left

	for (i, s) in enumerate(specs)
		# --- BUTTON ---
		if s isa UIButtonSpec
			lbl = (s.icon !== nothing) ? with_icon(s.icon; text = s.label) : s.label

			btn = Makie.Button(gl[1, i]; label = lbl, s.attrs...)
			dict[Symbol(:btn_, i)] = btn

			Makie.on(btn.clicks) do _
				Base.@async begin
					try
						s.action(ctx, uifig, btn)
					catch e
						@error "Widget error" exception=(e, catch_backtrace())
						action_set_status!(ctx, "Error: $(e)")
					end
				end
			end

			# --- TOGGLE ---
		elseif s isa UIToggleSpec
			# Container for [Label | Toggle] to keep them grouped in the toolbar slot
			sub = Makie.GridLayout(gl[1, i])

			# 1. Label
			Makie.Label(sub[1, 1], s.label, halign = :right)

			# 2. Toggle
			tgl = Makie.Toggle(sub[1, 2]; active = s.active, s.attrs...)
			dict[Symbol(:tgl_, i)] = tgl

			Makie.on(tgl.active) do val
				Base.@async begin
					try
						s.action(ctx, uifig, val)
					catch e
						@error "Widget error" exception=(e, catch_backtrace())
						action_set_status!(ctx, "Error: $(e)")
					end
				end
			end

			# Tweak subgrid spacing
			Makie.colgap!(sub, 4)
		end
	end
	return dict
end

function build_statusbar!(uifig::UIFigure, ::Val{:status}, ctx::UIContext)
	haskey(uifig.slots, :status) || return
	gl = uifig.slots[:status]
	txt = (ctx.status !== nothing) ? ctx.status : Makie.Observable("")
	Makie.Label(gl[1, 1], txt, halign = :left, fontsize = 12)
end

function build_legend!(uifig::UIFigure, ::Val{:legend}, panels::Vector{UIPanel})
	haskey(uifig.slots, :legend) || return

	seen = Set{String}()
	elements = Any[]
	labels = String[]

	for p in panels, plt in p.plots
		if hasproperty(plt, :label)
			lbl = plt.label[]
			if lbl !== nothing && !isempty(lbl) && !(lbl in seen)
				push!(elements, plt)
				push!(labels, lbl)
				push!(seen, lbl)
			end
		end
	end

	if !isempty(elements)
		slot_spec = get_slot_spec(uifig, :legend)
		slot_attrs = slot_spec !== nothing ? slot_spec.attrs : (;)

		Makie.Legend(uifig.slots[:legend][1, 1], elements, labels; slot_attrs...)
	end
end

function display!(ctx::UIContext, assem::UIPlot)
	if ctx.interactive && ctx.window !== nothing
		display(ctx.window, assem.uifig.figure)
	else
		BackendHandler.renderfig(assem.uifig.figure)
	end
end