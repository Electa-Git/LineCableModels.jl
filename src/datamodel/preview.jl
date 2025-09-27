using Makie, Colors
using Printf
using Dates
using Statistics


# _is_interactive_backend() = nameof(Makie.current_backend()) in (:GLMakie, :WGLMakie)
_is_interactive_backend() = current_backend_symbol() in (:gl, :wgl)
_is_static_backend() = current_backend_symbol() == :cairo
_is_gl_backend() = current_backend_symbol() == :gl


# finite & nonnegative
_valid_finite(x, y) = isfinite(x) && isfinite(y)


# Tunables (bands & palettes)
# ----------------------------
const RHO_MIN       = 1e-9      # for legend floor
const RHO_METAL_MAX = 1e-6
const RHO_SEMIMETAL = 1e-4
const RHO_SEMI_MAX  = 1e3
const RHO_LEAKY_MAX = 1e8
const RHO_MAX       = 1e10      # for legend ceiling

const METAL_GRADIENT = [
	RGB(0.92, 0.90, 0.86), # warm-silver (copper-ish)
	RGB(0.89, 0.89, 0.89), # neutral silver
	RGB(0.86, 0.89, 0.92), # cool-silver (aluminium-ish)
	RGB(0.70, 0.72, 0.75),  # slightly darker metal
]

const SEMIMETAL_GRADIENT = [
	RGB(0.70, 0.72, 0.75),   # gray
	RGB(0.80, 0.75, 0.65),    # sand/bronze hint
]

const SEMICON_GRADIENT = [
	RGB(1.00, 0.83, 0.40),   # light amber
	RGB(0.85, 0.55, 0.18),    # dark amber-brown
]

const LEAKY_GRADIENT = [
	RGB(0.42, 0.55, 0.15),   # olive/earthy
	RGB(0.13, 0.13, 0.13),    # charcoal
]

const INSULATOR_GRADIENT = [
	RGB(0.07, 0.07, 0.07),   # near-black (keep >0 so overlays remain visible)
	RGB(0.00, 0.00, 0.00),
]

# Overlays
const MU_OVERLAY_GRADIENT  = [RGB(0.20, 0.50, 0.95), RGB(0.56, 0.00, 0.91)]  # blue → indigo
const EPS_OVERLAY_GRADIENT = [RGB(0.00, 0.85, 0.70), RGB(0.00, 0.55, 0.90)]  # teal → cyan


# Linear interpolation across a list of colors in [0,1]
# robust gradient (no reinterpret)
_interpolate_gradient(colors::Vector{<:Colorant}, t::Real) = begin
	n = length(colors);
	n >= 2 || throw(ArgumentError("Need ≥ 2 colors"))
	tc = clamp(Float64(t), 0, 1)
	x  = tc * (n - 1)
	i  = clamp(floor(Int, x) + 1, 1, n - 1)
	f  = x - (i - 1)
	c1 = RGB(colors[i]);
	c2 = RGB(colors[i+1])
	RGB(
		(1 - f) * red(c1) + f * red(c2),
		(1 - f) * green(c1) + f * green(c2),
		(1 - f) * blue(c1) + f * blue(c2),
	)
end

# Log normalization helper: map v∈[a,b] (log10) → t∈[0,1]
_lognorm(v, a, b) = begin
	va = clamp(v, min(a, b), max(a, b))
	(log10(va) - log10(a)) / (log10(b) - log10(a))
end

_overlay(a::Colors.RGBA, b::Colors.RGBA) = begin
	a1, a2 = alpha(a), alpha(b)
	out_a = a2 + a1*(1 - a2)
	out_a == 0 && return Colors.RGBA(0, 0, 0, 0)
	r = (red(b)*a2 + red(a)*a1*(1 - a2)) / out_a
	g = (green(b)*a2 + green(a)*a1*(1 - a2)) / out_a
	b_ = (blue(b)*a2 + blue(a)*a1*(1 - a2)) / out_a
	Colors.RGBA(r, g, b_, out_a)
end

# Clamp lightness to keep overlays visible on "black"
function _ensure_min_lightness(c::RGB, Lmin::Float64 = 0.07)
	hsl = HSL(c)
	L = max(hsl.l, Lmin)
	rgb = RGB(HSL(hsl.h, hsl.s, L))
	return rgb
end

# ----------------------------
# Base color controlled by ρ
# ----------------------------
function _base_color_from_rho(ρ::Real)::RGB
	if !isfinite(ρ)
		return INSULATOR_GRADIENT[end]
	elseif ρ ≤ RHO_METAL_MAX
		t = _lognorm(ρ, 1e-8, RHO_METAL_MAX)
		return _interpolate_gradient(METAL_GRADIENT, t)
	elseif ρ ≤ RHO_SEMIMETAL
		t = _lognorm(ρ, RHO_METAL_MAX, RHO_SEMIMETAL)
		return _interpolate_gradient(SEMIMETAL_GRADIENT, t)
	elseif ρ ≤ RHO_SEMI_MAX
		t = _lognorm(ρ, RHO_SEMIMETAL, RHO_SEMI_MAX)
		return _interpolate_gradient(SEMICON_GRADIENT, t)
	elseif ρ ≤ RHO_LEAKY_MAX
		t = _lognorm(ρ, RHO_SEMI_MAX, RHO_LEAKY_MAX)
		return _interpolate_gradient(LEAKY_GRADIENT, t)
	else
		t = _lognorm(min(ρ, RHO_MAX), RHO_LEAKY_MAX, RHO_MAX)
		return _ensure_min_lightness(_interpolate_gradient(INSULATOR_GRADIENT, t), 0.07)
	end
end

# ----------------------------
# Overlays (μr & εr)
# ----------------------------
# μr in [1, 300] → alpha up to ~0.5, stronger on dark bases
function _mu_overlay(base::RGB, μr::Real)::Colors.RGBA
	μn = clamp((_lognorm(max(μr, 1.0), 1.0, 300.0)), 0, 1)
	tint = _interpolate_gradient(MU_OVERLAY_GRADIENT, μn)
	L = HSL(base).l
	α = 0.50 * μn * (0.6 + 0.4*(1 - L)) # reduce on bright silver, boost on dark
	Colors.RGBA(tint.r, tint.g, tint.b, α)
end

# εr in [1, 1000] → alpha up to ~0.6 on insulators, ~0.2 on metals
function _eps_overlay(base::RGB, εr::Real, ρ::Real)::Colors.RGBA
	εn = clamp((_lognorm(max(εr, 1.0), 1.0, 1000.0)), 0, 1)
	tint = _interpolate_gradient(EPS_OVERLAY_GRADIENT, εn)
	# weight more if it's an insulator/leaky (so it shows on dark)
	band_weight = ρ > RHO_SEMI_MAX ? 1.0 : (ρ > RHO_METAL_MAX ? 0.6 : 0.35)
	L = HSL(base).l
	α = (0.20 + 0.40*band_weight) * εn * (0.55 + 0.45*(1 - L))
	Colors.RGBA(tint.r, tint.g, tint.b, α)
end

"""
	get_material_color_makie(material_props; mu_scale=1.0, eps_scale=1.0)

Piecewise ρ→base color (metals→silver, semiconductors→amber, etc.) with
blue/purple magnetic overlay (μr) and teal/cyan permittivity overlay (εr).
`mu_scale` and `eps_scale` scale overlay strength (1.0 = default).
"""
function get_material_color_makie(material_props; mu_scale = 1.0, eps_scale = 1.0)
	ρ  = to_nominal(material_props.rho)
	εr = to_nominal(material_props.eps_r)
	μr = to_nominal(material_props.mu_r)

	base = _base_color_from_rho(ρ) |> c -> _ensure_min_lightness(c, 0.07)

	# Compose overlays
	mu  = _mu_overlay(base, μr);
	mu  = Colors.RGBA(mu.r, mu.g, mu.b, clamp(alpha(mu)*mu_scale, 0, 1))
	eps = _eps_overlay(base, εr, ρ);
	eps = Colors.RGBA(eps.r, eps.g, eps.b, clamp(alpha(eps)*eps_scale, 0, 1))

	out = _overlay(Colors.RGBA(base.r, base.g, base.b, 1.0), mu)
	out = _overlay(out, eps)
	return out
end

function show_material_scale(; size = (800, 400), backend = nothing)
	# if backend !== nothing
	# 	_use_makie_backend(backend)
	# end
	ensure_backend!(backend === nothing ? :cairo : backend)

	fig = Figure(size = size)

	# sampling density for smooth bars
	N = 1024

	# --- ρ colorbar (log scale by ticks/limits) -------------------------------
	ρmin_log, ρmax_log = log10(RHO_MIN), log10(RHO_MAX)
	# sample uniformly in log(ρ) so the bar matches your piecewise mapping
	cm_ρ = begin
		cols = Vector{RGBA}(undef, N)
		for i in 1:N
			t = (i - 1) / (N - 1)
			ρ = 10^(ρmin_log + t * (ρmax_log - ρmin_log))
			c = _base_color_from_rho(ρ)
			cols[i] = RGBA(c.r, c.g, c.b, 1.0)
		end
		cols
	end

	cb_ρ = Colorbar(fig[1, 1];
		colormap = cm_ρ,
		limits   = (ρmin_log, ρmax_log),   # we encode log(ρ) in limits/ticks
		vertical = false,
		label    = "Base color by resistivity ρ [Ω·m] (log scale)",
	)

	# label ticks at meaningful boundaries
	edges = [RHO_MIN, 1e-8, 1e-7, RHO_METAL_MAX, RHO_SEMIMETAL, RHO_SEMI_MAX,
		1e4, 1e6, RHO_LEAKY_MAX, RHO_MAX]
	cb_ρ.ticks = (log10.(edges), string.(edges))

	# --- μr overlay colorbar (blue→indigo on mid-gray) ------------------------
	μmin, μmax = 1.0, 300.0
	base_mid = RGB(0.5, 0.5, 0.5)
	cm_μ = begin
		cols = Vector{RGBA}(undef, N)
		for i in 1:N
			t = (i - 1) / (N - 1)
			μ = 10^(log10(μmin) + t * (log10(μmax) - log10(μmin)))
			o = _mu_overlay(base_mid, μ)
			out = _overlay(RGBA(base_mid.r, base_mid.g, base_mid.b, 1.0), o)
			cols[i] = out
		end
		cols
	end

	cb_μ = Colorbar(fig[2, 1];
		colormap = cm_μ,
		limits   = (μmin, μmax),
		vertical = false,
		label    = "Magnetic overlay μᵣ (blue→indigo)",
	)
	cb_μ.ticks = (
		[1, 2, 5, 10, 20, 50, 100, 200, 300],
		string.([1, 2, 5, 10, 20, 50, 100, 200, 300]),
	)

	# --- εr overlay colorbar (teal→cyan on dark base) -------------------------
	εmin, εmax = 1.0, 1000.0
	base_dark = RGB(0.10, 0.10, 0.10)
	cm_ε = begin
		cols = Vector{RGBA}(undef, N)
		for i in 1:N
			t = (i - 1) / (N - 1)
			ε = 10^(log10(εmin) + t * (log10(εmax) - log10(εmin)))
			o = _eps_overlay(base_dark, ε, RHO_MAX + 1)  # treat as strong insulator
			out = _overlay(RGBA(base_dark.r, base_dark.g, base_dark.b, 1.0), o)
			cols[i] = out
		end
		cols
	end

	cb_ε = Colorbar(fig[3, 1];
		colormap = cm_ε,
		limits   = (εmin, εmax),
		vertical = false,
		label    = "Permittivity overlay εᵣ (teal→cyan)",
	)
	cb_ε.ticks = ([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000],
		string.([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]))

	renderfig(fig)
	return fig
end


#################################
# Geometry helpers (polygons)   #
#################################
# polygons (Float32 points; filter non-finite)
function _annulus_poly(rin::Real, rex::Real, x0::Real, y0::Real; N::Int = 256)
	N ≥ 32 || throw(ArgumentError("N too small for a smooth annulus"))
	θo = range(0, 2π; length = N);
	θi = reverse(θo)
	xo = x0 .+ rex .* cos.(θo);
	yo = y0 .+ rex .* sin.(θo)
	xi = x0 .+ rin .* cos.(θi);
	yi = y0 .+ rin .* sin.(θi)
	px = vcat(xo, xi, xo[1]);
	py = vcat(yo, yi, yo[1])
	pts = Makie.Point2f.(px, py)
	filter(p -> _valid_finite(p[1], p[2]), pts)
end

function _circle_poly(r::Real, x0::Real, y0::Real; N::Int = 128)
	θ = range(0, 2π; length = N)
	x = x0 .+ r .* cos.(θ);
	y = y0 .+ r .* sin.(θ)
	pts = Makie.Point2f.(vcat(x, x[1]), vcat(y, y[1]))
	filter(p -> _valid_finite(p[1], p[2]), pts)
end

#############################
# Layer -> Makie primitives #
#############################
function _plot_layer_makie!(ax, layer, label::String;
	x0::Real = 0.0, y0::Real = 0.0, display_legend::Bool = true,
	legend_sink::Union{Nothing, Tuple} = nothing,
)

	if layer isa WireArray
		rwire = to_nominal(layer.radius_wire)
		nW = layer.num_wires
		lay_r = nW == 1 ? 0.0 : to_nominal(layer.radius_in)
		color = get_material_color_makie(layer.material_props)

		coords = calc_wirearray_coords(nW, rwire, to_nominal(lay_r), C = (x0, y0))

		plots = Any[]
		handle = nothing
		for (i, (x, y)) in enumerate(coords)
			poly = Makie.poly!(ax, _circle_poly(rwire, x, y);
				color = color,
				strokecolor = :black,
				strokewidth = 0.5,
				label = (i==1 && display_legend) ? label : "")
			push!(plots, poly)
			if i==1 && display_legend
				handle = poly
			end
		end

		# Legend sink: push one entry per layer. If sink has 3rd slot, store the group.
		if legend_sink !== nothing && display_legend && handle !== nothing
			push!(legend_sink[1], handle)
			push!(legend_sink[2], label)
			if length(legend_sink) >= 3
				push!(legend_sink[3], plots)   # group = all wires in this layer
			end
			if length(legend_sink) >= 4
				push!(legend_sink[4], to_nominal(layer.material_props.rho))  # <-- rho key
			end
		end
		return plots
	end

	if layer isa Strip || layer isa Tubular || layer isa Semicon || layer isa Insulator
		rin = to_nominal(layer.radius_in)
		rex = to_nominal(layer.radius_ext)
		color = get_material_color_makie(layer.material_props)

		poly = Makie.poly!(ax, _annulus_poly(rin, rex, x0, y0);
			color = color,
			label = display_legend ? label : "")

		if legend_sink !== nothing && display_legend
			push!(legend_sink[1], poly)
			push!(legend_sink[2], label)
			if length(legend_sink) >= 3
				push!(legend_sink[3], [poly])
			end
			if length(legend_sink) >= 4
				push!(legend_sink[4], NaN)      # not a wirearray
			end
		end
		return (poly,)
	end

	if layer isa ConductorGroup
		plots = Any[]
		first_label = true
		for sub in layer.layers
			append!(
				plots,
				_plot_layer_makie!(ax, sub,
					first_label ? lowercase(string(nameof(typeof(layer)))) : "";
					x0 = x0, y0 = y0, display_legend = display_legend,
					legend_sink = legend_sink),
			)
			first_label = false
		end
		return plots
	end

	if layer isa Sector
		vertices = layer.vertices
		# Convert vertices to Makie.Point2f format with offset
		makie_points = [Makie.Point2f(v[1] + x0, v[2] + y0) for v in vertices]
		# Ensure polygon is closed by adding first point at the end if needed
		if length(makie_points) > 0 && makie_points[1] != makie_points[end]
			push!(makie_points, makie_points[1])
		end
		
		color = get_material_color_makie(layer.material_props)

		poly = Makie.poly!(ax, makie_points;
			color = color,
			strokecolor = :black,
			strokewidth = 0.5,
			label = display_legend ? label : "")

		if legend_sink !== nothing && display_legend
			push!(legend_sink[1], poly)
			push!(legend_sink[2], label)
			if length(legend_sink) >= 3
				push!(legend_sink[3], [poly])
			end
			if length(legend_sink) >= 4
				push!(legend_sink[4], NaN)      # not a wirearray
			end
		end
		return (poly,)
	end

	if layer isa SectorInsulator
		outer_vertices = [(v[1] + x0, v[2] + y0) for v in layer.outer_vertices]
		# Convert to Makie.Point2f format
		outer_points = [Makie.Point2f(v[1], v[2]) for v in outer_vertices]
		# Ensure polygon is closed
		if length(outer_points) > 0 && outer_points[1] != outer_points[end]
			push!(outer_points, outer_points[1])
		end

		# (Not used for now) The inner boundary is the conductor's vertices. It must be reversed for the hole to be drawn correctly.
		inner_vertices = [(v[1] + x0, v[2] + y0) for v in layer.inner_sector.vertices]
		inner_points = [Makie.Point2f(v[1], v[2]) for v in inner_vertices]
		# Ensure inner polygon is closed
		if length(inner_points) > 0 && inner_points[1] != inner_points[end]
			push!(inner_points, inner_points[1])
		end
		color = get_material_color_makie(layer.material_props)
		# Create a shape with a hole by passing the outer boundary and holes as a vector of vectors
		polygon_with_hole = Makie.Polygon(outer_points, [inner_points])
		poly = Makie.poly!(ax, polygon_with_hole;
			color = color,
			strokecolor = :black,
			strokewidth = 0.5,
			label = display_legend ? label : "")

		if legend_sink !== nothing && display_legend
			push!(legend_sink[1], poly)
			push!(legend_sink[2], label)
			if length(legend_sink) >= 3
				push!(legend_sink[3], [poly])
			end
			if length(legend_sink) >= 4
				push!(legend_sink[4], NaN)      # not a wirearray
			end
		end
		return (poly,)
	end

	@warn "Unknown layer type $(typeof(layer)); skipping"
	return ()
end

function apply_default_theme!()
	bg = _is_static_backend() ? :white : :gray90
	set_theme!(backgroundcolor = bg, fonts = (; icons = ICON_TTF))
end



###############################################
# CableDesign cross-section (Makie version)   #
###############################################
function preview(design::CableDesign;
	x_offset::Real = 0.0,
	y_offset::Real = 0.0,
	backend::Union{Nothing, Symbol} = nothing,
	size::Tuple{Int, Int} = (800, 600),
	display_plot::Bool = true,
	display_legend::Bool = true,
	display_id::Bool = false,
	axis = nothing,
	legend_sink::Union{Nothing, Tuple{Vector{Any}, Vector{String}}} = nothing,
	display_colorbars::Bool = true,
	side_frac::Real = 0.26,     # ~26% right column
)

	ensure_backend!(backend)

	# backgroundcolor = (_is_static_backend() ? :white : :gray90)
	# set_theme!(backgroundcolor = backgroundcolor)
	apply_default_theme!()

	fig =
		isnothing(axis) ? Makie.Figure(size = size, figure_padding = (10, 10, 10, 10)) :
		nothing

	# ── 2 columns: left = main axis, right = container (button + legend + bars)
	local ax
	local side
	if isnothing(axis)
		ax   = Makie.Axis(fig[1, 1], aspect = Makie.DataAspect())
		side = fig[1, 2] = Makie.GridLayout()          # single container on the right
		Makie.colsize!(fig.layout, 1, Makie.Relative(1 - side_frac))
		Makie.colsize!(fig.layout, 2, Makie.Relative(side_frac))
		Makie.rowsize!(fig.layout, 1, Makie.Relative(1.0))

		ax.xlabel = "y [m]"
		ax.ylabel = "z [m]"

		ax.title =
			display_id ? "Cable design preview: $(design.cable_id)" :
			"Cable design preview"

		avail_w = size[1] * (1 - side_frac)
		avail_h = size[2]
		s = floor(Int, min(avail_w, avail_h)*0.9)
		Makie.colsize!(fig.layout, 1, Makie.Fixed(s))
		Makie.rowsize!(fig.layout, 1, Makie.Fixed(s))
	else
		ax = axis
		side = nothing
	end

	# legend sink
	local own_legend = false
	local sink = legend_sink
	if sink === nothing && display_legend
		sink = (Any[], String[], Vector{Vector{Any}}(), Float64[])   # handles, labels, groups, rho_keys
		own_legend = true
	end

	let r = try
			to_nominal(design.components[end].insulator_group.radius_ext)
		catch
			NaN
		end
		if isfinite(r) && r > 0
			Makie.poly!(ax, _circle_poly(r, x_offset, y_offset);
				color = :white,
				strokecolor = :transparent)
		end
	end

	# draw layers
	for comp in design.components
		for layer in comp.conductor_group.layers
			_plot_layer_makie!(ax, layer, lowercase(string(nameof(typeof(layer))));
				x0 = x_offset, y0 = y_offset,
				display_legend = display_legend, legend_sink = sink)
		end
		for layer in comp.insulator_group.layers
			_plot_layer_makie!(ax, layer, lowercase(string(nameof(typeof(layer))));
				x0 = x_offset, y0 = y_offset,
				display_legend = display_legend, legend_sink = sink)
		end
	end


	# Right column: stack button, legend, colorbars
	if isnothing(axis)
		row_idx = 1

		if _is_interactive_backend()
			# Reset button at top
			_add_reset_button!(side[row_idx, 1], ax, fig)
			row_idx += 1
			_add_save_svg_button!(
				side[row_idx, 1], design;
				display_id = display_id,
				display_legend = display_legend,
				display_colorbars = display_colorbars,
				side_frac = side_frac,
				size = size,
				base = design.cable_id,
			)
			row_idx += 1
		end

		# Legend (optional)
		if display_legend && own_legend
			handles = sink[1]
			labels  = sink[2]
			groups  = length(sink) >= 3 ? sink[3] : [[h] for h in handles]
			rhos    = length(sink) >= 4 ? sink[4] : fill(NaN, length(handles))

			# Merge consecutive wirearray entries with equal rho
			merged_handles = Any[]
			merged_labels  = String[]
			merged_groups  = Vector{Any}[]  # Vector{Vector{Any}}

			i = 1
			while i <= length(handles)
				h = handles[i];
				l = labels[i];
				g = groups[i];
				ρ = rhos[i]
				if l == "wirearray" && isfinite(ρ)
					j = i + 1
					merged_g = Vector{Any}(g)
					while j <= length(handles) &&
							  labels[j] == "wirearray" &&
							  isfinite(rhos[j]) &&
							  isapprox(ρ, rhos[j]; rtol = 1e-6, atol = 0.0)
						append!(merged_g, groups[j])
						j += 1
					end
					push!(merged_handles, h)        # keep first handle for the group
					push!(merged_labels, l)         # keep label "wirearray"
					push!(merged_groups, merged_g)  # all wires across merged layers
					i = j
				else
					push!(merged_handles, h)
					push!(merged_labels, l)
					push!(merged_groups, g)
					i += 1
				end
			end

			# Build legend with merged entries
			leg = Makie.Legend(
				side[row_idx, 1],
				merged_handles,
				merged_labels,
				padding = (6, 6, 6, 6),
				halign = :center,
				valign = :top,
			)

			# Clicking one entry toggles its whole merged group
			for (h, grp) in zip(merged_handles, merged_groups)
				Makie.on(h.visible) do v
					for p in grp
						p === h && continue
						p.visible[] = v
					end
				end
			end

			row_idx += 1
		end

		# Colorbars (optional)
		if display_colorbars
			# read actual ranges (helper you already have)
			ρmin, ρmax, μmin, μmax, εmin, εmax = _collect_material_ranges(design)

			cbgrid = side[row_idx, 1] = Makie.GridLayout()
			_build_colorbars!(cbgrid; ρmin, ρmax, μmin, μmax, εmin, εmax)

		end
	end

	if display_plot && isnothing(axis) && !is_in_testset()
		resize_to_layout!(fig)
		n = next_fignum()
		scr =
			_is_gl_backend() ?
			gl_screen("Fig. $(n) – CableDesign preview: $(design.cable_id)") :
			nothing
		if scr === nothing
			renderfig(fig)
		else
			display(scr, fig)
		end

	end
	return fig, ax
end

function preview(system::LineCableSystem;
	earth_model = nothing,
	zoom_factor = nothing,
	backend::Union{Nothing, Symbol} = nothing,
	size::Tuple{Int, Int} = (800, 600),
	display_plot::Bool = true,
	display_id::Bool = false,
	axis = nothing,
	display_legend::Bool = true,
	display_colorbars::Bool = true,
	side_frac::Real = 0.26,
)

	ensure_backend!(backend)
	# backgroundcolor = (_is_static_backend() ? :white : :gray90)

	# set_theme!(backgroundcolor = backgroundcolor)
	apply_default_theme!()

	fig =
		isnothing(axis) ? Makie.Figure(size = size, figure_padding = (10, 10, 10, 10)) :
		nothing

	# Layout: left = main axis, right = legend/colorbars (only if we own the axis)
	local ax
	local side
	if isnothing(axis)
		ax   = Makie.Axis(fig[1, 1], aspect = Makie.DataAspect())
		side = fig[1, 2] = Makie.GridLayout()
		Makie.colsize!(fig.layout, 1, Makie.Relative(1 - side_frac))
		Makie.colsize!(fig.layout, 2, Makie.Relative(side_frac))
		Makie.rowsize!(fig.layout, 1, Makie.Relative(1.0))

		ax.xlabel = "y [m]"
		ax.ylabel = "z [m]"

		ax.title =
			display_id ? "Cable system cross-section: $(system.system_id)" :
			"Cable system cross-section"

		# Make the plotting canvas square if we own the axis
		avail_w = size[1] * (1 - side_frac)
		avail_h = size[2]
		s = floor(Int, min(avail_w, avail_h)*0.9)
		Makie.colsize!(fig.layout, 1, Makie.Fixed(s))
		Makie.rowsize!(fig.layout, 1, Makie.Fixed(s))
	else
		ax = axis
		side = nothing
	end

	# Air/earth interface
	Makie.hlines!(ax, [0.0], color = :black, linewidth = 1.5)

	# Compute barycentered, square view from cable bounding box
	x0s = Float64[to_nominal(c.horz) for c in system.cables]
	y0s = Float64[to_nominal(c.vert) for c in system.cables]
	radii = Float64[
		(comp = last(c.design_data.components);
			max(to_nominal(comp.conductor_group.radius_ext),
				to_nominal(comp.insulator_group.radius_ext)))
		for c in system.cables
	]
	cx = isempty(x0s) ? 0.0 : mean(x0s)
	cy = isempty(y0s) ? -1.0 : mean(y0s)
	x_min = isempty(x0s) ? -1.0 : minimum(x0s .- radii)
	x_max = isempty(x0s) ? 1.0 : maximum(x0s .+ radii)
	y_min = isempty(y0s) ? -1.0 : minimum(y0s .- radii)
	y_max = isempty(y0s) ? 1.0 : maximum(y0s .+ radii)
	half_x = max(x_max - cx, cx - x_min)
	half_y = max(y_max - cy, cy - y_min)
	base_halfspan = max(half_x, half_y)
	base_halfspan = base_halfspan > 0 ? base_halfspan : 1.0
	pad_factor = 1.05
	zf = zoom_factor === nothing ? 1.5 : Float64(zoom_factor)
	halfspan = base_halfspan * pad_factor * zf
	x_limits = (cx - halfspan, cx + halfspan)
	y_limits = (cy - halfspan, cy + halfspan)

	# Expanded fill extents beyond visible region
	BUFFER_FILL = 5.0
	x_fill =
		(x_limits[1] - 0.5*halfspan - BUFFER_FILL, x_limits[2] + 0.5*halfspan + BUFFER_FILL)
	y_fill_min = y_limits[1] - 0.5*halfspan - BUFFER_FILL

	# Build legend entries only for earth layers
	earth_handles = Any[]
	earth_labels  = String[]

	# Plot earth layers if provided and horizontal (vertical_layers == false)
	if !isnothing(earth_model) && getfield(earth_model, :vertical_layers) == false
		cumulative_depth = 0.0
		# Skip air layer (index 1). Iterate finite-thickness layers; stop on Inf.
		for (i, layer) in enumerate(earth_model.layers[2:end])
			# Compute color using the same material convention
			# Adapt EarthLayer base_* fields to material_props (rho, eps_r, mu_r)
			mat = (;
				rho = layer.base_rho_g,
				eps_r = layer.base_epsr_g,
				mu_r = layer.base_mur_g,
			)
			fillcol = get_material_color_makie(mat)
			# Slight transparency for fill
			fillcol = Makie.RGBA(fillcol.r, fillcol.g, fillcol.b, 0.25)

			if isinf(layer.t)
				# Semi-infinite: fill from current depth down to far below visible
				ytop = cumulative_depth  # bottom of previous finite layer
				ybot = y_fill_min        # push well below visible range
			else
				# Finite thickness: update cumulative and compute band extents
				t = to_nominal(layer.t)
				ytop = cumulative_depth
				ybot = cumulative_depth - t
				cumulative_depth = ybot
			end

			xs = (x_fill[1], x_fill[2], x_fill[2], x_fill[1])
			ys = (ytop, ytop, ybot, ybot)

			# Filled band and a colored interface line
			poly = Makie.poly!(
				ax,
				collect(Makie.Point2f.(xs, ys)),  # ensure a Vector, not a Tuple
				color = fillcol,
				strokecolor = :transparent,
				label = "",
			)
			Makie.hlines!(ax, [ybot], color = fillcol, linewidth = 1.0)

			if display_legend && isnothing(axis)
				push!(earth_handles, poly)
				push!(earth_labels, "Earth layer $(i)")
			end
		end
	end

	# Draw each cable onto the same axis (no legend for cable components)
	for cable in system.cables
		x0 = to_nominal(cable.horz)
		y0 = to_nominal(cable.vert)
		# Reuse the design-level preview on our axis
		preview(
			cable.design_data;
			x_offset = x0,
			y_offset = y0,
			backend = backend,
			size = size,
			display_plot = false,
			display_legend = false,
			axis = ax,
		)
	end

	# Set limits only when we own the axis (square extents)
	if isnothing(axis)
		Makie.xlims!(ax, x_limits...)
		Makie.ylims!(ax, y_limits...)
	end

	# Right-column: buttons, earth-only legend and optional colorbars
	if isnothing(axis)
		row_idx = 1
		if _is_interactive_backend()
			_add_reset_button!(side[row_idx, 1], ax, fig)
			row_idx += 1
			_add_save_svg_button!(
				side[row_idx, 1], system;
				earth_model = earth_model,
				zoom_factor = zoom_factor,
				display_legend = display_legend,
				display_colorbars = display_colorbars,
				side_frac = side_frac,
				display_id = display_id,
				size = size,
				base = system.system_id,
			)
			row_idx += 1
		end

		if display_legend && !isempty(earth_handles)
			Makie.Legend(
				side[row_idx, 1],
				earth_handles,
				earth_labels,
				padding = (6, 6, 6, 6),
				halign = :center,
				valign = :top,
			)
			row_idx += 1
		end

		if display_colorbars
			ρmin, ρmax, μmin, μmax, εmin, εmax = _collect_earth_ranges(earth_model)
			cbgrid = side[row_idx, 1] = Makie.GridLayout()
			_build_colorbars!(
				cbgrid;
				ρmin,
				ρmax,
				μmin,
				μmax,
				εmin,
				εmax,
				alpha_global = 0.25,
				showμminmax = false,
				showεminmax = false,
			)
		end
	end


	if display_plot && isnothing(axis) && !is_in_testset()
		resize_to_layout!(fig)
		n = next_fignum()
		scr =
			_is_gl_backend() ?
			gl_screen("Fig. $(n) – LineCableSystem preview: $(system.system_id)") :
			nothing
		if scr === nothing
			renderfig(fig)
		else
			display(scr, fig)
		end

	end

	return fig, ax
end

# Add a save-to-SVG button to a grid cell, re-rendering with Cairo backend
function _add_save_svg_button!(parent_cell, system;
	earth_model = nothing,
	zoom_factor = nothing,
	display_id::Bool,
	display_legend::Bool,
	display_colorbars::Bool,
	side_frac::Real,
	size::Tuple{Int, Int},
	base::String = "preview",
	save_dir::AbstractString = pwd(),
)
	btn = Makie.Button(
		parent_cell,
		label = with_icon(MI_SAVE; text = "Save SVG"),
		halign = :center,
		valign = :top,
		width = Makie.Auto(),
	)
	Makie.on(btn.clicks) do _
		@async begin
			orig_label = btn.label[]
			btn.label[] = "Saving…"
			orig_color = hasproperty(btn, :buttoncolor) ? btn.buttoncolor[] : nothing
			try
				# _use_makie_backend(:cairo)
				ensure_backend!(:cairo)
				if system isa CableDesign
					fig, _ = preview(
						system;
						display_legend = display_legend,
						display_colorbars = display_colorbars,
						display_plot = false,
						size = size,
						display_id = display_id,
						backend = :cairo,
						side_frac = side_frac,
					)
				elseif system isa LineCableSystem
					fig, _ = preview(
						system;
						earth_model = earth_model,
						zoom_factor = zoom_factor,
						display_id = display_id,
						display_legend = display_legend,
						display_colorbars = display_colorbars,
						display_plot = false,
						size = size,
						backend = :cairo,
						side_frac = side_frac,
					)
				end
				ts = Dates.format(Dates.now(), "yyyymmdd-HHMMSS")
				file = joinpath(save_dir, "$(base)_$ts.svg")
				Makie.save(file, fig)
				btn.label[] = "Saved ✓"
				@info "Saved figure to $(file)"
				hasproperty(btn, :buttoncolor) &&
					(btn.buttoncolor[] = Makie.RGBA(0.15, 0.65, 0.25, 1.0))
				sleep(1.2)
			catch e
				@error "Save failed: $(typeof(e)): $(e)"
				btn.label[] = "Failed ✗"
				hasproperty(btn, :buttoncolor) &&
					(btn.buttoncolor[] = Makie.RGBA(0.80, 0.20, 0.20, 1.0))
				sleep(1.6)
			finally
				if orig_color !== nothing
					btn.buttoncolor[] = orig_color
				end
				btn.label[] = orig_label
			end
		end
	end
	return btn
end

# Add a reset button to a grid cell, wired to reset axis limits
function _add_reset_button!(parent_cell, ax, fig)
	btn = Makie.Button(
		parent_cell,
		label = with_icon(MI_REFRESH; text = "Reset view"),
		halign = :center,
		valign = :top,
		width = Makie.Relative(1.0),
	)
	Makie.on(btn.clicks) do _
		reset_limits!(ax)
		# resize_to_layout!(fig)

	end
	return btn
end

function _build_colorbars!(cbgrid::Makie.GridLayout;
	ρmin::Real, ρmax::Real, μmin::Real, μmax::Real, εmin::Real, εmax::Real,
	cb_bar_h::Int = 12, alpha_global::Real = 1.0, showρminmax::Bool = true,
	showμminmax::Bool = true, showεminmax::Bool = true,
)
	Makie.colsize!(cbgrid, 1, Makie.Fixed(2))

	function _nice(x)
		axv = abs(x)
		axv == 0 && return "0"
		(axv ≥ 1e-3 && axv < 1e4) ? @sprintf("%.4g", x) : @sprintf("%.1e", x)
	end

	N = 256
	idx=1

	# ρ bar sampled in log-space between actual min/max
	if showρminmax
		cm_ρ = let cols = Vector{Makie.RGBA}(undef, N)
			lo, hi = log10(ρmin), log10(ρmax)
			for i in 1:N
				t = (i-1)/(N-1)
				ρ = 10^(lo + t*(hi - lo))
				c = _base_color_from_rho(ρ)
				cols[i] = Makie.RGBA(c.r, c.g, c.b, 1*alpha_global)
			end;
			cols
		end
		Makie.Label(cbgrid[idx, 1], L"\rho"; halign = :left, fontsize = 16)
		Makie.Colorbar(cbgrid[idx, 2]; colormap = cm_ρ, limits = (0.0, 1.0),
			vertical = false,
			ticks = ([0.0, 1.0], [_nice(ρmin), _nice(ρmax)]),
			labelvisible = false, height = cb_bar_h)
		idx += 1
	end

	# μr bar overlay on mid-gray
	if showμminmax
		base_mid = Makie.RGB(0.5, 0.5, 0.5)
		cm_μ = let cols = Vector{Makie.RGBA}(undef, N)
			lo, hi = log10(μmin), log10(μmax)
			for i in 1:N
				t = (i-1)/(N-1)
				μ = 10^(lo + t*(hi - lo))
				o = _mu_overlay(base_mid, μ)
				cols[i] = _overlay(
					Makie.RGBA(base_mid.r, base_mid.g, base_mid.b, 1),
					o*alpha_global,
				)
			end;
			cols
		end
		Makie.Label(cbgrid[idx, 1], L"\mu_{r}"; halign = :left, fontsize = 16)
		Makie.Colorbar(cbgrid[idx, 2]; colormap = cm_μ, limits = (0.0, 1.0),
			vertical = false,
			ticks = ([0.0, 1.0], [_nice(μmin), _nice(μmax)]),
			labelvisible = false, height = cb_bar_h)
		idx += 1
	end

	if showεminmax
		# εr bar overlay on dark
		base_dark = Makie.RGB(0.10, 0.10, 0.10)
		cm_ε = let cols = Vector{Makie.RGBA}(undef, N)
			lo, hi = log10(εmin), log10(εmax)
			for i in 1:N
				t = (i-1)/(N-1)
				ε = 10^(lo + t*(hi - lo))
				o = _eps_overlay(base_dark, ε, RHO_MAX + 1)
				cols[i] = _overlay(
					Makie.RGBA(base_dark.r, base_dark.g, base_dark.b, 1),
					o*alpha_global,
				)
			end;
			cols
		end
		Makie.Label(cbgrid[idx, 1], L"\varepsilon_{r}"; halign = :left, fontsize = 16)
		Makie.Colorbar(cbgrid[idx, 2]; colormap = cm_ε, limits = (0.0, 1.0),
			vertical = false,
			ticks = ([0.0, 1.0], [_nice(εmin), _nice(εmax)]),
			labelvisible = false, height = cb_bar_h)
	end

	return cbgrid
end

# collect actual property ranges from the design (finite values only)
function _collect_material_ranges(design::CableDesign)
	rhos  = Float64[]
	mus   = Float64[]
	epses = Float64[]

	_push_props!(layer) = begin
		ρ  = try
			to_nominal(layer.material_props.rho)
		catch
			NaN
		end
		μr = try
			to_nominal(layer.material_props.mu_r)
		catch
			NaN
		end
		εr = try
			to_nominal(layer.material_props.eps_r)
		catch
			NaN
		end
		isfinite(ρ) && push!(rhos, ρ)
		isfinite(μr) && push!(mus, μr)
		isfinite(εr) && push!(epses, εr)
		nothing
	end

	for comp in design.components
		for L in comp.conductor_group.layers
			if L isa ConductorGroup
				for s in L.layers
					_push_props!(s);
				end
			else
				_push_props!(L)
			end
		end
		for L in comp.insulator_group.layers
			_push_props!(L)
		end
	end

	ρmin = isempty(rhos) ? RHO_MIN : minimum(rhos)
	ρmax = isempty(rhos) ? RHO_MAX : maximum(rhos)
	μmin = isempty(mus) ? 1.0 : max(1.0, minimum(mus))
	μmax = isempty(mus) ? 300.0 : maximum(mus)
	εmin = isempty(epses) ? 1.0 : max(1.0, minimum(epses))
	εmax = isempty(epses) ? 1000.0 : maximum(epses)
	ρmax == ρmin && (ρmax = nextfloat(ρmax))
	μmax == μmin && (μmax += 1e-6)
	εmax == εmin && (εmax += 1e-6)

	return ρmin, ρmax, μmin, μmax, εmin, εmax
end

function _collect_earth_ranges(earth_model)
	rhos = Float64[];
	mus = Float64[];
	epses = Float64[]
	if !isnothing(earth_model)
		for layer in earth_model.layers[2:end]
			ρ = try
				to_nominal(layer.base_rho_g)
			catch
				NaN
			end
			μr = try
				to_nominal(layer.base_mur_g)
			catch
				NaN
			end
			εr = try
				to_nominal(layer.base_epsr_g)
			catch
				NaN
			end
			isfinite(ρ) && push!(rhos, ρ)
			isfinite(μr) && push!(mus, μr)
			isfinite(εr) && push!(epses, εr)
		end
	end
	ρmin = isempty(rhos) ? RHO_MIN : minimum(rhos)
	ρmax = isempty(rhos) ? RHO_MAX : maximum(rhos)
	μmin = isempty(mus) ? 1.0 : max(1.0, minimum(mus))
	μmax = isempty(mus) ? 300.0 : maximum(mus)
	εmin = isempty(epses) ? 1.0 : max(1.0, minimum(epses))
	εmax = isempty(epses) ? 1000.0 : maximum(epses)
	return ρmin, ρmax, μmin, μmax, εmin, εmax
end
