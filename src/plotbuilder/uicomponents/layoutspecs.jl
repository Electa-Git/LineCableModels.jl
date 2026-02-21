# -------------------------
# The Architect: make_layout
# -------------------------

function make_layout(::Val{layout}) where {layout}
	error("make_layout(Val(:$layout)) not implemented")
end

function make_layout(::Val{:single})
	# 1. ROOT CONFIGURATION (The global visual effect)
	# We define a container for :root to apply gaps and padding.
	root = UIContainerSpec(
		:root, nothing, (1, 1); # Position ignored for root
		layout = (;
			rowgap = GRID_ROW_GAP,
			colgap = LEGEND_GAP, # The gap between Canvas and Legend
			alignmode = Makie.Outside(FIG_PADDING...), # (L, R, B, T)
		),
	)

	# 2. SLOTS
	slots = [
		# Toolbar: Rigid height, internal spacing for buttons
		UISlotSpec(:toolbar, :root, (1, 1);
			layout = (;
				height = CTLBAR_HEIGHT,
				tellheight = true,
				colgap = CTLBAR_GAP # Spacing between buttons
			),
		),

		# Canvas: Takes available space
		UISlotSpec(:canvas, :root, (2, 1);
			layout = (;
				alignmode = Makie.Inside() # Standard plot behavior
			),
		),

		# Status: Rigid height
		UISlotSpec(:status, :root, (3, 1);
			layout = (;
				height = STATUSBAR_HEIGHT,
				tellheight = true,
			),
		),

		# Legend: Fixed width column
		UISlotSpec(:legend, :root, (1:3, 2);
			layout = (;
				width = LEGEND_WIDTH, # Enforce width at slot level too for safety
				tellwidth = true,
				alignmode = Makie.Inside(),
			),
			attrs = (; valign = :top),
		),
	]

	# 3. ROOT SIZING (Structural constraints)
	rs = Dict(
		:root => Any[
			Makie.Fixed(CTLBAR_HEIGHT),
			Makie.Relative(1.0),
			Makie.Fixed(STATUSBAR_HEIGHT),
		],
	)

	cs = Dict(:root => Any[
		Makie.Relative(1.0),
		Makie.Fixed(LEGEND_WIDTH),
	])

	return UILayoutSpec(:single, [root], slots, rs, cs)
end

make_layout(::Val{:grid}) = make_layout(Val(:single))