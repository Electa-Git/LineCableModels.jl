# -------------------------
# The Architect: make_layout
# -------------------------

function make_layout(::Val{layout}) where {layout}
	error("make_layout(Val(:$layout)) not implemented")
end

const _LEGEND_WIDTH = 150
const _TOOLBAR_H = 36
const _STATUS_H = 20

function make_layout(::Val{:windows})
	# 3x2 Grid

	slots = [
		# Toolbar
		# layout: Container behavior (height, padding)
		UISlotSpec(:toolbar, :root, (1, 1);
			layout = (;
				height = _TOOLBAR_H,
				tellheight = true,
				alignmode = Makie.Outside(5),
			),
		),

		# Canvas
		# layout: Standard plot area
		UISlotSpec(:canvas, :root, (2, 1);
			layout = (; alignmode = Makie.Inside()),
		),

		# Status
		UISlotSpec(:status, :root, (3, 1);
			layout = (;
				height = _STATUS_H,
				tellheight = true,
				alignmode = Makie.Outside(2),
			),
		),

		# Legend
		# layout: Grid behavior (padding from edge)
		# attrs:  Legend behavior (anchor to top)
		UISlotSpec(:legend, :root, (1:3, 2);
			layout = (; alignmode = Makie.Outside(5)),
			attrs  = (; valign = :top),
		),
	]

	rs = Dict(
		:root => Any[
			Makie.Fixed(_TOOLBAR_H),
			Makie.Relative(1.0),
			Makie.Fixed(_STATUS_H),
		],
	)

	cs = Dict(:root => Any[
		Makie.Relative(1.0),
		Makie.Fixed(_LEGEND_WIDTH),
	])

	return UILayoutSpec(:windows, UIContainerSpec[], slots, rs, cs)
end

make_layout(::Val{:grid}) = make_layout(Val(:windows))