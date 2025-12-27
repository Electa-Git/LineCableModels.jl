# Trait surface
function layout_spec(::Type{S}, ::Val{layout}) where {S, layout}
	return layout_spec(Val(layout))
end

# -------------------------
# Layout-only fallbacks
# -------------------------

# :windows -> toolbar / canvas / status (single view)
function layout_spec(::Val{:windows})
	containers = UIContainerSpec[]
	slots = UISlotSpec[
		UISlotSpec(:toolbar, :root, (1, 1); props = (; rowsize = Makie.Auto(0))),
		UISlotSpec(:canvas, :root, (2, 1); props = (; rowsize = Makie.Auto())),
		UISlotSpec(:status, :root, (3, 1); props = (; rowsize = Makie.Auto(0))),
	]
	return UILayoutSpec(:windows, containers, slots)
end

function layout_spec(::Val{:grid})
	containers = UIContainerSpec[]
	slots = UISlotSpec[
		UISlotSpec(:toolbar, :root, (1, 1); props = (; rowsize = Makie.Auto(0))),
		UISlotSpec(:panelgrid, :root, (2, 1); props = (; rowsize = Makie.Auto())),
		UISlotSpec(:status, :root, (3, 1); props = (; rowsize = Makie.Auto(0))),
	]
	return UILayoutSpec(:grid, containers, slots)
end

# Last resort
function layout_spec(::Val{layout}) where {layout}
	error("layout_spec(Val($(layout))) not implemented")
end
