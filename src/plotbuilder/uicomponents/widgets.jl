# -------------------------
# Traits (Declarative)
# -------------------------

function controls_default(::Type{S}, ctx::UIContext, page, panels) where {S}
	!ctx.interactive && return UIWidgetSpec[]

	return UIWidgetSpec[
		UIButtonSpec(
			"",
			MI_REFRESH,
			(c, a, o) -> action_refresh(a),
		),
		UIButtonSpec(
			"",
			MI_SAVE,
			(c, a, o) -> action_export_svg!(c, a, "plot_export.svg"),
		),
	]
end

function controls_custom(::Type{S}, ctx, page, panels) where {S}
	return UIWidgetSpec[]
end

# -------------------------
# The Architect: make_widgets
# -------------------------

function make_widgets(::Type{S}, ctx::UIContext, page, panels) where {S}
	w = controls_default(S, ctx, page, panels)
	append!(w, controls_custom(S, ctx, page, panels))
	return w
end