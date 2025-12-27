function controls_default(
	ctx::UIContext,
	::Type{S},
	pbfig::PageSpec,
	panels::Vector{UIPanel},
) where {S}
	return Dict{Symbol, Any}()
end

function controls_extra(
	::Type{S},
	ctx::UIContext,
	pbfig::PageSpec,
	panels::Vector{UIPanel},
) where {S}
	return Dict{Symbol, Any}()
end

function build_widgets!(
	uifig::UIFigure,
	ctx::UIContext,
	::Type{S},
	pbfig::PageSpec,
	panels::Vector{UIPanel};
	kwargs...,
) where {S}
	w = controls_default(ctx, S, pbfig, panels)
	merge!(w, controls_extra(S, ctx, pbfig, panels))
	return w
end
