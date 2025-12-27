function action_set_status!(ctx::UIContext, msg::AbstractString)
	if ctx.status !== nothing
		ctx.status[] = String(msg)
	end
	return nothing
end

function action_refresh(axis)
	Makie.autolimits!(axis)
	return nothing
end

"""
	action_export_svg!(ctx, assembly, path; kwargs...)

Required semantics:
- switch to Cairo
- re-render the same PB payload under Cairo
- save SVG
- restore backend
"""
function action_export_svg!(
	ctx::UIContext,
	assem::PlotAssembly,
	path::AbstractString;
	kwargs...,
)
	BackendHandler.with_backend(:cairo) do
		r = RenderSpec(assem.spec, PageSpec[assem.pbfig])
		a2 = first(render(r; backend = :cairo, display = false))
		Makie.save(path, a2.uifig.figure; kwargs...)
	end
	action_set_status!(ctx, "Saved SVG to $(path)")
	return nothing
end
