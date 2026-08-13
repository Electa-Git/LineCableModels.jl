function action_set_status!(ctx::UIContext, msg::AbstractString)
	if ctx.status !== nothing
		ctx.status[] = String(msg)
	end
	return nothing
end

# Adapting signature to match build! call (ctx, uifig, btn) -> needs mapping to Assembly?
# Ideally, we pass the Assembly. 
# We can rely on the fact that UIPlot is the container.
# Let's assume the action receives (ctx, assem::UIPlot, widget) 
# The build! loop in pipeline.jl needs to wrap this. 
# Correction: In pipeline.jl, we can't fully bind UIPlot because it's being built.
# However, `uifig` contains everything graphical.
# `action_refresh` needs access to panels. `uifig` does NOT store panels directly (UIPlot does).
# We should store panels in uifig or return to pipeline to bind actions AFTER UIPlot creation.
# Let's fix pipeline.jl logic in the next iteration or use a workaround here.
# Workaround: action_refresh takes `uifig` and assumes it can find axes? 
# No, `uifig` has `containers`. We can iterate `uifig.slots[:canvas].content`.

function action_refresh(uifig::UIFigure)
	# Iterate over axes in the canvas slot
	canvas = uifig.slots[:canvas]
	for content in canvas.content
		if content.content isa Makie.Axis
			Makie.autolimits!(content.content)
		end
	end
	return nothing
end

# Signature overload for compatibility if called with UIPlot
action_refresh(assem::UIPlot) = action_refresh(assem.uifig)

function action_export_svg!(ctx::UIContext, assem_or_fig, path::AbstractString)
	# We need the spec and page to re-render.
	# If we only have uifig, we are stuck.
	# The widgets need the UIPlot. 
	# FIX: The widgets must be wired up AFTER UIPlot is created in pipeline.jl.
	# See note in pipeline.jl.
	error("Export requires full UIPlot assembly context.")
end

function action_export_svg!(ctx::UIContext, assem::UIPlot, path::AbstractString)
	action_set_status!(ctx, "Exporting to $path...")

	BackendHandler.with_backend(:cairo) do
		# THEME SWITCHING: Force interactive=false for export style
		export_theme = make_theme(ctx; interactive = false)

		Makie.with_theme(export_theme) do
			# Reconstruct RenderSpec from the UIPlot data
			r = RenderSpec(assem.spec, PageSpec[assem.page])

			# Render headless
			new_assems = render(r; backend = :cairo, display = false)

			if !isempty(new_assems)
				target_fig = new_assems[1].uifig.figure
				Makie.save(path, target_fig)
			end
		end
	end

	action_set_status!(ctx, "Saved SVG to $path")
	return nothing
end
