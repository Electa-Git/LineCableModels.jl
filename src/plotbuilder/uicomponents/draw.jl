# -------------------------
# The Artist: draw!
# -------------------------

function draw!(axis, s::SeriesSpec; kwargs...)
	draw!(Val(s.kind), axis, s; kwargs...)
end

function draw!(::Val{kind}, axis, s::SeriesSpec; kwargs...) where {kind}
	@warn "Unsupported plot kind :$kind"
	return Any[]
end

function draw!(::Val{:line}, axis, s::SeriesSpec; kwargs...)
	(s.xdata === nothing || s.ydata === nothing) && return Any[]

	plots = Any[]
	if s.ydata isa AbstractMatrix
		for k in 1:size(s.ydata, 2)
			# PBSeries/SeriesSpec only supports one label.
			# We label the first trace for the legend.
			lbl = (k==1) ? s.label : nothing
			p = Makie.lines!(axis, s.xdata, view(s.ydata, :, k); label = lbl, kwargs...)
			push!(plots, p)
		end
	else
		p = Makie.lines!(axis, s.xdata, s.ydata; label = s.label, kwargs...)
		push!(plots, p)
	end
	return plots
end

function draw!(::Val{:scatter}, axis, s::SeriesSpec; kwargs...)
	(s.xdata === nothing || s.ydata === nothing) && return Any[]

	plots = Any[]
	if s.ydata isa AbstractMatrix
		for k in 1:size(s.ydata, 2)
			lbl = (k==1) ? s.label : nothing
			p = Makie.scatter!(axis, s.xdata, view(s.ydata, :, k); label = lbl, kwargs...)
			push!(plots, p)
		end
	else
		p = Makie.scatter!(axis, s.xdata, s.ydata; label = s.label, kwargs...)
		push!(plots, p)
	end
	return plots
end

function draw!(::Val{:heatmap}, axis, s::SeriesSpec; kwargs...)
	(s.xdata === nothing || s.ydata === nothing || s.zdata === nothing) && return Any[]
	p = Makie.heatmap!(axis, s.xdata, s.ydata, s.zdata; kwargs...)
	return Any[p]
end