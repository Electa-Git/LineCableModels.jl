function draw!(axis, s::SeriesSpec; kwargs...)
	return draw!(Val(s.kind), axis, s; kwargs...)
end

function draw!(::Val{kind}, axis, s::SeriesSpec; kwargs...) where {kind}
	error("draw!(Val($(kind)), ...) not implemented")
end

# -------------------------
# Minimal primitives
# -------------------------

function draw!(::Val{:line}, axis, s::SeriesSpec; kwargs...)
	x = s.xdata
	y = s.ydata
	x === nothing && error("SeriesSpec(:line): xdata is nothing")
	y === nothing && error("SeriesSpec(:line): ydata is nothing")

	plots = Any[]
	if y isa AbstractVector
		push!(plots, Makie.lines!(axis, x, y; label = s.label, kwargs...))
	elseif y isa AbstractMatrix
		# columns as separate traces
		for k in 1:size(y, 2)
			lab = (k == 1) ? s.label : nothing
			push!(plots, Makie.lines!(axis, x, view(y, :, k); label = lab, kwargs...))
		end
	else
		error("SeriesSpec(:line): unsupported ydata type $(typeof(y))")
	end
	return plots
end

function draw!(::Val{:scatter}, axis, s::SeriesSpec; kwargs...)
	x = s.xdata
	y = s.ydata
	x === nothing && error("SeriesSpec(:scatter): xdata is nothing")
	y === nothing && error("SeriesSpec(:scatter): ydata is nothing")

	plots = Any[]
	if y isa AbstractVector
		push!(plots, Makie.scatter!(axis, x, y; label = s.label, kwargs...))
	elseif y isa AbstractMatrix
		for k in 1:size(y, 2)
			lab = (k == 1) ? s.label : nothing
			push!(plots, Makie.scatter!(axis, x, view(y, :, k); label = lab, kwargs...))
		end
	else
		error("SeriesSpec(:scatter): unsupported ydata type $(typeof(y))")
	end
	return plots
end

function draw!(::Val{:heatmap}, axis, s::SeriesSpec; kwargs...)
	x = s.xdata
	y = s.ydata
	z = s.zdata

	x === nothing && error("SeriesSpec(:heatmap): xdata is nothing")
	y === nothing && error("SeriesSpec(:heatmap): ydata is nothing")
	z === nothing && error("SeriesSpec(:heatmap): zdata is nothing")

	# Expect x::Vector, y::Vector, z::Matrix
	(y isa AbstractVector) ||
		error("SeriesSpec(:heatmap): expected ydata Vector, got $(typeof(y))")
	(z isa AbstractMatrix) ||
		error("SeriesSpec(:heatmap): expected zdata Matrix, got $(typeof(z))")

	return Any[Makie.heatmap!(axis, x, y, z; kwargs...)]
end
