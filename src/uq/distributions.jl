"""
Freedman–Diaconis rule to guesstimate number of bins for histogram

For samples x:

IQR = interquartile range = q75 - q25

bin width
h = 2 * IQR / N^(1/3)
number of bins ~ (max(x) - min(x)) / h
"""
function _auto_nbins(x::AbstractVector{<:Real};
	nbins_min::Int = 10,
	nbins_max::Int = 200,
)
	n = length(x)
	n == 0 && error("Empty sample set.")

	xs = sort(float.(x))
	xmin, xmax = xs[1], xs[end]
	span = xmax - xmin

	# degenerate span: all samples equal (or numerically so)
	if span <= 0 || !isfinite(span)
		return nbins_min
	end

	q25 = quantile(xs, 0.25)
	q75 = quantile(xs, 0.75)
	iqr = q75 - q25

	# iqr ~ 0 → data essentially degenerate → fallback
	if iqr <= 0 || !isfinite(iqr)
		return clamp(ceil(Int, sqrt(n)), nbins_min, nbins_max)
	end

	h = 2 * iqr / n^(1/3)

	# h tiny or broken → fallback
	if h <= 0 || !isfinite(h)
		return clamp(ceil(Int, sqrt(n)), nbins_min, nbins_max)
	end

	raw = span / h

	# If raw bin count is insane, just clamp **before** converting to Int
	if !isfinite(raw) || raw <= nbins_min
		return nbins_min
	elseif raw >= nbins_max
		return nbins_max
	else
		return ceil(Int, raw)
	end
end



# Build a piecewise-constant PDF from samples
function _pdf_from_hist(x::AbstractVector{<:Real}; nbins::Union{Int, Nothing} = nothing)
	n = length(x)
	n == 0 && error("Empty sample set.")

	nb = isnothing(nbins) ? _auto_nbins(x) : nbins

	h = fit(Histogram, float.(x); nbins = nb, closed = :left)
	edges = collect(h.edges[1])
	widths = diff(edges)
	dens = h.weights ./ (n .* widths)

	return LineParametersPDF(edges, dens)  # ctor re-normalizes area
end

# Density at x0 
@inline function (hp::LineParametersPDF)(x0::Real)
	i = searchsortedlast(hp.edges, float(x0))
	(i < 1 || i >= length(hp.edges)) && return 0.0
	return hp.dens[i]
end