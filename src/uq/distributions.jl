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

"""
Finds the index `i` of the bin `[edges[i], edges[i+1])` that `x` falls into.
Returns 0 if `x` is out of bounds.
Handles the right-most edge `x == edges[end]` correctly.
"""
function _binsearch(d::LineParametersPDF, x::Real)
	if x < d.edges[1] || x > d.edges[end]
		return 0 # Out of bounds
	end

	# Handle the maximum edge case, which searchsortedlast fucks up
	if x == d.edges[end]
		return length(d.dens) # Belongs to the last bin
	end

	# searchsortedlast finds the largest index i s.t. edges[i] <= x
	# This is exactly the bin index we need.
	i = searchsortedlast(d.edges, x)

	# This should be redundant given the initial check, but belt and suspenders.
	(i < 1 || i > length(d.dens)) && return 0

	return i
end

"""
Computes the stable integral of x^k over [a, b]
Returns: (b^(k+1) - a^(k+1)) / (k+1)
"""
function _stable_pow_integral(a::T, b::T, k::Int) where {T <: Real}
	n = k + 1
	h = b - a # width

	# If width is effectively zero, integral is zero
	if h == 0
		return zero(T)
	end

	# Use the stable factored form: (b-a)/n * sum(a^j * b^(n-1-j) for j=0..n-1)
	# n-1 = k
	s = zero(T)
	@inbounds for j in 0:k
		s += a^j * b^(k - j)
	end

	return s * h / n
end

"""
Computes raw moments m_1...m_K in a *single pass*.
Returns a Vector m where m[k] = E[X^k].
"""
function _raw_moments(d::LineParametersPDF{T}, K::Int) where {T <: Real}
	e = d.edges
	dens = d.dens

	# m[k] will hold the k-th raw moment
	m = zeros(T, K)

	@inbounds for i in 1:length(dens)
		# Skip bins with zero density.
		d_i = dens[i]
		d_i == 0 && continue

		a = e[i]
		b = e[i+1]

		# Calculate all moments 1..K for this bin
		for k in 1:K
			# This is the numerically stable integral of x^k from a to b,
			# which is (b^(k+1) - a^(k+1)) / (k+1).
			integral_term = _stable_pow_integral(a, b, k)

			# Add this bin's contribution to the k-th moment
			m[k] += d_i * integral_term
		end
	end
	return m
end

# ─────────────────────────────────────────────────────────────────────────────
# Distributions.jl API
# ─────────────────────────────────────────────────────────────────────────────

# --- Basic properties ---

Distributions.minimum(d::LineParametersPDF) = d.edges[1]
Distributions.maximum(d::LineParametersPDF) = d.edges[end]

function Distributions.insupport(d::LineParametersPDF, x::Real)
	# Is it in the bounds? This isn't rocket science.
	return d.edges[1] <= x <= d.edges[end]
end

# --- PDF / LOGPDF ---

function Distributions.pdf(d::LineParametersPDF{T}, x::Real) where {T}
	i = _binsearch(d, x)
	# If index is 0 (out of bounds), density is 0. Otherwise, look it up.
	return i == 0 ? zero(T) : d.dens[i]
end

function Distributions.logpdf(d::LineParametersPDF{T}, x::Real) where {T}
	p = Distributions.pdf(d, x)
	# Don't try to log(0). It's -Inf.
	return p > 0 ? log(p) : -Inf
end

# --- CDF (Cumulative Distribution Function) ---

function Distributions.cdf(d::LineParametersPDF{T}, x::Real) where {T}
	if x < Distributions.minimum(d)
		return zero(T)
	end
	if x >= Distributions.maximum(d)
		return one(T)
	end

	i_x = _binsearch(d, x) # The bin that x is currently in
	widths = diff(d.edges)

	# 1. Sum the area of all *full* bins before the current one
	area_full_bins = sum(
		(d.dens[j] * widths[j] for j in 1:(i_x-1));
		init = zero(T),
	)

	# 2. Add the partial area of the current bin
	area_partial_bin = d.dens[i_x] * (x - d.edges[i_x])

	return area_full_bins + area_partial_bin
end

# --- Quantile (Inverse CDF) ---

"""
Pre-calculates cumulative probabilities for efficient sampling.
This is what `sampler` should actually be doing.
"""
struct LineParametersPDFSampler{T <: Real} <:
	   Distributions.Sampleable{Distributions.Univariate, Distributions.Continuous}
	d::LineParametersPDF{T}
	cum_probs::Vector{T} # Cumulative probability at the *end* of each bin
end

function Distributions.sampler(d::LineParametersPDF)
	widths = diff(d.edges)
	bin_probs = d.dens .* widths
	cum_probs = cumsum(bin_probs)

	# Ensure the last value is exactly 1.0 to avoid float rounding
	# errors when sampling u=1.0
	cum_probs[end] = 1.0

	return LineParametersPDFSampler(d, cum_probs)
end

function Distributions.quantile(s::LineParametersPDFSampler{T}, q::Real) where {T}
	# This is the actual inverse-CDF logic.
	d = s.d

	if q <= 0
		return Distributions.minimum(d)
	end
	if q >= 1
		return Distributions.maximum(d)
	end

	# Find the first bin `i` where the cumulative probability >= q
	i = findfirst(p -> p >= q, s.cum_probs)
	# This should never be nothing thanks to the q >= 1 check, but...
	if i === nothing
		return Distributions.maximum(d)
	end

	# Get probability accumulated *before* this bin
	q_prev = (i == 1) ? zero(T) : s.cum_probs[i-1]

	# How much more probability do we need *from this bin*?
	q_needed = q - q_prev

	# If density is zero, any width is fine, just return the start edge.
	# Avoids a 0/0 NaN.
	if d.dens[i] <= 0
		return d.edges[i]
	end

	# Calculate the partial width into this bin
	# width = probability / density
	width_needed = q_needed / d.dens[i]

	return d.edges[i] + width_needed
end

# `quantile(d, q)` will be slow as it builds the sampler *every time*.
# This is the price you pay for a stateless distribution object.
function Distributions.quantile(d::LineParametersPDF, q::Real)
	return Distributions.quantile(Distributions.sampler(d), q)
end

# --- RAND (Random Sampling) ---

# Use the efficient sampler-based method
function Base.rand(rng::AbstractRNG, s::LineParametersPDFSampler)
	u = rand(rng, eltype(s.cum_probs)) # Get a random probability 0..1
	return Distributions.quantile(s, u) # Run it through the inverse CDF
end

# This one will be called if you just do `rand(d)`
function Base.rand(rng::AbstractRNG, d::LineParametersPDF)
	# This is inefficient as fuck. It builds the sampler on every. single. draw.
	# But it's what the Distributions.jl API expects as a fallback.
	# Use `rand(rng, sampler(d))` for batch sampling.
	s = Distributions.sampler(d)
	return Base.rand(rng, s)
end

# Get all moments up to k, then return the k-th
Distributions.moment(d::LineParametersPDF, k::Integer) = _raw_moments(d, k)[k]

function Distributions.mean(d::LineParametersPDF{T}) where {T}
	# E[X] = ∫ x * p(x) dx
	return Distributions.moment(d, 1)
end

function Distributions.var(d::LineParametersPDF{T}) where {T}
	# Var(X) = E[X^2] - (E[X])^2
	# E[X^2] = ∫ x^2 * p(x) dx
	# For bin i, integral is d.dens[i] * ∫(from e_i to e_{i+1}) x^2 dx
	m = _raw_moments(d, 2)
	m1 = m[1]
	m2 = m[2]

	v = m2 - m1^2
	# Handle floating point noise. Variance cannot be negative.
	return v < 0 ? zero(T) : v
end

Distributions.std(d::LineParametersPDF) = sqrt(Distributions.var(d))

function Distributions.skewness(d::LineParametersPDF{T}) where {T}
	m = _raw_moments(d, 3)
	m1, m2, m3 = m[1], m[2], m[3]

	μ = m1
	μ2 = m2 - μ^2 # Variance

	# Degenerate case: variance is zero. Return NaN.
	if μ2 <= eps(T) # Use machine epsilon for float comparison
		return T(NaN)
	end

	μ3 = m3 - 3*μ*m2 + 2*μ^3
	return μ3 / μ2^(3/2)
end

function Distributions.kurtosis(d::LineParametersPDF{T}) where {T}
	# Pass `true` for excess kurtosis (subtracts 3)
	return Distributions.kurtosis(d, true)
end

function Distributions.kurtosis(d::LineParametersPDF{T}, excess::Bool) where {T}
	m = _raw_moments(d, 4)
	m1, m2, m3, m4 = m[1], m[2], m[3], m[4]

	μ = m1
	μ2 = m2 - μ^2 # Variance

	# Degenerate case: variance is zero. Return NaN.
	if μ2 <= eps(T)
		return T(NaN)
	end

	μ4 = m4 - 4*μ*m3 + 6*μ^2*m2 - 3*μ^4

	kurt = μ4 / μ2^2
	return excess ? (kurt - 3) : kurt
end

function Distributions.mode(d::LineParametersPDF)
	# Returns *a* mode. The distribution is multi-modal
	# if the max density spans multiple (or disjoint) bins.
	# We'll just return the midpoint of the *first* bin with max density.

	max_dens, i = findmax(d.dens)
	return (d.edges[i] + d.edges[i+1]) / 2
end

function Distributions.modes(d::LineParametersPDF{T}) where {T}
	# The "modes" are technically *intervals*, not points.
	# This is a pain in the ass.
	# We'll just return the midpoints of all bins with max density.

	max_dens = maximum(d.dens)
	# Find all bins that are numerically close to the max
	indices = findall(p -> p ≈ max_dens, d.dens)

	return [(d.edges[i] + d.edges[i+1]) / 2 for i in indices]
end
