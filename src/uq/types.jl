# ─────────────────────────────────────────────────────────────────────────────
# Piecewise-constant PDF as a ContinuousUnivariateDistribution
# ─────────────────────────────────────────────────────────────────────────────
struct LineParametersPDF{T <: Real} <: ContinuousUnivariateDistribution
	edges::Vector{T}   # length B+1, sorted ascending
	dens::Vector{T}    # length B, area ≈ 1
end

# Your constructor (modified slightly to remove redundant collect)
function LineParametersPDF(
	edges::AbstractVector{T},
	dens::AbstractVector{T},
) where {T <: Real}
	length(edges) == length(dens) + 1 ||
		throw(ArgumentError("edges must have length length(dens)+1"))

	e = Vector(edges) # Just convert once
	d = Vector(dens)

	# Ensure increasing edges
	issorted(e) || throw(ArgumentError("edges must be sorted ascending"))

	# Normalize area to 1 (robust against floating crap)
	widths = diff(e)
	area = dot(d, widths) # Cleaner than sum(d .* widths)
	area <= zero(T) && throw(ArgumentError("non-positive total area in LineParametersPDF"))

	# Only normalize if it's not already 1 (avoids float division noise)
	if !(area ≈ 1.0)
		d ./= area
	end

	return LineParametersPDF{T}(e, d)
end

struct LineParametersMCSummary{U <: Real}
	"Frequencies \\[Hz\\]."
	f::Vector{U}

	"Statistics tensors for R, L, C, G \\[per element and frequency\\]."
	stats::NamedTuple{
		(:R, :L, :C, :G),
		Tuple{
			Array{NamedTuple, 3},  # R[i,j,k]
			Array{NamedTuple, 3},  # L[i,j,k]
			Array{NamedTuple, 3},  # C[i,j,k]
			Array{NamedTuple, 3},  # G[i,j,k]
		},
	}

	"Empirical PDFs per R, L, C, G entry or `nothing` if not requested."
	pdf::Union{
		Nothing,
		NamedTuple{
			(:R, :L, :C, :G),
			Tuple{
				Array{LineParametersPDF{U}, 3},
				Array{LineParametersPDF{U}, 3},
				Array{LineParametersPDF{U}, 3},
				Array{LineParametersPDF{U}, 3},
			},
		},
	}

	"Optional Monte Carlo samples of R/L/C/G."
	samples::Union{
		Nothing,
		NamedTuple{
			(:R, :L, :C, :G),
			Tuple{
				Array{U, 4},  # R[i,j,k,trial]
				Array{U, 4},  # L
				Array{U, 4},  # C
				Array{U, 4},  # G
			},
		},
	}

	"Frequency-dependent LineParameters with Measurement-valued entries."
	measurements::LineParameters{Complex{Measurement{U}}, U}
end

struct CableDesignMCSummary{U <: Real}
	"Statistics for R, L, C (each is a NamedTuple from the mc stats kernel)."
	stats::NamedTuple{
		(:R, :L, :C),
		Tuple{NamedTuple, NamedTuple, NamedTuple},
	}

	"Empirical PDFs for R, L, C or `nothing` if not requested."
	pdf::Union{
		Nothing,
		NamedTuple{
			(:R, :L, :C),
			Tuple{LineParametersPDF{U}, LineParametersPDF{U}, LineParametersPDF{U}},
		},
	}

	"Optional raw samples of R, L, C (each a Vector of length ntrials)."
	samples::Union{
		Nothing,
		NamedTuple{
			(:R, :L, :C),
			Tuple{Vector{U}, Vector{U}, Vector{U}},
		},
	}

	"Measurements for R, L, C (mean ± std)."
	measurements::Vector{Measurement{U}}
end
