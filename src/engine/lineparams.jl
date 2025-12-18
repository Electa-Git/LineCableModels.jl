abstract type LineParamsDomain end
struct PhaseDomain <: LineParamsDomain end
struct ModalDomain <: LineParamsDomain end

struct SeriesImpedance{T} <: AbstractArray{T, 3}
	values::Array{T, 3}   # n×n×nfreq, units: Ω/m
end

struct ShuntAdmittance{T} <: AbstractArray{T, 3}
	values::Array{T, 3}   # n×n×nfreq, units: S/m
end

"""
$(TYPEDEF)

Represents the frequency-dependent line parameters (series impedance and shunt admittance matrices) for a cable or line system.

$(TYPEDFIELDS)
"""
struct LineParameters{T <: COMPLEXSCALAR, U <: REALSCALAR, D <: LineParamsDomain}
	"Series impedance matrices \\[Ω/m\\]."
	Z::SeriesImpedance{T}
	"Shunt admittance matrices \\[S/m\\]."
	Y::ShuntAdmittance{T}
	"Frequencies \\[Hz\\]."
	f::Vector{U}

	@doc """
	$(TYPEDSIGNATURES)

	Constructs a [`LineParameters`](@ref) instance.

	# Arguments

	- `Z`: Series impedance matrices \\[Ω/m\\].
	- `Y`: Shunt admittance matrices \\[S/m\\].
	- `f`: Frequencies \\[Hz\\].

	# Returns

	- A [`LineParameters`](@ref) object with prelocated impedance and admittance matrices for a given frequency range.

	# Examples

	```julia
	params = $(FUNCTIONNAME)(Z, Y, f)
	```
	"""
	function LineParameters(
		::Type{D},
		Z::SeriesImpedance{T},
		Y::ShuntAdmittance{T},
		f::AbstractVector{U},
	) where {D <: LineParamsDomain, T <: COMPLEXSCALAR, U <: REALSCALAR}
		size(Z, 1) == size(Z, 2) || throw(DimensionMismatch("Z must be square"))
		size(Y, 1) == size(Y, 2) || throw(DimensionMismatch("Y must be square"))
		size(Z, 3) == size(Y, 3) == length(f) ||
			throw(DimensionMismatch("Z and Y must have same dimensions (n×n×nfreq)"))
		new{T, U, D}(Z, Y, Vector{U}(f))
	end

	# Backward-compatible constructor: defaults to PhaseDomain
	LineParameters(
		Z::SeriesImpedance{T},
		Y::ShuntAdmittance{T},
		f::AbstractVector{U},
	) where {T <: COMPLEXSCALAR, U <: REALSCALAR} =
		LineParameters(PhaseDomain, Z, Y, f)
end

SeriesImpedance(A::AbstractArray{T, 3}) where {T} = SeriesImpedance{T}(Array(A))
ShuntAdmittance(A::AbstractArray{T, 3}) where {T} = ShuntAdmittance{T}(Array(A))

# --- Outer convenience constructors -------------------------------------------

"""
$(TYPEDSIGNATURES)

Construct from 3D arrays and frequency vector. Arrays are wrapped
into `SeriesImpedance` and `ShuntAdmittance` automatically.
"""
LineParameters(
	::Type{D},
	Z::AbstractArray{Tc, 3},
	Y::AbstractArray{Tc, 3},
	f::AbstractVector{U},
) where {D <: LineParamsDomain, Tc <: COMPLEXSCALAR, U <: REALSCALAR} =
	LineParameters(D, SeriesImpedance(Z), ShuntAdmittance(Y), f)


# Backward-compatible constructor: defaults to PhaseDomain
LineParameters(
	Z::AbstractArray{Tc, 3},
	Y::AbstractArray{Tc, 3},
	f::AbstractVector{U},
) where {Tc <: COMPLEXSCALAR, U <: REALSCALAR} =
	LineParameters(PhaseDomain, Z, Y, f)


# """
# $(TYPEDSIGNATURES)

# Backward-compatible constructor without frequencies. A dummy equally-spaced
# `Vector{BASE_FLOAT}` is used with length `size(Z,3)`.
# """
# function LineParameters(
# 	Z::AbstractArray{Tc, 3},
# 	Y::AbstractArray{Tc, 3},
# ) where {Tc <: COMPLEXSCALAR}
# 	nfreq = size(Z, 3)
# 	(size(Y, 3) == nfreq) || throw(DimensionMismatch("Z and Y must have same nfreq"))
# 	# Provide a placeholder frequency vector to preserve legacy call sites
# 	f = collect(BASE_FLOAT.(1:nfreq))
# 	return LineParameters(SeriesImpedance(Z), ShuntAdmittance(Y), f)
# end

# --- Tiny domain extractors ---------------------------------------------------
@inline domain(::Type{<:LineParameters{T, U, D}}) where {T, U, D <: LineParamsDomain} = D
@inline domain(lp::LineParameters) = domain(typeof(lp))
