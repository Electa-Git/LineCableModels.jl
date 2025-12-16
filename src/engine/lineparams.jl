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
struct LineParameters{T <: COMPLEXSCALAR, U <: REALSCALAR}
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
		Z::SeriesImpedance{T},
		Y::ShuntAdmittance{T},
		f::AbstractVector{U},
	) where {T <: COMPLEXSCALAR, U <: REALSCALAR}
		size(Z, 1) == size(Z, 2) || throw(DimensionMismatch("Z must be square"))
		size(Y, 1) == size(Y, 2) || throw(DimensionMismatch("Y must be square"))
		size(Z, 3) == size(Y, 3) == length(f) ||
			throw(DimensionMismatch("Z and Y must have same dimensions (n×n×nfreq)"))
		new{T, U}(Z, Y, Vector{U}(f))
	end
end

SeriesImpedance(A::AbstractArray{T, 3}) where {T} = SeriesImpedance{T}(Array(A))
ShuntAdmittance(A::AbstractArray{T, 3}) where {T} = ShuntAdmittance{T}(Array(A))

# --- Outer convenience constructors -------------------------------------------

"""
$(TYPEDSIGNATURES)

Construct from 3D arrays and frequency vector. Arrays are wrapped
into `SeriesImpedance` and `ShuntAdmittance` automatically.
"""
function LineParameters(
	Z::AbstractArray{Tc, 3},
	Y::AbstractArray{Tc, 3},
	f::AbstractVector{U},
) where {Tc <: COMPLEXSCALAR, U <: REALSCALAR}
	return LineParameters(SeriesImpedance(Z), ShuntAdmittance(Y), f)
end

"""
$(TYPEDSIGNATURES)

Backward-compatible constructor without frequencies. A dummy equally-spaced
`Vector{BASE_FLOAT}` is used with length `size(Z,3)`.
"""
function LineParameters(
	Z::AbstractArray{Tc, 3},
	Y::AbstractArray{Tc, 3},
) where {Tc <: COMPLEXSCALAR}
	nfreq = size(Z, 3)
	(size(Y, 3) == nfreq) || throw(DimensionMismatch("Z and Y must have same nfreq"))
	# Provide a placeholder frequency vector to preserve legacy call sites
	f = collect(BASE_FLOAT.(1:nfreq))
	return LineParameters(SeriesImpedance(Z), ShuntAdmittance(Y), f)
end
