"""
	LineCableModels.Utils

The [`Utils`](@ref) module provides utility functions for the  [`LineCableModels.jl`](index.md) package. This module includes functions for handling measurements, numerical comparisons, and other common tasks.

# Overview

- Provides general constants used throughout the package.
- Includes utility functions for numerical comparisons and handling measurements.
- Contains functions to compute uncertainties and bounds for measurements.

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module Utils

# Export public API
export resolve_T, coerce_to_T, is_headless, is_in_testset, display_path
export set_verbosity!

export to_nominal,
	to_certain,
	percent_to_uncertain,
	bias_to_uncertain,
	to_upper,
	to_lower,
	percent_error

export _to_σ, _bessel_diff, symtrans!, line_transpose!

# Module-specific dependencies
using ..Commons
using ..UncertainBessels: besselk
using Measurements: Measurement, value, uncertainty, measurement, ±, Measurements, result
using Statistics
using Plots
using LinearAlgebra

"""
$(TYPEDSIGNATURES)

Returns the nominal (deterministic) value of inputs that may contain
`Measurements.Measurement` numbers, recursively handling Complex and arrays.

# Arguments

- `x`: Input value which can be a `Measurement` type or any other type.

# Returns

- Measurement → its `value`
- Complex → `complex(to_nominal(real(z)), to_nominal(imag(z)))`
- AbstractArray → broadcasts `to_nominal` elementwise
- Anything else → returned unchanged

# Examples

```julia
using Measurements

$(FUNCTIONNAME)(1.0)  # Output: 1.0
$(FUNCTIONNAME)(5.2 ± 0.3)  # Output: 5.2
```
"""
to_nominal(x::Measurement) = value(x)
to_nominal(z::Complex) = complex(to_nominal(real(z)), to_nominal(imag(z)))
to_nominal(A::AbstractArray) = to_nominal.(A)
to_nominal(x) = x

"""
$(TYPEDSIGNATURES)

Converts a measurement to a value with zero uncertainty, retaining the numeric type `Measurement`.

# Arguments

- `value`: Input value that may be a `Measurement` type or another type.

# Returns

- If input is a `Measurement`, returns the same value with zero uncertainty; otherwise returns the original value unchanged.

# Examples

```julia
x = 5.0 ± 0.1
result = $(FUNCTIONNAME)(x)  # Output: 5.0 ± 0.0

y = 10.0
result = $(FUNCTIONNAME)(y)  # Output: 10.0
```
"""
function to_certain(value)
	return value isa Measurement ? (Measurements.value(value) ± 0.0) : value
end

"""
$(TYPEDSIGNATURES)

Converts a value to a measurement with uncertainty based on percentage.

# Arguments

- `val`: The nominal value.
- `perc`: The percentage uncertainty (0 to 100).

# Returns

- A `Measurement` type with the given value and calculated uncertainty.

# Examples

```julia
using Measurements

$(FUNCTIONNAME)(100.0, 5)  # Output: 100.0 ± 5.0
$(FUNCTIONNAME)(10.0, 10)  # Output: 10.0 ± 1.0
```
"""
function percent_to_uncertain(val, perc) #perc from 0 to 100
	measurement(val, (perc * val) / 100)
end

"""
$(TYPEDSIGNATURES)

Computes the uncertainty of a measurement by incorporating systematic bias.

# Arguments

- `nominal`: The deterministic nominal value (Float64).
- `measurements`: A vector of `Measurement` values from the `Measurements.jl` package.

# Returns

- A new `Measurement` object representing the mean measurement value with an uncertainty that accounts for both statistical variation and systematic bias.

# Notes

- Computes the mean value and its associated uncertainty from the given `measurements`.
- Determines the **bias** as the absolute difference between the deterministic `nominal` value and the mean measurement.
- The final uncertainty is the sum of the standard uncertainty (`sigma_mean`) and the systematic bias.

# Examples
```julia
using Measurements

nominal = 10.0
measurements = [10.2 ± 0.1, 9.8 ± 0.2, 10.1 ± 0.15]
result = $(FUNCTIONNAME)(nominal, measurements)
println(result)  # Output: Measurement with adjusted uncertainty
```
"""
function bias_to_uncertain(nominal::Float64, measurements::Vector{<:Measurement})
	# Compute the mean value and uncertainty from the measurements
	mean_measurement = mean(measurements)
	mean_value = Measurements.value(mean_measurement)  # Central value
	sigma_mean = Measurements.uncertainty(mean_measurement)  # Uncertainty of the mean
	# Compute the bias (deterministic nominal value minus mean measurement)
	bias = abs(nominal - mean_value)
	return mean_value ± (sigma_mean + bias)
end

"""
$(TYPEDSIGNATURES)

Computes the upper bound of a measurement value.

# Arguments

- `m`: A numerical value, expected to be of type `Measurement` from the `Measurements.jl` package.

# Returns

- The upper bound of `m`, computed as `value(m) + uncertainty(m)` if `m` is a `Measurement`.
- `NaN` if `m` is not a `Measurement`.

# Examples

```julia
using Measurements

m = 10.0 ± 2.0
upper = $(FUNCTIONNAME)(m)  # Output: 12.0

not_a_measurement = 5.0
upper_invalid = $(FUNCTIONNAME)(not_a_measurement)  # Output: NaN
```
"""
function to_upper(m::Number)
	if m isa Measurement
		return Measurements.value(m) + Measurements.uncertainty(m)
	else
		return NaN
	end
end

"""
$(TYPEDSIGNATURES)

Computes the lower bound of a measurement value.

# Arguments

- `m`: A numerical value, expected to be of type `Measurement` from the `Measurements.jl` package.

# Returns

- The lower bound, computed as `value(m) - uncertainty(m)` if `m` is a `Measurement`.
- `NaN` if `m` is not a `Measurement`.

# Examples

```julia
using Measurements

m = 10.0 ± 2.0
lower = $(FUNCTIONNAME)(m)  # Output: 8.0

not_a_measurement = 5.0
lower_invalid = $(FUNCTIONNAME)(not_a_measurement)  # Output: NaN
```
"""
function to_lower(m::Number)
	if m isa Measurement
		return Measurements.value(m) - Measurements.uncertainty(m)
	else
		return NaN
	end
end

"""
$(TYPEDSIGNATURES)

Computes the percentage uncertainty of a measurement.

# Arguments

- `m`: A numerical value, expected to be of type `Measurement` from the `Measurements.jl` package.

# Returns

- The percentage uncertainty, computed as `100 * uncertainty(m) / value(m)`, if `m` is a `Measurement`.
- `NaN` if `m` is not a `Measurement`.

# Examples

```julia
using Measurements

m = 10.0 ± 2.0
percent_err = $(FUNCTIONNAME)(m)  # Output: 20.0

not_a_measurement = 5.0
percent_err_invalid = $(FUNCTIONNAME)(not_a_measurement)  # Output: NaN
```
"""
function percent_error(m::Number)
	if m isa Measurement
		return 100 * Measurements.uncertainty(m) / Measurements.value(m)
	else
		return NaN
	end
end

@inline _nudge_float(x::AbstractFloat) = isfinite(x) && x == trunc(x) ? nextfloat(x) : x #redundant and I dont care

_coerce_args_to_T(args...) =
	any(x -> x isa Measurement, args) ? Measurement{BASE_FLOAT} : BASE_FLOAT

# Promote scalar to T if T is Measurement; otherwise take nominal if x is Measurement.
function _coerce_scalar_to_T(x, ::Type{T}) where {T}
	if T <: Measurement
		return x isa Measurement ? x : (zero(T) + x)
	else
		return x isa Measurement ? T(value(x)) : convert(T, x)
	end
end

# Arrays: promote/demote elementwise, preserving shape. Arrays NEVER decide T.
function _coerce_array_to_T(A::AbstractArray, ::Type{T}) where {T}
	if T <: Measurement
		return (eltype(A) === T) ? A : (A .+ zero(T))             # Real → Measurement(σ=0)
	elseif eltype(A) <: Measurement
		B = value.(A)                                             # Measurement → Real (nominal)
		return (eltype(B) === T) ? B : convert.(T, B)
	else
		return (eltype(A) === T) ? A : convert.(T, A)
	end
end

"""
$(TYPEDSIGNATURES)

Determines if the current execution environment is headless (without display capability).

# Returns

- `true` if running in a continuous integration environment or without display access.
- `false` otherwise when a display is available.

# Examples

```julia
if $(FUNCTIONNAME)()
	# Use non-graphical backend
	gr()
else
	# Use interactive backend
	plotlyjs()
end
```
"""
function is_headless()::Bool
	# 1. Check for common CI environment variables
	if get(ENV, "CI", "false") == "true"
		return true
	end

	# 2. Check if a display is available (primarily for Linux)
	if !haskey(ENV, "DISPLAY") && Sys.islinux()
		return true
	end

	# 3. Check for GR backend's specific headless setting
	if get(ENV, "GKSwstype", "") in ("100", "nul", "nil")
		return true
	end

	return false
end

function display_path(file_name)
	return is_headless() ? basename(file_name) : relpath(file_name)
end

"""
$(TYPEDSIGNATURES)

Checks if the code is running inside a `@testset` by checking if `Test` is loaded
in the current session and then calling `get_testset_depth()`.
"""
function is_in_testset()
	# Start with the current module
	current_module = @__MODULE__

	# Walk up the module tree (e.g., from the sandbox to Main)
	while true
		if isdefined(current_module, :Test) &&
		   isdefined(current_module.Test, :get_testset_depth)
			# Found the Test module, check the test set depth
			return current_module.Test.get_testset_depth() > 0
		end

		# Move to the parent module
		parent = parentmodule(current_module)
		if parent === current_module # Reached the top (Main)
			break
		end
		current_module = parent
	end

	return false
end



"""
Apply `f` to every square block of `M` defined by `map`, in-place.

- `M`: n×n or n×n×nf (any eltype).
- `map`: length-n; equal ids => same block (non-contiguous ok).
- `f`: function like `f(B::AbstractMatrix, args...) -> k×k Matrix`.
- `args...`: extra positional args passed to `f`.
- `slice_positions`: positions in `args` that should be indexed as `args[i][idx]`
   per block (useful for `phase_map`).

Returns `M`.
"""
function block_transform!(M,
	map::AbstractVector{<:Integer},
	f::F,
	args...; slice_positions = Int[]) where {F}
	n = size(M, 1)
	(size(M, 2) == n && length(map) == n) || throw(ArgumentError("shape mismatch"))
	groups = unique(map)                     # preserve first-seen order
	blocks = [findall(==(g), map) for g in groups]

	# helper to build per-block args (slice selected ones)
	make_args(idx) =
		ntuple(i -> (i in slice_positions ? args[i][idx] : args[i]), length(args))

	if ndims(M) == 2
		for idx in blocks
			Bv = @view M[idx, idx]
			R = f(Matrix(Bv), make_args(idx)...)  #  f decides what to do
			size(R) == size(Bv) || throw(ArgumentError("f must return $(size(Bv))"))
			@inbounds Bv .= R
		end
	elseif ndims(M) == 3
		_, _, nf = size(M)
		for k in 1:nf
			for idx in blocks
				Bv = @view M[idx, idx, k]
				R = f(Matrix(Bv), make_args(idx)...)
				size(R) == size(Bv) || throw(ArgumentError("f must return $(size(Bv))"))
				@inbounds Bv .= R
			end
		end
	else
		throw(ArgumentError("M must be 2D or 3D"))
	end
	return M
end

# Non-mutating
block_transform(M, cmap, f, args...; slice_positions = Int[]) =
	block_transform!(copy(M), cmap, f, args...; slice_positions = slice_positions)


# Reciprocity  symmetrization — in place
symtrans!(A) = (A .= 0.5 .* (A .+ transpose(A)); A)

# Reciprocity symmetrization (power lines want transpose, not adjoint)
symtrans(A) = (A .+ transpose(A)) / 2


# Circulant projection (N×N), least-squares fit: C[i,j] = c[(j-i) mod N]
function line_transpose!(A::AbstractMatrix)
	n = size(A, 1);
	n == size(A, 2) || throw(ArgumentError("square"))
	c = similar(diag(A))  # length n

	# Average wrap-diagonals (use mod to avoid negatives)
	@inbounds for k in 0:(n-1)
		s = zero(eltype(A))
		for i in 1:n
			j = 1 + mod(i-1 + k, n)
			s += A[i, j]
		end
		c[k+1] = s / n
	end
	# Write back circulant matrix
	@inbounds for i in 1:n, j in 1:n
		A[i, j] = c[mod1(j - i + 1, n)] #c[1+mod(j-i, n)]
	end
	return A
end


function isdiag_approx(A; rtol = 1e-8, atol = 1e-8)
	isapprox(A, Diagonal(diag(A)); rtol = rtol, atol = atol)
end

function offdiag_ratio(A)
	n = size(A, 1)
	n == size(A, 2) || throw(ArgumentError("square"))
	T = real(float(eltype(A)))
	dmax = zero(T)
	odmax = zero(T)
	@inbounds for j in 1:n
		dj = abs(A[j, j])
		dmax = dj > dmax ? dj : dmax
		for i in 1:n
			i == j && continue
			v = abs(A[i, j])
			odmax = v > odmax ? v : odmax
		end
	end
	return odmax / max(dmax, eps(T))
end

isdiag_rel(A; τ = 1e-4) = offdiag_ratio(A) ≤ τ

function issymmetric_approx(A; rtol = 1e-8, atol = 1e-8)
	size(A, 1) == size(A, 2) || return false
	return isapprox(A, transpose(A); rtol = rtol, atol = atol)
end


@inline _to_σ(ρ) = isinf(ρ) ? zero(ρ) : (iszero(ρ) ? inv(zero(ρ)) : inv(ρ))

@inline function _bessel_diff(γs, d::T, D::T) where {T}
	zmax = max(abs(γs)*d, abs(γs)*D)
	return isapprox(to_nominal(zmax), 0.0, atol = TOL) ? log(D/d) :
		   (besselk(0, γs*d) - besselk(0, γs*D))
end

include("logging.jl")
include("typecoercion.jl")
include("macros.jl")

end # module Utils
