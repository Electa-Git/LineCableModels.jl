"""
	LineCableModels.UncertainBessels

Uncertainty-aware wrappers for Bessel functions.

[`UncertainBessels`](@ref) lifts selected functions from `SpecialFunctions` so they accept
`Measurement` and `Complex{Measurement}` inputs. The wrapper evaluates the
underlying function at the nominal complex argument and propagates uncertainty
via first-order finite differences using the four partial derivatives ``\\frac{\\partial \\mathrm{Re} \\, f}{\\partial x}, \\frac{\\partial \\mathrm{Re} \\, f}{\\partial y}, \\frac{\\partial \\mathrm{Im} \\, f}{\\partial x}, \\frac{\\partial \\mathrm{Im} \\, f}{\\partial y}`` with ``x = \\mathrm{Re}(z)`` and ``y = \\mathrm{Im}(z)``. No new Bessel algorithms are implemented: for plain numeric inputs, results and numerical behaviour are those of
`SpecialFunctions`.

Numerical scaling (as defined by `SpecialFunctions`) is supported for the
“x” variants (e.g. `besselix`, `besselkx`, `besseljx`, …) to improve stability
for large or complex arguments. In particular, the modified functions use
exponential factors to temper growth along ``\\mathrm{Re}(z)`` (e.g. ``I_\\nu`` and ``K_\\nu``);
other scaled variants follow conventions in `SpecialFunctions` and DLMF guidance
for complex arguments. See [NIST:DLMF](@cite) and [6897971](@cite).

# Overview

- Thin, uncertainty-aware wrappers around `SpecialFunctions` (`besselj`, `bessely`,
  `besseli`, `besselk`, `besselh`) and their scaled counterparts (`…x`).
- For `Complex{Measurement}` inputs, uncertainty is propagated using the 4-component
  gradient with respect to ``\\mathrm{Re}(z)`` and ``\\mathrm{Im}(z)``.
- For `Measurement` (real) inputs, a 1-D finite-difference derivative is used.
- No change in semantics for `Real`/`Complex` inputs: calls delegate to `SpecialFunctions`.

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)

# Usage

```julia
# do not import SpecialFunctions directly
using LineCableModels.UncertainBessels 
z = complex(1.0, 1.0 ± 0.5)
J0_cpl = besselj(0, z) 			# Complex{Measurement}
J0_nom = besselj(0, value(z)) 	# nominal comparison
I1 = besselix(1, z) 			# scaled I1 with uncertainty
```

# Numerical notes

- Scaled modified Bessels remove large exponential factors along ``\\mathrm{Re}(z)`` (e.g., ``I_\\nu`` and ``K_\\nu`` are scaled by opposite signs of ``|\\mathrm{Re}(z)|``), improving conditioning. Scaled forms for the other families follow the definitions in `SpecialFunctions` and DLMF.
- Uncertainty propagation is first order (linearization at the nominal point).
  Large uncertainties or strong nonlinearity may reduce accuracy.

# See also

- [`LineCableModels.Engine.InternalImpedance`](@ref)
- [`LineCableModels.Engine.EarthImpedance`](@ref)
"""
module UncertainBessels

# Module-specific dependencies
using ..Commons
using Calculus: Calculus
using SpecialFunctions: SpecialFunctions
using Measurements: Measurements, Measurement

export besselix, besselkx, besseljx, besselyx, besselhx
export besseli, besselk, besselj, bessely, besselh

# Complex argument with measurement parts
@inline function _lift_complex_measurement(f, ν, ẑ::Complex{<:Measurement})
	return Measurements.result(
		f(ν, Measurements.value(ẑ)),
		vcat(
			Calculus.gradient(
				x -> real(f(ν, complex(x[1], x[2]))),
				[reim(Measurements.value(ẑ))...],
			),
			Calculus.gradient(
				x -> imag(f(ν, complex(x[1], x[2]))),
				[reim(Measurements.value(ẑ))...],
			),
		),
		ẑ,
	)
end

# Real argument with measurement
@inline function _lift_real_measurement(f, ν, x::Measurements.Measurement)
	x0 = Measurements.value(x)
	y0 = f(ν, x0)
	dy = Calculus.derivative(t -> f(ν, t), x0)
	return Measurements.result(y0, (dy,), x)
end


# Complex inputs with uncertainty
@inline besselix(ν, z::Complex{<:Measurements.Measurement{T}}) where {T <: AbstractFloat} =
	_lift_complex_measurement(SpecialFunctions.besselix, ν, z)
@inline besselkx(ν, z::Complex{<:Measurements.Measurement{T}}) where {T <: AbstractFloat} =
	_lift_complex_measurement(SpecialFunctions.besselkx, ν, z)
@inline besseljx(ν, z::Complex{<:Measurements.Measurement{T}}) where {T <: AbstractFloat} =
	_lift_complex_measurement(SpecialFunctions.besseljx, ν, z)
@inline besselyx(ν, z::Complex{<:Measurements.Measurement{T}}) where {T <: AbstractFloat} =
	_lift_complex_measurement(SpecialFunctions.besselyx, ν, z)
@inline besselhx(ν, z::Complex{<:Measurements.Measurement{T}}) where {T <: AbstractFloat} =
	_lift_complex_measurement(SpecialFunctions.besselhx, ν, z)
@inline besselj(ν, z::Complex{<:Measurements.Measurement{T}}) where {T <: AbstractFloat} =
	_lift_complex_measurement(SpecialFunctions.besselj, ν, z)
@inline bessely(ν, z::Complex{<:Measurements.Measurement{T}}) where {T <: AbstractFloat} =
	_lift_complex_measurement(SpecialFunctions.bessely, ν, z)
@inline besseli(ν, z::Complex{<:Measurements.Measurement{T}}) where {T <: AbstractFloat} =
	_lift_complex_measurement(SpecialFunctions.besseli, ν, z)
@inline besselk(ν, z::Complex{<:Measurements.Measurement{T}}) where {T <: AbstractFloat} =
	_lift_complex_measurement(SpecialFunctions.besselk, ν, z)
@inline besselh(ν, z::Complex{<:Measurements.Measurement{T}}) where {T <: AbstractFloat} =
	_lift_complex_measurement(SpecialFunctions.besselh, ν, z)

# Real inputs with uncertainty
@inline besselix(ν, x::Measurements.Measurement{T}) where {T <: AbstractFloat} =
	_lift_real_measurement(SpecialFunctions.besselix, ν, x)
@inline besselkx(ν, x::Measurements.Measurement{T}) where {T <: AbstractFloat} =
	_lift_real_measurement(SpecialFunctions.besselkx, ν, x)
@inline besseljx(ν, x::Measurements.Measurement{T}) where {T <: AbstractFloat} =
	_lift_real_measurement(SpecialFunctions.besseljx, ν, x)
@inline besselyx(ν, x::Measurements.Measurement{T}) where {T <: AbstractFloat} =
	_lift_real_measurement(SpecialFunctions.besselyx, ν, x)
@inline besselhx(ν, x::Measurements.Measurement{T}) where {T <: AbstractFloat} =
	_lift_real_measurement(SpecialFunctions.besselhx, ν, x)
@inline besselj(ν, x::Measurements.Measurement{T}) where {T <: AbstractFloat} =
	_lift_real_measurement(SpecialFunctions.besselj, ν, x)
@inline bessely(ν, x::Measurements.Measurement{T}) where {T <: AbstractFloat} =
	_lift_real_measurement(SpecialFunctions.bessely, ν, x)
@inline besseli(ν, x::Measurements.Measurement{T}) where {T <: AbstractFloat} =
	_lift_real_measurement(SpecialFunctions.besseli, ν, x)
@inline besselk(ν, x::Measurements.Measurement{T}) where {T <: AbstractFloat} =
	_lift_real_measurement(SpecialFunctions.besselk, ν, x)
@inline besselh(ν, x::Measurements.Measurement{T}) where {T <: AbstractFloat} =
	_lift_real_measurement(SpecialFunctions.besselh, ν, x)

# Plain Float/Complex fallbacks 
@inline besselix(ν, z::T) where {T <: AbstractFloat} = SpecialFunctions.besselix(ν, z)
@inline besselkx(ν, z::T) where {T <: AbstractFloat} = SpecialFunctions.besselkx(ν, z)
@inline besseljx(ν, z::T) where {T <: AbstractFloat} = SpecialFunctions.besseljx(ν, z)
@inline besselyx(ν, z::T) where {T <: AbstractFloat} = SpecialFunctions.besselyx(ν, z)
@inline besselhx(ν, z::T) where {T <: AbstractFloat} = SpecialFunctions.besselhx(ν, z)
@inline besselj(ν, z::T) where {T <: AbstractFloat} = SpecialFunctions.besselj(ν, z)
@inline bessely(ν, z::T) where {T <: AbstractFloat} = SpecialFunctions.bessely(ν, z)
@inline besseli(ν, z::T) where {T <: AbstractFloat} = SpecialFunctions.besseli(ν, z)
@inline besselk(ν, z::T) where {T <: AbstractFloat} = SpecialFunctions.besselk(ν, z)
@inline besselh(ν, z::T) where {T <: AbstractFloat} = SpecialFunctions.besselh(ν, z)

@inline besselix(ν, z::Complex{T}) where {T <: AbstractFloat} =
	SpecialFunctions.besselix(ν, z)
@inline besselkx(ν, z::Complex{T}) where {T <: AbstractFloat} =
	SpecialFunctions.besselkx(ν, z)
@inline besseljx(ν, z::Complex{T}) where {T <: AbstractFloat} =
	SpecialFunctions.besseljx(ν, z)
@inline besselyx(ν, z::Complex{T}) where {T <: AbstractFloat} =
	SpecialFunctions.besselyx(ν, z)
@inline besselhx(ν, z::Complex{T}) where {T <: AbstractFloat} =
	SpecialFunctions.besselhx(ν, z)
@inline besselj(ν, z::Complex{T}) where {T <: AbstractFloat} =
	SpecialFunctions.besselj(ν, z)
@inline bessely(ν, z::Complex{T}) where {T <: AbstractFloat} =
	SpecialFunctions.bessely(ν, z)
@inline besseli(ν, z::Complex{T}) where {T <: AbstractFloat} =
	SpecialFunctions.besseli(ν, z)
@inline besselk(ν, z::Complex{T}) where {T <: AbstractFloat} =
	SpecialFunctions.besselk(ν, z)
@inline besselh(ν, z::Complex{T}) where {T <: AbstractFloat} =
	SpecialFunctions.besselh(ν, z)

end # module UncertainBessels
