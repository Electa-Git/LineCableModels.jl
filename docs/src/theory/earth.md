# Earth return

```@meta
CurrentModule = LineCableModels.Engine
```

## The physical model

The air above the interface `y = 0` and the earth below it are homogeneous half-spaces, each
with its resistivity ``\rho``, relative permittivity ``\varepsilon_r`` and relative
permeability ``\mu_r``. Each cable is a source of circumferential currents on its exterior
circle, which lies wholly in one medium. The formula closes the circumference currents of all
cables completely through the matrix `L` of its state, and the fields vary along the line as
``e^{-\Gamma x}``, with a prescribed longitudinal wavenumber ``\Gamma``, zero by default.

The formula gives two quantities per pair of cables, from their exterior surfaces: the
earth-return impedance `Zg` and the earth-return potential coefficients `Pg`. They are the
matrices ``\mathbf Z_{\mathrm{env}}`` and ``\mathbf P_{\mathrm{env}}`` of the
[matrix formulation](matrix.md), added to every pair of conductors of each pair of cables.

The fields of each medium are spectral integrals over the real line of the transverse
wavenumber. The formula evaluates them on a complex ray, which [`ray!`](@ref) chooses with the
breakpoints of the quadrature, and [`integrate`](@ref) evaluates each integral by adaptive
quadrature. The option `integral` of the formula selects the method and its controls,
`integral = (method = :quad, rtol, atol, maxevals)`.

## Validity of the default formula

The default formula, `:unified`, holds within two ranges.

- **Prescribed Γ.** A nonzero Γ lies between the propagation constants of the two media. It
  is in range at a frequency where Im γ_air ≤ Im Γ ≤ Im γ_earth and Re Γ ≤ Re γ_earth. Γ = 0,
  the default, always lies in range.
- **Thick receiver.** The formula represents each receiving conductor by the mean field on
  its exterior circumference. A conductor that is thick compared with the transverse
  wavelength of its medium receives a voltage that departs from the full-field value. The
  error grows as (κ_m r_p)², where κ_m is the outgoing root of γ_m² − Γ² in the medium of
  conductor p and r_p is its exterior radius, the jacket for an insulated cable. The range
  is |κ_m r_p| ≤ 0.1. Emission from a thick conductor is accurate.

The calibration below is indicative. It compares the default formula with a finite-element
reference for two bare copper conductors buried at 1 m depth and 2 m apart, in earth of
1000 Ω·m and relative permittivity 12, at 100 MHz with Γ = 0, where |κ_earth| ≈ 7.26 1/m.

| Radii (m) | κr | Scaled Y error | P error, thick receiver | P error, thin receiver |
| --- | ---: | ---: | ---: | ---: |
| 0.002 / 0.002 | 0.015 | 0.74 % | none | ≤ 1.0 % (all entries) |
| 0.01 / 0.01 | 0.073 | 1.08 % | none | not reported |
| 0.0425 / 0.002 | 0.31 / 0.015 | 1.55 % | 11.0 % (mutual) | 2.0 % (mutual) |
| 0.0425 / 0.0425 | 0.31 | 5.18 % | 11 % mutual, 4.5 % self | none |

The scaled Y error is max|ΔY| / max|diag Y_ref|. P = jω Y⁻¹, with one row per receiver. The
reference itself has an error floor of about 0.7 %, the first row. The calibration covers
buried conductors only. In air at Γ = 0, κ ≈ k₀, so a 20 mm radius reaches κr ≈ 0.13 at
300 MHz, a regime outside the calibration.

## The properties of the air and the earth

The law of the slot `earth_properties` gives the resistivity, the relative permittivity and the
relative permeability of the air and of the earth at each frequency. The one law of this
release, `:constant`, gives the static values of each layer at every frequency.

## Earth-return formulas

```@autodocs
Modules = [EarthReturn]
Order = [:module, :constant, :type, :function]
Public = true
Private = false
```

```@docs
defaults(::Type{EarthReturn.Formula{:unified}})
validate(options::FormulationOptions{EarthReturn.Formula{:unified}})
validate(problem::LineParametersProblem, formula::EarthReturn.Formula{:unified})
initialize_buffers(formula::EarthReturn.Formula{:unified}, ::LineParametersProblem, input::CoaxialInput{T}, buffers) where {T}
prepare!(formula::EarthReturn.Formula{:unified}, buffers, input::CoaxialInput, media, s, k)
impedance!(formula::EarthReturn.Formula{:unified}, Zg::AbstractMatrix, buffers, input, media, s, k)
potential!(formula::EarthReturn.Formula{:unified}, Pg::AbstractMatrix, buffers, input, media, s, k)
```

## Spectral integrals

```@docs
SpectralIntegral
integrate
ray!
defaults(::typeof(integrate), ::Val{:quad})
validate(section::NamedTuple, ::typeof(integrate), ::Val{:quad})
select(default::Pair{typeof(integrate), Symbol}, given::NamedTuple)
initialize_buffers(::Type{SpectralIntegral}, ::Val{:quad}, ::Type{T}, buffers) where {T}
```

## Earth properties

```@autodocs
Modules = [EarthProperties]
Order = [:module, :constant, :type, :function]
Public = true
Private = false
```

```@docs
properties(::EarthProperties.Formula{:constant}, layer::EarthProperties.EarthLayer{T}, f::T) where {T <: Real}
```

[Back to Contents](contents.md)
