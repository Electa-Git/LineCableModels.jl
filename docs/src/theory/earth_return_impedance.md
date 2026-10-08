# Earth return impedance

```@meta
CurrentModule = LineCableModels.Engine.EarthImpedance
```

Earth-return formulas supply the self and mutual series terms associated
with the media outside each conductor or cable boundary. Their assumptions
specify the conductor arrangement (overhead, buried, or mixed), the earth
structure (homogeneous or layered), and the treatment of displacement current
and longitudinal propagation. The [matrix formulation](matrix_formulation.md) determines
where each contribution enters.

- [Source coefficients of the default formula in two half-spaces](@ref LineCableModels.Engine.EarthAdmittance.source_coefficients(::Union{LineCableModels.Engine.EarthImpedance.Formula{:unified}, LineCableModels.Engine.EarthAdmittance.Formula{:unified}}, ::Union{Val{:self}, Val{:mutual}}, ::Union{Val{1}, Val{2}}, ::Union{Val{1}, Val{2}}, ::Any, ::Any))
- [Carson homogeneous-earth overhead correction integral](external-impedance/1926/homogeneous-earth-overhead-integral/Carson1926.md)
- [Gary approximation for overhead wires using complex depth](@ref earth_impedance(::Formula{:gary1976}, ::Val{:self}, ::Val{1}, ::Val{1}, ::Any, ::Any, ::Any))
- [Lucca impedance for a mixed pair in homogeneous earth](@ref earth_impedance(::Formula{:lucca1994}, ::Val{:mutual}, ::Val{1}, ::Val{2}, ::Any, ::Any, ::Any))
- [Pollaczek generalized induction coefficients for overhead, buried, and mixed conductors](external-impedance/1926/homogeneous-earth-generalized-induction-green-function/Pollaczek1926.md)
- [Saad homogeneous-earth underground closed form](@ref earth_impedance(::Formula{:saad1996}, ::Val{:self}, ::Val{2}, ::Val{2}, ::Any, ::Any, ::Any))
- [Wise high-frequency overhead displacement-current integral](external-impedance/1934/homogeneous-earth-overhead-displacement-current-integral/Wise1934.md)
- [Wedepohl-Wilcox underground low-order impedance](@ref earth_impedance(::Formula{:wedepohl1973}, ::Val{:self}, ::Val{2}, ::Val{2}, ::Any, ::Any, ::Any))
- [Xue underground earth-return impedance](external-impedance/2018/complete-field-and-quasi-tem-underground/Xue2018.md)
- [Ametani mixed-pair exponential-image approximation](external-impedance/2009/homogeneous-earth-mixed-exponential-image/Ametani2009.md)

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

The default coefficients require the [complete-current matrix calculation](@ref LineCableModels.Engine.earth!(::Union{LineCableModels.Engine.EarthImpedance.Formula{:unified}, LineCableModels.Engine.EarthAdmittance.Formula{:unified}}, ::Any, ::Any)) before selecting physical exterior entries.

[Back to Contents](contents.md)
