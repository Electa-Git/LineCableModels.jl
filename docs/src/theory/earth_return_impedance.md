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

The default coefficients require the [complete-current matrix calculation](@ref LineCableModels.Engine.earth!(::Union{LineCableModels.Engine.EarthImpedance.Formula{:unified}, LineCableModels.Engine.EarthAdmittance.Formula{:unified}}, ::Any, ::Any)) before selecting physical exterior entries.

[Back to Contents](contents.md)
