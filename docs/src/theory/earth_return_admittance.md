# Earth return admittance

```@meta
CurrentModule = LineCableModels.Engine.EarthAdmittance
```

These formulations describe the external electric field through potential
coefficients or source-defined scalar admittances. The complete potential
matrix, including the local insulation contribution, is assembled before
conversion to shunt admittance; scalar reciprocals do not define a
multiconductor admittance matrix. Geometry, voltage reference, and
propagation assumptions remain those of each source.

- [Default two-half-space source-potential coefficients](@ref source_potential_coefficient(::Union{Formula{:unified}, Val{:unified}}, ::Union{Val{:self}, Val{:mutual}}, ::Val, ::Val, ::Any, ::Any, ::Any))
- [Ideal-earth electrostatic image potential coefficients](@ref earth_potential_coefficient(::Formula{:ideal}, ::Val{:self}, ::Val{1}, ::Val{1}, ::Any, ::Any, ::Any))
- [Pollaczek underground earth-return admittance](external-admittance/1926/homogeneous-earth-generalized-induction-green-function/Pollaczek1926.md)
- [Wise homogeneous-earth overhead potential coefficient](external-admittance/1948/homogeneous-earth-overhead-potential-coefficient/Wise1948.md)
- [Xue underground earth-return admittance](external-admittance/2018/complete-field-and-quasi-tem-underground/Xue2018.md)

The default coefficients require the [complete-current matrix calculation](@ref LineCableModels.Engine.earth!(::Union{LineCableModels.Engine.EarthImpedance.Formula{:unified}, LineCableModels.Engine.EarthAdmittance.Formula{:unified}}, ::Any, ::Any)) before selecting physical exterior entries.

[Back to Contents](contents.md)
