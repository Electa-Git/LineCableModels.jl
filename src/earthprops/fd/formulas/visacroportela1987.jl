assumptions(::Val{:VisacroPortela1987}) = (normalization_frequency = 100.0,)

description(::Formula{:VisacroPortela1987}) =
    "Visacro–Portela 1987 (empirical soil dispersion)"

#=
Evaluate the Visacro–Portela empirical soil relation.

The reference conductivity is in S/m and the reference frequency is in Hz.
Returned resistivity uses Ω·m and relative permittivity is dimensionless.

# Reference

S. Visacro and C. M. Portela, *Soil Permittivity and Conductivity Behavior on
Frequency Range of Transient Phenomena in Electric Power Systems*, International
Symposium on High Voltage Engineering, 1987.
=#
function (::Functor{:VisacroPortela1987})(
        material::EarthMaterial{T},
        frequency::T,
        values::NamedTuple
) where {T <: Real}
    conductivity_reference = inv(material.rho)
    relative_permittivity = convert(T, 2.34e6) *
                            conductivity_reference^convert(T, 0.535) *
                            frequency^convert(T, -0.597)
    conductivity = conductivity_reference *
                   (frequency / convert(T, values.normalization_frequency))^convert(T, 0.072)
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

# Return the stable discovery identifier.
:VisacroPortela1987
