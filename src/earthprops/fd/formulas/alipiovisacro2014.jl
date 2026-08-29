"Return the fitted assumptions of the Alipio–Visacro causal soil model."
assumptions(::Val{:AlipioVisacro2014}) = (
    exponent = 0.54,
    epsilon_infinity = 12.0,
    scale = 1.26,
    conductivity_exponent = -0.73
)

description(::Formula{:AlipioVisacro2014}) =
    "Alipio–Visacro 2014 (causal soil dispersion)"

#=
Evaluate the Alipio–Visacro causal soil relation.

The reference conductivity is converted from S/m to mS/m for the published
fit. The returned resistivity is converted back to Ω·m and relative
permittivity is dimensionless.

# Reference

R. Alipio and S. Visacro, *Modeling the Frequency Dependence of Electrical
Parameters of Soil*, IEEE Transactions on Electromagnetic Compatibility,
56(5), 2014. DOI: 10.1109/TEMC.2014.2313977.
=#
function (::Functor{:AlipioVisacro2014})(
        material::EarthMaterial{T},
        frequency::T,
        values::NamedTuple
) where {T <: Real}
    gamma = convert(T, values.exponent)
    epsilon_infinity = convert(T, values.epsilon_infinity)
    scale = convert(T, values.scale)
    conductivity_exponent = convert(T, values.conductivity_exponent)
    thousand = convert(T, 1000)
    million = convert(T, 1e6)
    conductivity_reference = thousand / material.rho
    dispersion = scale * conductivity_reference^conductivity_exponent
    epsilon0 = vacuum_permittivity(frequency)
    relative_permittivity = epsilon_infinity +
                            tan((one(frequency) * π) * gamma / 2) *
                            convert(T, 1e-3) * conductivity_reference * dispersion *
                            frequency^(gamma - one(gamma)) /
                            (
                                2 * (one(frequency) * π) * epsilon0 *
                                convert(T, 10)^(6 * gamma)
                            )
    conductivity = conductivity_reference *
                   (one(frequency) + dispersion * (frequency / million)^gamma)
    return EarthMaterial{T}(
        thousand / conductivity,
        relative_permittivity,
        material.mu_r
    )
end

# Return the stable discovery identifier.
:AlipioVisacro2014
