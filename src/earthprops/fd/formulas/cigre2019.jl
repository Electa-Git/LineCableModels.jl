"Return the fitted assumptions recommended by CIGRE Technical Brochure 781."
assumptions(::Val{:CIGRE2019}) = (
    epsilon_infinity = 12.0,
    epsilon_scale = 9.5e4,
    epsilon_conductivity_exponent = 0.27,
    epsilon_frequency_exponent = -0.46,
    conductivity_scale = 4.7e-6,
    conductivity_frequency_exponent = 0.54
)

description(::Formula{:CIGRE2019}) =
    "CIGRE WG C4.33 2019 (recommended soil dispersion)"

#=
Evaluate the frequency-dependent soil relation recommended by CIGRE WG C4.33.

Input resistivity and returned resistivity use Ω·m; conductivity is evaluated
internally in S/m and relative permittivity is dimensionless.

# Reference

CIGRE WG C4.33, *Impact of Soil-Parameter Frequency Dependence on the Response
of Grounding Electrodes and on the Lightning Performance of Electrical
Systems*, Technical Brochure 781, 2019.
=#
function (::Functor{:CIGRE2019})(
        material::EarthMaterial{T},
        frequency::T,
        values::NamedTuple
) where {T <: Real}
    conductivity_reference = inv(material.rho)
    conductivity_exponent = convert(T, values.epsilon_conductivity_exponent)
    relative_permittivity = convert(T, values.epsilon_infinity) +
                            convert(T, values.epsilon_scale) *
                            conductivity_reference^conductivity_exponent *
                            frequency^convert(T, values.epsilon_frequency_exponent)
    conductivity = conductivity_reference +
                   convert(T, values.conductivity_scale) *
                   conductivity_reference^conductivity_exponent *
                   frequency^convert(T, values.conductivity_frequency_exponent)
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

# Return the stable discovery identifier.
:CIGRE2019
