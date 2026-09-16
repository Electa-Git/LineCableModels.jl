"Return the fitted parameters recommended by CIGRE Technical Brochure 781."
assumptions(::Val{:cigre2019}) = (
    epsilon_infinity = 12.0,
    epsilon_scale = 9.5e4,
    epsilon_conductivity_exponent = 0.27,
    epsilon_frequency_exponent = -0.46,
    conductivity_scale = 4.7e-6,
    conductivity_frequency_exponent = 0.54
)

"""
$(TYPEDSIGNATURES)

CIGRE WG C4.33 recommended empirical soil-dispersion relation.

**Expression.** With ``\\sigma_0=1/\\rho_0``,

```math
\\varepsilon_r(f)=12+9.5\\times10^4\\sigma_0^{0.27}f^{-0.46},
\\qquad
\\sigma(f)=\\sigma_0+4.7\\times10^{-6}\\sigma_0^{0.27}f^{0.54}.
```

The fitted constants are exposed through the formula parameters.

**Reference.** CIGRE WG C4.33, Technical Brochure 781 (2019).
"""
function description(::Type{<:Formula{:cigre2019}}; compact::Bool=false)
    compact ? "CIGRE" : "CIGRE WG C4.33 recommended soil dispersion (2019)"
end

function earth_material(
        ::Formula{:cigre2019}, material::EarthMaterial{T}, frequency::T,
        values::NamedTuple, options::FormulationOptions, workspace
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

formulation_options(::FormulaMethod{<:Formula{:cigre2019}, typeof(earth_material)}) = FormulationOptions()

:cigre2019
