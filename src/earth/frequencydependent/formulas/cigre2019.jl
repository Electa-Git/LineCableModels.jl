# Construct `:cigre2019` with the fitted parameters recommended by CIGRE Technical
# Brochure 781 as parameter defaults.
function Formula{:cigre2019}(; parameters::NamedTuple = (;),
        options::Union{NamedTuple, FormulationOptions} = FormulationOptions())
    defaults = (epsilon_infinity = 12.0, epsilon_scale = 9.5e4,
        epsilon_conductivity_exponent = 0.27, epsilon_frequency_exponent = -0.46,
        conductivity_scale = 4.7e-6, conductivity_frequency_exponent = 0.54)
    return Formula{:cigre2019}(defaults, parameters, options)
end

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

function earth_material(formula::Formula{:cigre2019}, functor, workspace)
    (; material, frequency) = functor.input
    T = typeof(frequency)
    parameters = formula.parameters
    conductivity_reference = inv(material.rho)
    conductivity_exponent = convert(T, parameters.epsilon_conductivity_exponent)
    relative_permittivity = convert(T, parameters.epsilon_infinity) +
                            convert(T, parameters.epsilon_scale) *
                            conductivity_reference^conductivity_exponent *
                            frequency^convert(T, parameters.epsilon_frequency_exponent)
    conductivity = conductivity_reference +
                   convert(T, parameters.conductivity_scale) *
                   conductivity_reference^conductivity_exponent *
                   frequency^convert(T, parameters.conductivity_frequency_exponent)
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

formulation_options(::Expression{<:Formula{:cigre2019}, typeof(earth_material)}) = FormulationOptions()

:cigre2019
