"Return the fitted parameters of the Portela soil-dispersion relation."
assumptions(::Val{:portela1999}) = (beta = 0.1, exponent = 0.72)

function validate(selected::Formula{:portela1999})
    exponent = selected.parameters.exponent
    isinteger(exponent) && isodd(exponent) && throw(ArgumentError(
        "Portela exponent must not be an odd integer (tangent pole)"))
    return selected
end

"""
$(TYPEDSIGNATURES)

Portela causal power-law soil-dispersion model.

**Expression.** For ``\\omega=2\\pi f``, ``\\beta=0.1``, and
``\\alpha=0.72``,

```math
\\sigma(f)=\\sigma_0+\\beta10^{-6}\\omega^\\alpha,
\\qquad
\\varepsilon_r(f)=\\frac{\\beta10^{-6}\\tan(\\pi\\alpha/2)
\\omega^{\\alpha-1}}{\\varepsilon_0}.
```

**Reference.** C. M. Portela, “Measurement and Modeling of Soil
Electromagnetic Behavior,” *IEEE International Symposium on Electromagnetic
Compatibility*, 1004–1009, 1999.
"""
function description(::Type{<:Formula{:portela1999}}; compact::Bool=false)
    compact ? "Portela" : "Portela power-law soil dispersion (1999)"
end

function earth_material(
        ::Formula{:portela1999}, material::EarthMaterial{T}, frequency::T,
        values::NamedTuple, options::NamedTuple, workspace
) where {T <: Real}
    beta = convert(T, values.beta)
    exponent = convert(T, values.exponent)
    angular_frequency = 2 * (one(frequency) * π) * frequency
    fitted_scale = beta * convert(T, 1e-6)
    conductivity = inv(material.rho) + fitted_scale * angular_frequency^exponent
    relative_permittivity = fitted_scale *
                            tan((one(frequency) * π) * exponent / 2) *
                            angular_frequency^(exponent - one(exponent)) /
                            vacuum_permittivity(frequency)
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

computation_options(::FormulaMethod{<:Formula{:portela1999}, typeof(earth_material)}) = (;)

:portela1999
