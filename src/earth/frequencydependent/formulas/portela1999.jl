# Construct `:portela1999` with the fitted parameters of the Portela soil-dispersion
# relation as parameter defaults.
function Formula{:portela1999}(; parameters::NamedTuple = (;),
        options::Union{NamedTuple, FormulationOptions} = FormulationOptions())
    defaults = (beta = 0.1, exponent = 0.72)
    return Formula{:portela1999}(defaults, parameters, options)
end

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

For complex frequency ``s`` in 1/s, the analytic continuation is

```math
\\kappa(s)=\\sigma_0+\\frac{\\beta10^{-6}}{\\cos(\\pi\\alpha/2)}s^\\alpha.
```

The authors' causal model ``\\kappa(s)`` is conjugate-symmetric.
The returned real pair satisfies
``\\varepsilon_{r,s}=\\operatorname{Im}\\kappa(s)/(\\operatorname{Im}s\\,\\varepsilon_0)``
and ``1/\\rho_s=\\operatorname{Re}\\kappa(s)-\\operatorname{Re}s\\,\\varepsilon_0\\varepsilon_{r,s}``.

**Reference.** C. M. Portela, “Measurement and Modeling of Soil
Electromagnetic Behavior,” *IEEE International Symposium on Electromagnetic
Compatibility*, 1004–1009, 1999.
"""
function description(::Type{<:Formula{:portela1999}}; compact::Bool=false)
    compact ? "Portela" : "Portela power-law soil dispersion (1999)"
end

function earth_material(formula::Formula{:portela1999}, functor, workspace)
    (; material, frequency, s) = functor.input
    T = typeof(frequency)
    parameters = formula.parameters
    beta = convert(T, parameters.beta)
    exponent = convert(T, parameters.exponent)
    epsilon0 = vacuum_permittivity(T)
    fitted_scale = beta * convert(T, 1e-6)
    κ = inv(material.rho) +
        (fitted_scale / cos((one(T) * π) * exponent / 2)) * s^exponent
    relative_permittivity = imag(κ) / (imag(s) * epsilon0)
    conductivity = real(κ) - real(s) * epsilon0 * relative_permittivity
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

formulation_options(::Expression{<:Formula{:portela1999}, typeof(earth_material)}) = FormulationOptions()

:portela1999
