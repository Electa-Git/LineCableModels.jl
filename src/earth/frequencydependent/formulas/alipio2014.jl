# Construct `:alipio2014` with the fitted parameters of the Alipio-Visacro causal soil
# model as parameter defaults.
function Formula{:alipio2014}(; parameters::NamedTuple = (;),
        options::Union{NamedTuple, FormulationOptions} = FormulationOptions())
    defaults = (exponent = 0.54, epsilon_infinity = 12.0, scale = 1.26,
        conductivity_exponent = -0.73)
    return Formula{:alipio2014}(defaults, parameters, options)
end

function validate(selected::Formula{:alipio2014})
    exponent = selected.parameters.exponent
    isinteger(exponent) && isodd(exponent) && throw(ArgumentError(
        "Alipio–Visacro exponent must not be an odd integer (tangent pole)"))
    return selected
end

"""
$(TYPEDSIGNATURES)

Causal power law fitted by Alipio and Visacro to measured soil dispersion.
The relation evaluates conductivity and relative permittivity from the static
reference resistivity. Ametani is not involved in this model.

**Expression.** With ``\\widehat\\sigma_0=1000/\\rho_0`` in mS/m,
``\\gamma=0.54``, and ``D=1.26\\widehat\\sigma_0^{-0.73}``,

```math
\\widehat\\sigma(f)=\\widehat\\sigma_0
\\left[1+D\\left(\\frac{f}{10^6}\\right)^\\gamma\\right],
```

```math
\\varepsilon_r(f)=12+
\\tan\\left(\\frac{\\pi\\gamma}{2}\\right)
\\frac{10^{-3}\\widehat\\sigma_0D f^{\\gamma-1}}
{2\\pi\\varepsilon_0\\cdot10^{6\\gamma}}.
```

With ``\\sigma_0=1/\\rho_0`` and ``\\varepsilon_\\infty=12``,
the analytic continuation at complex frequency ``s`` in 1/s is

```math
\\kappa(s)=\\sigma_0+s\\varepsilon_0\\varepsilon_\\infty+\\frac{\\sigma_0D}{\\cos(\\pi\\gamma/2)}\\left(\\frac{s}{2\\pi10^6}\\right)^\\gamma.
```

The authors' causal model ``\\kappa(s)`` is conjugate-symmetric.
The returned real pair satisfies
``\\varepsilon_{r,s}=\\operatorname{Im}\\kappa(s)/(\\operatorname{Im}s\\,\\varepsilon_0)``
and ``1/\\rho_s=\\operatorname{Re}\\kappa(s)-\\operatorname{Re}s\\,\\varepsilon_0\\varepsilon_{r,s}``.

"""
function description(::Type{<:Formula{:alipio2014}}; compact::Bool=false)
    compact ? "Alipio" : "Alipio–Visacro causal soil dispersion (2014)"
end

function earth_material(formula::Formula{:alipio2014}, functor, workspace)
    (; material, frequency, s) = functor.input
    T = typeof(frequency)
    parameters = formula.parameters
    gamma = convert(T, parameters.exponent)
    epsilon_infinity = convert(T, parameters.epsilon_infinity)
    scale = convert(T, parameters.scale)
    conductivity_exponent = convert(T, parameters.conductivity_exponent)
    thousand = convert(T, 1000)
    million = convert(T, 1e6)
    conductivity_reference = thousand / material.rho
    dispersion = scale * conductivity_reference^conductivity_exponent
    epsilon0 = vacuum_permittivity(typeof(frequency))
    conductivity_reference /= thousand
    κ = conductivity_reference + s * epsilon0 * epsilon_infinity +
        (conductivity_reference * dispersion / cos((one(T) * π) * gamma / 2)) *
        (s / (2 * (one(T) * π) * million))^gamma
    relative_permittivity = imag(κ) / (imag(s) * epsilon0)
    conductivity = real(κ) - real(s) * epsilon0 * relative_permittivity
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

formulation_options(::Expression{<:Formula{:alipio2014}, typeof(earth_material)}) = FormulationOptions()

:alipio2014
