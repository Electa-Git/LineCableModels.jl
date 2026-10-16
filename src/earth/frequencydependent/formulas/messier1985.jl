# Construct `:messier1985` with the high-frequency permittivity parameter of the Messier
# relation as parameter defaults.
function Formula{:messier1985}(; parameters::NamedTuple = (;),
        options::Union{NamedTuple, FormulationOptions} = FormulationOptions())
    defaults = (epsilon_infinity = 8.0,)
    return Formula{:messier1985}(defaults, parameters, options)
end

function validate(selected::Formula{:messier1985})
    selected.parameters.epsilon_infinity >= 0 || throw(ArgumentError(
        "Messier epsilon_infinity must be nonnegative for the real square roots"))
    return selected
end

"""
$(TYPEDSIGNATURES)

Messier square-root dispersive soil model.

**Expression.** With ``\\sigma_0=1/\\rho_0`` and
``\\varepsilon_\\infty=8``,

```math
\\varepsilon_r(f)=\\varepsilon_\\infty+
\\sqrt{\\frac{\\sigma_0\\varepsilon_\\infty}{\\pi f\\varepsilon_0}},
\\qquad
\\sigma(f)=\\sigma_0+
\\sqrt{4\\pi f\\sigma_0\\varepsilon_0\\varepsilon_\\infty}.
```

For complex frequency ``s`` in 1/s, the analytic continuation is

```math
\\kappa(s)=\\sigma_0+s\\varepsilon_0\\varepsilon_\\infty+2\\sqrt{s\\sigma_0\\varepsilon_0\\varepsilon_\\infty}.
```

The authors' causal model ``\\kappa(s)`` is conjugate-symmetric.
The returned real pair satisfies
``\\varepsilon_{r,s}=\\operatorname{Im}\\kappa(s)/(\\operatorname{Im}s\\,\\varepsilon_0)``
and ``1/\\rho_s=\\operatorname{Re}\\kappa(s)-\\operatorname{Re}s\\,\\varepsilon_0\\varepsilon_{r,s}``.

**Reference.** M. Messier, *Another Soil Conductivity Model*, JAYCOR,
Santa Barbara, 1985.
"""
function description(::Type{<:Formula{:messier1985}}; compact::Bool=false)
    compact ? "Messier" : "Messier square-root soil dispersion (1985)"
end

function earth_material(formula::Formula{:messier1985}, functor, workspace)
    (; material, frequency, s) = functor.input
    T = typeof(frequency)
    parameters = formula.parameters
    conductivity_reference = inv(material.rho)
    epsilon_infinity = convert(T, parameters.epsilon_infinity)
    epsilon0 = vacuum_permittivity(typeof(frequency))
    κ = conductivity_reference + s * epsilon0 * epsilon_infinity +
        2 * sqrt(s * conductivity_reference * epsilon0 * epsilon_infinity)
    relative_permittivity = imag(κ) / (imag(s) * epsilon0)
    conductivity = real(κ) - real(s) * epsilon0 * relative_permittivity
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

formulation_options(::Expression{<:Formula{:messier1985}, typeof(earth_material)}) = FormulationOptions()

:messier1985
