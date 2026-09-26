"Return the fitted parameters of the Alipio–Visacro causal soil model."
assumptions(::Val{:alipio2014}) = (
    exponent = 0.54,
    epsilon_infinity = 12.0,
    scale = 1.26,
    conductivity_exponent = -0.73
)

function validate(selected::Formula{:alipio2014})
    exponent = selected.parameters.exponent
    isinteger(exponent) && isodd(exponent) && throw(ArgumentError(
        "Alipio–Visacro exponent must not be an odd integer (tangent pole)"))
    return selected
end

"""
$(TYPEDSIGNATURES)

Alipio–Visacro causal power-law soil dispersion fitted from measured data.
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
"""
function description(::Type{<:Formula{:alipio2014}}; compact::Bool=false)
    compact ? "Alipio" : "Alipio–Visacro causal soil dispersion (2014)"
end

function earth_material(
        ::Formula{:alipio2014}, material::EarthMaterial{T}, frequency::T,
        values::NamedTuple, options::FormulationOptions, workspace
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
                            (2 * (one(frequency) * π) * epsilon0 *
                             convert(T, 10)^(6 * gamma))
    conductivity = conductivity_reference *
                   (one(frequency) + dispersion * (frequency / million)^gamma)
    return EarthMaterial{T}(thousand / conductivity, relative_permittivity, material.mu_r)
end

formulation_options(::FormulaMethod{<:Formula{:alipio2014}, typeof(earth_material)}) = FormulationOptions()

:alipio2014
