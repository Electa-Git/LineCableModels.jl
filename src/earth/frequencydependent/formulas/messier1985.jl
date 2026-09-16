"Return the high-frequency permittivity parameter of the Messier relation."
assumptions(::Val{:messier1985}) = (epsilon_infinity = 8.0,)

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

**Reference.** M. Messier, *Another Soil Conductivity Model*, JAYCOR,
Santa Barbara, 1985.
"""
function description(::Type{<:Formula{:messier1985}}; compact::Bool=false)
    compact ? "Messier" : "Messier square-root soil dispersion (1985)"
end

function earth_material(
        ::Formula{:messier1985}, material::EarthMaterial{T}, frequency::T,
        values::NamedTuple, options::FormulationOptions, workspace
) where {T <: Real}
    conductivity_reference = inv(material.rho)
    epsilon_infinity = convert(T, values.epsilon_infinity)
    epsilon0 = vacuum_permittivity(frequency)
    pi_typed = one(frequency) * π
    relative_permittivity = epsilon_infinity + sqrt(
        conductivity_reference * epsilon_infinity / (pi_typed * frequency * epsilon0)
    )
    conductivity = conductivity_reference + sqrt(
        4 * pi_typed * frequency * conductivity_reference * epsilon0 * epsilon_infinity
    )
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

formulation_options(::FormulaMethod{<:Formula{:messier1985}, typeof(earth_material)}) = FormulationOptions()

:messier1985
