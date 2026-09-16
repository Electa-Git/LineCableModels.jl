const LONGMIRE_SMITH_COEFFICIENTS = (
    3.4e6, 2.74e5, 2.58e4, 3.38e3, 5.26e2, 1.33e2, 2.72e1,
    1.25e1, 4.8, 2.17, 0.98, 0.392, 0.173
)

"Return the high-frequency and relaxation parameters of Longmire–Smith."
assumptions(::Val{:longmire1975}) = (
    epsilon_infinity = 5.0,
    corner_scale = 125.0,
    corner_exponent = 0.8312
)

function validate(selected::Formula{:longmire1975})
    selected.parameters.corner_scale > 0 || throw(ArgumentError(
        "Longmire–Smith corner_scale must be positive"))
    return selected
end

"""
$(TYPEDSIGNATURES)

Longmire–Smith thirteen-term dielectric-relaxation soil model.

**Expression.** With ``\\sigma_0=1/\\rho_0``, base corner
``f_c=(125\\sigma_0)^{0.8312}``, tabulated coefficients ``a_n``, and
``f_n=10^{n-1}f_c``,

```math
\\varepsilon_r(f)=5+\\sum_{n=1}^{13}
\\frac{a_n}{1+(f/f_n)^2},
```

```math
\\sigma(f)=\\sigma_0+2\\pi f\\varepsilon_0
\\sum_{n=1}^{13}\\frac{a_n(f/f_n)}{1+(f/f_n)^2}.
```

**Reference.** C. L. Longmire and K. S. Smith, *A Universal Impedance for
Soils*, Defense Nuclear Agency, 1975.
"""
function description(::Type{<:Formula{:longmire1975}}; compact::Bool=false)
    compact ? "Longmire" : "Longmire–Smith 13-term dielectric relaxation (1975)"
end

function earth_material(
        ::Formula{:longmire1975}, material::EarthMaterial{T}, frequency::T,
        values::NamedTuple, options::FormulationOptions, workspace
) where {T <: Real}
    conductivity_reference = inv(material.rho)
    corner = (convert(T, values.corner_scale) * conductivity_reference)^convert(T, values.corner_exponent)
    permittivity_sum = zero(frequency)
    conductivity_sum = zero(frequency)
    decade = one(frequency)
    ten = convert(T, 10)
    @inbounds for coefficient in LONGMIRE_SMITH_COEFFICIENTS
        corner_frequency = corner * decade
        ratio = frequency / corner_frequency
        denominator = one(frequency) + ratio^2
        typed_coefficient = convert(T, coefficient)
        permittivity_sum += typed_coefficient / denominator
        conductivity_sum += typed_coefficient * ratio / denominator
        decade *= ten
    end
    relative_permittivity = convert(T, values.epsilon_infinity) + permittivity_sum
    conductivity = conductivity_reference +
                   2 * (one(frequency) * π) * frequency *
                   vacuum_permittivity(frequency) * conductivity_sum
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

formulation_options(::FormulaMethod{<:Formula{:longmire1975}, typeof(earth_material)}) = FormulationOptions()

:longmire1975
