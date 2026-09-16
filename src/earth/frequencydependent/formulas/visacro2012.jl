"Return the lower-frequency boundary of the Visacro–Alipio soil relation."
assumptions(::Val{:visacro2012}) = (frequency_boundary = 100.0,)

function validate(selected::Formula{:visacro2012})
    selected.parameters.frequency_boundary > 0 || throw(ArgumentError(
        "Visacro–Alipio frequency_boundary must be positive [Hz]"))
    return selected
end

"""
$(TYPEDSIGNATURES)

Visacro–Alipio empirical causal soil-dispersion fit with a 100 Hz lower
frequency boundary.

**Expression.** With ``f_e=\\max(f,100)`` and ``\\sigma_0=1/\\rho_0``,

```math
\\varepsilon_r(f)=1.3+7.6\\times10^3 f_e^{-0.4},
\\qquad
\\sigma(f)=\\sigma_0+1.2\\times10^{-6}\\sigma_0^{0.27}(f_e-100)^{0.65}.
```

**Reference.** S. Visacro and R. Alipio, *IEEE Transactions on Power
Delivery*, 27(2), 2012.
"""
function description(::Type{<:Formula{:visacro2012}}; compact::Bool=false)
    compact ? "Visacro" : "Visacro–Alipio empirical soil dispersion (2012)"
end

function earth_material(
        ::Formula{:visacro2012}, material::EarthMaterial{T}, frequency::T,
        values::NamedTuple, options::FormulationOptions, workspace
) where {T <: Real}
    frequency_boundary = convert(T, values.frequency_boundary)
    evaluated_frequency = frequency < frequency_boundary ?
                          frequency_boundary : frequency
    conductivity_reference = inv(material.rho)
    relative_permittivity = convert(T, 1.3) +
                            convert(T, 7.6e3) * evaluated_frequency^convert(T, -0.4)
    conductivity = conductivity_reference +
                   convert(T, 1.2e-6) * conductivity_reference^convert(T, 0.27) *
                   (evaluated_frequency - frequency_boundary)^convert(T, 0.65)
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

formulation_options(::FormulaMethod{<:Formula{:visacro2012}, typeof(earth_material)}) = FormulationOptions()

:visacro2012
