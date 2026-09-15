"Return the normalization parameter of the Visacro–Portela soil relation."
assumptions(::Val{:visacro1987}) = (normalization_frequency = 100.0,)

function validate(selected::Formula{:visacro1987})
    selected.parameters.normalization_frequency > 0 || throw(ArgumentError(
        "Visacro–Portela normalization_frequency must be positive [Hz]"))
    return selected
end

"""
$(TYPEDSIGNATURES)

Visacro–Portela empirical power laws for soil conductivity and permittivity.

**Expression.** With ``\\sigma_0=1/\\rho_0``,

```math
\\varepsilon_r(f)=2.34\\times10^6\\sigma_0^{0.535}f^{-0.597},
\\qquad
\\sigma(f)=\\sigma_0\\left(\\frac{f}{100}\\right)^{0.072}.
```

**Reference.** S. Visacro and C. M. Portela, *International Symposium on
High Voltage Engineering*, 1987.
"""
function description(::Type{<:Formula{:visacro1987}}; compact::Bool=false)
    compact ? "Visacro" : "Visacro–Portela empirical soil dispersion (1987)"
end

function earth_material(
        ::Formula{:visacro1987}, material::EarthMaterial{T}, frequency::T,
        values::NamedTuple, options::NamedTuple, workspace
) where {T <: Real}
    conductivity_reference = inv(material.rho)
    relative_permittivity = convert(T, 2.34e6) *
                            conductivity_reference^convert(T, 0.535) *
                            frequency^convert(T, -0.597)
    conductivity = conductivity_reference *
                   (frequency / convert(T, values.normalization_frequency))^convert(T, 0.072)
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

computation_options(::FormulaMethod{<:Formula{:visacro1987}, typeof(earth_material)}) = (;)

:visacro1987
