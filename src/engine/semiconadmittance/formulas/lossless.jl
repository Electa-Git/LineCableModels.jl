"""
$(TYPEDSIGNATURES)

**Identification.** Lossless semiconducting-screen approximation retaining
displacement current while suppressing conduction and polarization loss.

**Expression.**

```math
\\kappa=j\\omega\\varepsilon_0\\varepsilon_r.
```

This is the lossless specialization of the standard frequency-domain
constitutive relation.
"""
function description(::Type{<:Formula{:lossless}}; compact::Bool=false)
    compact ? "Lossless" : "Lossless semiconducting-screen admittivity"
end

"""
$(TYPEDSIGNATURES)

Evaluate lossless semiconducting-screen admittivity:

```math
\\kappa=j\\omega\\varepsilon_0\\varepsilon_r.
```

# Arguments

- `material`: semiconducting material and relative permittivity.
- `frequency`: evaluation frequency \\[Hz\\].
- `temperature`: operating temperature \\[°C\\].
- `values`: explicit physical model parameters.
- `options`: normalized numerical sections for this contribution.
- `workspace`: optional computation workspace supplying reusable numerical buffers.

# Returns

- Complex lossless admittivity \\[S/m\\].
"""
@inline function semicon_material(
        ::Formula{:lossless}, material::Material{T}, frequency::T,
        temperature::T, values::NamedTuple, options::FormulationOptions, workspace
) where {T <: Real}
    ε₀ = one(T) * 88541878128 * (one(T) * 10)^(-22)
    ω = 2 * (one(T) * π) * frequency
    return complex(zero(T), ω) * ε₀ * material.eps_r
end

formulation_options(::FormulaMethod{<:Formula{:lossless}, typeof(semicon_material)}) = FormulationOptions()

:lossless
