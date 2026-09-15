"""
$(TYPEDSIGNATURES)

**Identification.** Lossless cable-insulation approximation retaining
displacement current while suppressing conduction and polarization loss.

**Expression.**

```math
\\kappa=j\\omega\\varepsilon_0\\varepsilon_r.
```

This is the lossless specialization of the standard frequency-domain
constitutive relation.
"""
function description(::Type{<:Formula{:lossless}}; compact::Bool=false)
    compact ? "Lossless" : "Lossless cable-insulation admittivity"
end

"""
$(TYPEDSIGNATURES)

Evaluate lossless cable-insulation admittivity:

```math
\\kappa=j\\omega\\varepsilon_0\\varepsilon_r.
```

# Arguments

- `material`: Insulation material and relative permittivity.
- `frequency`: Evaluation frequency \\[Hz\\].
- `temperature`: Operating temperature \\[°C\\].
- `values`: Explicit physical/model parameters.
- `options`: Normalized numerical sections for this contribution.
- `workspace`: Optional execution resources.

# Returns

- Complex lossless admittivity \\[S/m\\].
"""
@inline function insulation_material(
        ::Val{:lossless}, material::Material{T}, frequency::T,
        temperature::T, values::NamedTuple, options::NamedTuple, workspace
) where {T <: Real}
    ε₀ = one(T) * 88541878128 * (one(T) * 10)^(-22)
    ω = 2 * (one(T) * π) * frequency
    return complex(zero(T), ω) * ε₀ * material.eps_r
end

computation_options(::FormulaMethod{:lossless, typeof(insulation_material)}) = (;)

:lossless
