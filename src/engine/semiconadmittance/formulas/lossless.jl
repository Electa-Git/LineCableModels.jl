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

- `material`: Semiconducting material and relative permittivity.
- `frequency`: Evaluation frequency \\[Hz\\].
- `temperature`: Operating temperature \\[°C\\].
- `values`: Explicit physical/model parameters.
- `options`: Normalized numerical sections for this contribution.
- `workspace`: Optional execution resources.

# Returns

- Complex lossless admittivity \\[S/m\\].
"""
@inline function semicon_material(
        ::Formula{:lossless}, material::Material{T}, frequency::T,
        temperature::T, values::NamedTuple, options::NamedTuple, workspace
) where {T <: Real}
    ε₀ = one(T) * 88541878128 * (one(T) * 10)^(-22)
    ω = 2 * (one(T) * π) * frequency
    return complex(zero(T), ω) * ε₀ * material.eps_r
end

computation_options(::FormulaMethod{<:Formula{:lossless}, typeof(semicon_material)}) = (;)

:lossless
