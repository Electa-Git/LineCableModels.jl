assumptions(::Val{:default}) = (;)

"""
$(TYPEDSIGNATURES)

**Identification.** Explicit default lossless cable-insulation approximation.

**Expression.** ``\\kappa=j\\omega\\varepsilon_0\\varepsilon_r``. Retain
displacement current and suppress both material conductivity and loss tangent.

**Reference.** Lossless dielectric approximation to Maxwell's constitutive
relation; this default is not an author-labelled loss model.
"""
description(::Formula{:default}) = "Default lossless cable-insulation admittivity"

"""
$(TYPEDSIGNATURES)

Evaluate the explicit default lossless dielectric relation:

```math
\\kappa=j\\omega\\varepsilon_0\\varepsilon_r.
```

Material conductivity and dielectric loss tangent are suppressed; displacement
current is retained. This choice does not suppress conductor or earth losses.

# Arguments

- `material`: Insulation material and relative permittivity.
- `frequency`: Evaluation frequency \\[Hz\\].
- `temperature`: Orchestration-supplied temperature \\[°C\\]; this relation
  applies no temperature correction.
- `values`: Empty assumption tuple.

# Returns

- Complex lossless admittivity \\[S/m\\].
"""
@inline function insulation_material(
        ::Val{:default}, material::Material{T}, frequency::T,
        temperature::T, values::NamedTuple
) where {T <: Real}
    ε₀ = one(T) * 88541878128 * (one(T) * 10)^(-22)
    ω = 2 * (one(T) * π) * frequency
    return complex(zero(T), ω) * ε₀ * material.eps_r
end

:default
