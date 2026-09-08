
"""
$(TYPEDSIGNATURES)

**Identification.** Explicit default lossless semiconducting-screen approximation.

**Expression.** ``\\kappa=j\\omega\\varepsilon_0\\varepsilon_r``. Retain
the screen geometry and permittivity but suppress conductivity and loss tangent.

**Reference.** Lossless dielectric approximation to Maxwell's constitutive
relation; request an explicit lossy law to represent semicon conduction.
"""
description(::Formula{:default}) = "Default lossless semiconducting-screen admittivity"

"""
$(TYPEDSIGNATURES)

Treat a semiconducting screen as a lossless dielectric under the default
constitutive choice:

```math
\\kappa_s=j\\omega\\varepsilon_0\\varepsilon_{r,s}.
```

This explicit approximation retains permittivity and suppresses conductivity
and dielectric loss tangent. Select a lossy relation explicitly to retain loss.

# Arguments

- `material`: Semiconducting material and relative permittivity.
- `frequency`: Evaluation frequency \\[Hz\\].
- `temperature`: Orchestration-supplied temperature \\[°C\\]; this relation
  applies no temperature correction.
- `values`: Empty physical-parameter tuple.
- `options`: Empty numerical sections for this equation.
- `workspace`: Optional execution resources.

# Returns

- Complex lossless admittivity \\[S/m\\].
"""
@inline function semicon_material(
        ::Val{:default}, material::Material{T}, frequency::T,
        temperature::T, values::NamedTuple, options::NamedTuple, workspace
) where {T <: Real}
    ε₀ = one(T) * 88541878128 * (one(T) * 10)^(-22)
    ω = 2 * (one(T) * π) * frequency
    return complex(zero(T), ω) * ε₀ * material.eps_r
end

computation_options(::FormulaMethod{:default, typeof(semicon_material)}) = (;)

:default
