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

- `functor`: the Functor of the evaluation point. Its input holds:
  - `material`: insulation material and relative permittivity.
  - `frequency`: evaluation frequency \\[Hz\\].
  - `temperature`: operating temperature \\[°C\\].
  - `options`: the normalized formulation options of the formula.
- `workspace`: optional computation workspace supplying reusable numerical buffers.

# Returns

- Complex lossless admittivity \\[S/m\\].
"""
@inline function insulation_material(::Formula{:lossless}, functor, workspace)
    (; material, frequency) = functor.input
    T = typeof(frequency)
    ε₀ = vacuum_permittivity(T)
    ω = 2 * (one(T) * π) * frequency
    return complex(zero(T), ω) * ε₀ * material.eps_r
end

formulation_options(::Expression{<:Formula{:lossless}, typeof(insulation_material)}) = FormulationOptions()

:lossless
