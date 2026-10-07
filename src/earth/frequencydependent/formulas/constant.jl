"""
$(TYPEDSIGNATURES)

**Identification.** Frequency-independent earth-material pass-through.

**Expression.**

```math
\\rho(f)=\\rho_0,\\qquad
\\varepsilon_r(f)=\\varepsilon_{r,0},\\qquad
\\mu_r(f)=\\mu_{r,0}.
```

This relation preserves the static material at every positive evaluation
frequency. It is the explicit equation selected by `:default`.
"""
function description(::Type{<:Formula{:constant}}; compact::Bool=false)
    compact ? "Constant" : "Constant frequency-independent earth material"
end

"""
$(TYPEDSIGNATURES)

Preserve the supplied static earth properties at the requested frequency.

# Arguments

- `functor`: the Functor of the evaluation point. Its input holds:
  - `material`: static earth material.
  - `frequency`: evaluation frequency \\[Hz\\].
  - `options`: the normalized formulation options of the formula.
- `workspace`: optional execution resources.

# Returns

- The unchanged `EarthMaterial`.
"""
function earth_material(::Formula{:constant}, functor, workspace)
    return functor.input.material
end

formulation_options(::Expression{<:Formula{:constant}, typeof(earth_material)}) = FormulationOptions()

:constant
