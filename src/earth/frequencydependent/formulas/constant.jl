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

- `material`: Static earth material.
- `frequency`: Evaluation frequency \\[Hz\\].
- `parameters`: Physical parameters of the selected relation.
- `options`: Normalized numerical sections for this contribution.
- `workspace`: Optional execution resources.

# Returns

- The unchanged `EarthMaterial`.
"""
function earth_material(
        ::Formula{:constant}, material::EarthMaterial, frequency::Real,
        values::NamedTuple, options::FormulationOptions, workspace
)
    material
end

formulation_options(::FormulaMethod{<:Formula{:constant}, typeof(earth_material)}) = FormulationOptions()

:constant
