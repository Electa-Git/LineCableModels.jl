"""
$(TYPEDSIGNATURES)

**Identification.** Generic lossy complex-admittivity representation of a
homogeneous semiconducting screen. Ohmic conduction, dielectric displacement,
and an optional polarization-loss contribution are evaluated together.

**Expression.** For resistivity ``\\rho``, real relative permittivity
``\\varepsilon_r``, and polarization loss tangent ``\\tan\\delta_p``,

```math
\\kappa=\\frac{1}{\\rho}+\\omega\\varepsilon_0\\varepsilon_r\\tan\\delta_p+
j\\omega\\varepsilon_0\\varepsilon_r.
```

The annular layer operator converts this material admittivity to
``Y=2\\pi\\kappa/\\ln(b/a)``. This is the standard frequency-domain
constitutive relation, not an author-specific empirical law. Ametani,
Miyamoto, and Nagaoka (2004), Eqs. (14)–(15), remain a useful application
reference for its zero-``\\tan\\delta_p`` semiconducting-screen specialization
and radial series assembly.
"""
function description(::Type{<:Formula{:lossy}}; compact::Bool=false)
    compact ? "Lossy" : "Lossy complex-admittivity semiconducting-screen model"
end

"""
$(TYPEDSIGNATURES)

Evaluate the standard lossy material admittivity for a semiconducting layer:

```math
\\kappa=\\frac{1}{\\rho}+\\omega\\varepsilon_0\\varepsilon_r\\tan\\delta_p+
j\\omega\\varepsilon_0\\varepsilon_r.
```

`material.tan_delta` represents polarization loss only; conduction is supplied
by `material.rho`. The common coaxial operator applies the annular geometry.

# Arguments

- `material`: Semiconducting material properties, including resistivity and
  relative permittivity.
- `frequency`: Evaluation frequency \\[Hz\\].
- `temperature`: Operating temperature \\[°C\\].
- `values`: Explicit physical/model parameters.
- `options`: Normalized numerical sections for this contribution.
- `workspace`: Optional computation workspace supplying reusable numerical buffers.

# Returns

- Complex screen admittivity \\[S/m\\].

# Notes

Ametani, Miyamoto, and Nagaoka (2004), DOI
10.1109/TPWRD.2003.822502, is retained as an application reference for the
semiconducting-screen specialization; the constitutive relation itself is
standard frequency-domain electromagnetism.
"""
@inline function semicon_material(
        ::Formula{:lossy},
        material::Material{T},
        frequency::T,
        temperature::T,
        values::NamedTuple, options::FormulationOptions, workspace
) where {T <: Real}
    ε₀ = one(T) * 88541878128 * (one(T) * 10)^(-22)
    ω = 2 * (one(T) * π) * frequency
    displacement = complex(zero(T), ω) * ε₀ * material.eps_r
    return conductivity(material.rho) + imag(displacement) * material.tan_delta +
           displacement
end

formulation_options(::FormulaMethod{<:Formula{:lossy}, typeof(semicon_material)}) = FormulationOptions()

:lossy
