"""
$(TYPEDSIGNATURES)

**Identification.** Generic lossy complex-admittivity representation of a
homogeneous cable-insulation layer. Ohmic conduction, dielectric displacement,
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
Miyamoto, and Nagaoka (2004), Eqs. (14)-(15), remain a useful cable-layer
application reference. The paper's main-insulation term is the lossless
specialization.
"""
function description(::Type{<:Formula{:lossy}}; compact::Bool=false)
    compact ? "Lossy" : "Lossy complex-admittivity cable-insulation model"
end

"""
$(TYPEDSIGNATURES)

Evaluate the standard lossy material admittivity for a cable-insulation layer:

```math
\\kappa=\\frac{1}{\\rho}+\\omega\\varepsilon_0\\varepsilon_r\\tan\\delta_p+
j\\omega\\varepsilon_0\\varepsilon_r.
```

`material.tan_delta` represents polarization loss only. Conduction is supplied
by `material.rho`. The common coaxial operator applies the annular geometry.

# Arguments

- `functor`: the Functor of the evaluation point. Its input holds:
  - `material`: insulation material properties, including resistivity and
    relative permittivity.
  - `frequency`: evaluation frequency \\[Hz\\].
  - `temperature`: operating temperature \\[°C\\].
  - `options`: the normalized formulation options of the formula.
- `workspace`: optional computation workspace supplying reusable numerical buffers.

# Returns

- Complex material admittivity \\[S/m\\].

# Notes

Ametani, Miyamoto, and Nagaoka (2004), DOI
10.1109/TPWRD.2003.822502, is retained as a cable-layer application
reference. The constitutive relation itself is standard frequency-domain
electromagnetism.
"""
@inline function insulation_material(::Formula{:lossy}, functor, workspace)
    (; material, frequency) = functor.input
    T = typeof(frequency)
    ε₀ = vacuum_permittivity(T)
    ω = 2 * (one(T) * π) * frequency
    displacement = complex(zero(T), ω) * ε₀ * material.eps_r
    return conductivity(material.rho) + imag(displacement) * material.tan_delta +
           displacement
end

formulation_options(::Expression{<:Formula{:lossy}, typeof(insulation_material)}) = FormulationOptions()

:lossy
