"""
$(TYPEDSIGNATURES)

**Identification.** Generic lossy complex-admittivity representation of a
homogeneous semiconducting screen. Ohmic conduction, dielectric displacement,
and an optional polarization-loss contribution are evaluated together.

**Expression.** For resistivity ``\\rho``, real relative permittivity
``\\varepsilon_r``, and polarization loss tangent ``\\tan\\delta_p``,

```math
\\kappa(j\\omega)=\\frac{1}{\\rho}+\\omega\\varepsilon_0\\varepsilon_r\\tan\\delta_p+j\\omega\\varepsilon_0\\varepsilon_r.
```

At σ ≠ 0, ``\\kappa(s)=1/\\rho+s\\varepsilon_0\\varepsilon_r``.
A constant loss tangent is not causal and has no value at a complex frequency.
A computation at σ ≠ 0 refuses a nonzero `material.tan_delta`.

The annular layer operator converts this material admittivity to
``Y=2\\pi\\kappa/\\ln(b/a)``. This is the standard frequency-domain
constitutive relation, not an author-specific empirical law. Ametani,
Miyamoto, and Nagaoka (2004), Eqs. (14)-(15), remain a useful application
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
\\kappa(j\\omega)=\\frac{1}{\\rho}+\\omega\\varepsilon_0\\varepsilon_r\\tan\\delta_p+j\\omega\\varepsilon_0\\varepsilon_r.
```

At σ ≠ 0, ``\\kappa(s)=1/\\rho+s\\varepsilon_0\\varepsilon_r``.
A constant loss tangent is not causal and has no value at a complex frequency.
A computation at σ ≠ 0 refuses a nonzero `material.tan_delta`.

`material.tan_delta` represents polarization loss only. Conduction is supplied
by `material.rho`. The common coaxial operator applies the annular geometry.

# Arguments

- `functor`: the Functor of the evaluation point. Its input holds:
  - `material`: semiconducting material properties, including resistivity and
    relative permittivity.
  - `frequency`: evaluation frequency \\[Hz\\].
  - `s`: complex frequency \\[1/s\\].
  - `temperature`: operating temperature \\[°C\\].
  - `options`: the normalized formulation options of the formula.
- `workspace`: optional computation workspace supplying reusable numerical buffers.

# Returns

- Complex screen admittivity \\[S/m\\].

# Notes

Ametani, Miyamoto, and Nagaoka (2004), DOI
10.1109/TPWRD.2003.822502, is retained as an application reference for the
semiconducting-screen specialization. The constitutive relation itself is
standard frequency-domain electromagnetism.
"""
@inline function semicon_material(formula::Formula{:lossy}, functor, workspace)
    (; material, frequency, s) = functor.input
    validate(s, formula, material)
    T = typeof(frequency)
    ε₀ = vacuum_permittivity(T)
    return conductivity(material.rho) + s * ε₀ * material.eps_r *
           complex(one(T), -material.tan_delta)
end

function validate(s::Complex, ::Formula{:lossy}, material::Material)
    iszero(real(s)) || iszero(material.tan_delta) || throw(ArgumentError(
        "a constant loss tangent is not causal and has no value at a complex frequency; at σ ≠ 0 set tan_delta = 0"))
    return s
end

formulation_options(::Expression{<:Formula{:lossy}, typeof(semicon_material)}) = FormulationOptions()

:lossy
