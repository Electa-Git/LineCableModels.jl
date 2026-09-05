assumptions(::Val{:Ametani2004}) = (;)

"""
$(TYPEDSIGNATURES)

**Identification.** Complex-permittivity representation of a concentric
semiconducting screen. Its conduction and displacement currents enter the same
radial dielectric network as the adjacent insulation layers.

**Expression.** For screen resistivity ``\\rho_s`` and permittivity
``\\varepsilon_s``,

```math
\\varepsilon_s^\\star=\\varepsilon_s+\\frac{1}{j\\omega\\rho_s},\\qquad
\\kappa_s=\\frac{1}{\\rho_s}+j\\omega\\varepsilon_s,
```

and an annular screen from ``a`` to ``b`` has

```math
Y_s=\\frac{2\\pi\\kappa_s}{\\ln(b/a)}.
```

**Reference.** A. Ametani, Y. Miyamoto, and N. Nagaoka, 2004, as reproduced
in A. Ametani, T. Ohno, and N. Nagaoka, *Cable System Transients: Theory,
Modeling and Simulation*, Wiley-IEEE Press, 2015, Eqs. 2.66–2.67.
"""
function description(::Formula{:Ametani2004})
    "Ametani semiconducting-screen admittance model (2004)"
end

"""
$(TYPEDSIGNATURES)

Retain semiconducting-screen conductivity and permittivity in Ametani's
complex-permittivity representation:

```math
\\varepsilon_s^\\star=\\varepsilon_s+\\frac{1}{j\\omega\\rho_s},
\\qquad
\\kappa_s=\\frac{1}{\\rho_s}+j\\omega\\varepsilon_s.
```

The common Coaxial Engine operator applies the annular geometry and combines
the semiconducting screen with adjacent dielectric layers radially in series.

# Arguments

- `material`: Static semiconducting-screen properties.
- `frequency`: Evaluation frequency \\[Hz\\].
- `temperature`: Operating temperature \\[°C\\].
- `values`: Formula assumptions.

# Returns

- Complex screen admittivity ``1/\\rho_s+j\\omega\\varepsilon_s`` \\[S/m\\].

# References

A. Ametani, Y. Miyamoto, and N. Nagaoka, 2004, as reproduced in Ametani,
Ohno, and Nagaoka (2015), Eqs. 2.66–2.67.
"""
@inline function semicon_material(
        ::Val{:Ametani2004},
        material::Material{T},
        frequency::T,
        temperature::T,
        values::NamedTuple
) where {T <: Real}
    ε₀ = one(T) * 88541878128 * (one(T) * 10)^(-22)
    ω = 2 * (one(T) * π) * frequency
    return conductivity(material.rho) +
           complex(zero(T), ω) * ε₀ * material.eps_r
end

:Ametani2004
