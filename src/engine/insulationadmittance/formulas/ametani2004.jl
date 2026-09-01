assumptions(::Val{:Ametani2004}) = (;)

"""
$(TYPEDSIGNATURES)

**Identification.** Complex-permittivity representation of a concentric
dielectric layer. Material conduction and displacement current enter the same
radial dielectric network.

**Expression.** For material resistivity ``\\rho`` and permittivity
``\\varepsilon``,

```math
\\varepsilon^\\star=\\varepsilon+\\frac{1}{j\\omega\\rho},\\qquad
\\kappa=\\frac{1}{\\rho}+j\\omega\\varepsilon,
```

and an annular layer from ``a`` to ``b`` has

```math
Y=\\frac{2\\pi\\kappa}{\\ln(b/a)}.
```

**Reference.** A. Ametani, Y. Miyamoto, and N. Nagaoka, 2004, as reproduced
in A. Ametani, T. Ohno, and N. Nagaoka, *Cable System Transients: Theory,
Modeling and Simulation*, Wiley-IEEE Press, 2015, Eqs. 2.66–2.67.
"""
function description(::Formula{:Ametani2004})
    "Ametani complex-permittivity dielectric admittance model (2004)"
end

"""
$(TYPEDSIGNATURES)

Evaluate Ametani's complex-permittivity constitutive relation:

```math
\\kappa=\\frac{1}{\\rho}+j\\omega\\varepsilon.
```

The common Coaxial Engine operator applies annular geometry and combines
adjacent dielectric layers radially in series.

# Arguments

- `material`: Static dielectric properties.
- `frequency`: Evaluation frequency \\[Hz\\].
- `temperature`: Operating temperature \\[°C\\].
- `values`: Formula assumptions.

# Returns

- Complex material admittivity ``1/\\rho+j\\omega\\varepsilon`` \\[S/m\\].

# References

A. Ametani, Y. Miyamoto, and N. Nagaoka, 2004, as reproduced in Ametani,
Ohno, and Nagaoka (2015), Eqs. 2.66–2.67.
"""
@inline function insulation_material(
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
