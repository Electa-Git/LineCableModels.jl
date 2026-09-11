
"""
$(TYPEDSIGNATURES)

**Identification.** Complex-permittivity representation of a concentric
semiconducting screen. Its conduction and displacement currents enter the same
radial dielectric network as the adjacent insulation layers.

**Expression.** For screen resistivity ``\\rho_s``, real permittivity
``\\varepsilon'_s`` and polarization loss tangent ``\\tan\\delta_p``,

```math
\\varepsilon_s^\\star=\\varepsilon'_s(1-j\\tan\\delta_p)+\\frac{1}{j\\omega\\rho_s},\\qquad
\\kappa_s=\\frac{1}{\\rho_s}+\\omega\\varepsilon'_s\\tan\\delta_p+j\\omega\\varepsilon'_s,
```

and an annular screen from ``a`` to ``b`` has

```math
Y_s=\\frac{2\\pi\\kappa_s}{\\ln(b/a)}.
```

With zero `tan_delta`, this is Ametani et al. (2004), Eq. (14); the total
radial admittance follows Eq. (15). Nonzero `tan_delta` supplies an additional
material polarization loss, excluding conduction already represented by `rho`.
The paper does not prescribe that input or a dielectric relaxation spectrum.

**Reference.** A. Ametani, Y. Miyamoto, and N. Nagaoka, “Semiconducting Layer
Impedance and its Effect on Cable Wave-Propagation and Transient Characteristics,”
IEEE Transactions on Power Delivery, 19(4), 1523–1531, 2004,
doi:10.1109/TPWRD.2003.822502, Eqs. (14)–(15).
"""
function description(::Formula{:Ametani2004})
    "Ametani semiconducting-screen admittance model (2004)"
end

"""
$(TYPEDSIGNATURES)

Retain semiconducting-screen conductivity and permittivity in Ametani's
complex-permittivity representation:

```math
\\varepsilon_s^\\star=\\varepsilon'_s(1-j\\tan\\delta_p)+\\frac{1}{j\\omega\\rho_s},
\\qquad
\\kappa_s=\\frac{1}{\\rho_s}+\\omega\\varepsilon'_s\\tan\\delta_p+j\\omega\\varepsilon'_s.
```

The common Coaxial Engine operator applies the annular geometry and combines
the semiconducting screen with adjacent dielectric layers radially in series.

# Arguments

- `material`: Static semiconducting-screen properties; `tan_delta` contains
  polarization loss only and excludes conduction specified by `rho`.
- `frequency`: Evaluation frequency \\[Hz\\].
- `temperature`: Operating temperature \\[°C\\].
- `values`: Explicit physical/model parameters.
- `options`: Normalized numerical sections for this contribution.
- `workspace`: Optional execution resources.

# Returns

- Complex screen admittivity \\[S/m\\].

# References

A. Ametani, Y. Miyamoto, and N. Nagaoka (2004),
doi:10.1109/TPWRD.2003.822502, Eqs. (14)–(15). The optional polarization-loss
parameter is a material input, not an additional empirical law from that paper.
"""
@inline function semicon_material(
        ::Val{:Ametani2004},
        material::Material{T},
        frequency::T,
        temperature::T,
        values::NamedTuple, options::NamedTuple, workspace
) where {T <: Real}
    ε₀ = one(T) * 88541878128 * (one(T) * 10)^(-22)
    ω = 2 * (one(T) * π) * frequency
    displacement = complex(zero(T), ω) * ε₀ * material.eps_r
    return conductivity(material.rho) + imag(displacement) * material.tan_delta +
           displacement
end

computation_options(::FormulaMethod{:Ametani2004, typeof(semicon_material)}) = (;)

:Ametani2004
