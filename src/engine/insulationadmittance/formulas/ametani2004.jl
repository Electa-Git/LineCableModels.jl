assumptions(::Val{:Ametani2004}) = (;)

"""
$(TYPEDSIGNATURES)

**Identification.** Complex-permittivity representation of a concentric
dielectric layer. Material conduction and displacement current enter the same
radial dielectric network.

**Expression.** For material resistivity ``\\rho``, real permittivity
``\\varepsilon'`` and polarization loss tangent ``\\tan\\delta_p``,

```math
\\varepsilon^\\star=\\varepsilon'(1-j\\tan\\delta_p)+\\frac{1}{j\\omega\\rho},\\qquad
\\kappa=\\frac{1}{\\rho}+\\omega\\varepsilon'\\tan\\delta_p+j\\omega\\varepsilon',
```

and an annular layer from ``a`` to ``b`` has

```math
Y=\\frac{2\\pi\\kappa}{\\ln(b/a)}.
```

The conductive complex-permittivity term and radial series combination follow
Ametani et al. (2004), Eqs. (14)–(15). Nonzero `material.tan_delta` supplies an
additional polarization-loss contribution; the paper does not prescribe this
parameter or a relaxation spectrum. It must exclude conduction already represented
by `rho`. Constant material inputs do not imply a fitted dispersion model.

**Reference.** A. Ametani, Y. Miyamoto, and N. Nagaoka, “Semiconducting Layer
Impedance and its Effect on Cable Wave-Propagation and Transient Characteristics,”
IEEE Transactions on Power Delivery, 19(4), 1523–1531, 2004,
doi:10.1109/TPWRD.2003.822502, Eqs. (14)–(15).
"""
function description(::Formula{:Ametani2004})
    "Ametani complex-permittivity dielectric admittance model (2004)"
end

"""
$(TYPEDSIGNATURES)

Evaluate Ametani's complex-permittivity constitutive relation:

```math
\\kappa=\\frac{1}{\\rho}+\\omega\\varepsilon'\\tan\\delta_p+j\\omega\\varepsilon'.
```

The common Coaxial Engine operator applies annular geometry and combines
adjacent dielectric layers radially in series.

# Arguments

- `material`: Static dielectric properties. `tan_delta` is the polarization
  contribution only, excluding conduction specified through `rho`.
- `frequency`: Evaluation frequency \\[Hz\\].
- `temperature`: Operating temperature \\[°C\\].
- `values`: Formula assumptions.

# Returns

- Complex material admittivity \\[S/m\\].

# References

A. Ametani, Y. Miyamoto, and N. Nagaoka (2004),
doi:10.1109/TPWRD.2003.822502, Eqs. (14)–(15). The optional polarization-loss
parameter is a material input, not an additional empirical law from that paper.
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
    displacement = complex(zero(T), ω) * ε₀ * material.eps_r
    return conductivity(material.rho) + imag(displacement) * material.tan_delta +
           displacement
end

:Ametani2004
