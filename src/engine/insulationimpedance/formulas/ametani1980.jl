assumptions(::Val{:Ametani1980}) = (;)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Insulation impedance |
| Geometry | Concentric annuli ``r_2<r<r_3``, ``r_4<r<r_5``, and ``r_6<r<r_7``. |
| Calculated quantities | Per-unit-length magnetic series contributions of the core–sheath, sheath–armor, and armor–exterior insulation regions |
| Earth structure | Not applicable. |
| Model and approximation | Not an analytical approximation within the stated coaxial component model; no expansion or truncation is printed. |
| Main source | A. Ametani (1980) |
| Citation key(s) | `:Ametani1980` |
| Evidence status | PDF page images checked |

**Expression.**

```math
Z_{ins}=\\frac{j\\omega\\mu_0\\mu_r}{2\\pi}\\ln\\frac{b}{a}.
```

The term vanishes when the annular region has zero thickness. It is assembled
with the conductor surface impedances to form the cable series-impedance
matrix.

**Reference.** A. Ametani, “A General Formulation of Impedance and Admittance
of Cables,” *IEEE Transactions on Power Apparatus and Systems*, PAS-99(3),
902–910, 1980. DOI: 10.1109/TPAS.1980.319718.
"""
function description(::Formula{:Ametani1980})
    "Ametani coaxial-insulation magnetic impedance (1980)"
end

"""
$(TYPEDSIGNATURES)

Calculate the longitudinal magnetic impedance of one concentric insulation
region as used in Ametani's single-core cable assembly:

```math
z_{ab}=\\frac{s\\mu_0\\mu_i}{2\\pi}\\ln\\frac{b}{a}.
```

# Arguments

- `r_in`: Inner insulation radius ``a`` \\[m\\].
- `r_ex`: Outer insulation radius ``b`` \\[m\\].
- `mu_r`: Relative insulation permeability ``\\mu_i`` \\[dimensionless\\].
- `s`: Complex angular frequency ``s=j\\omega`` \\[rad/s\\].
- `values`: Formula assumptions.

# Returns

- Longitudinal insulation impedance ``z_{ab}`` \\[Ω/m\\].

# Notes

Implements Ametani (1980) as reproduced in Ametani, Ohno, and Nagaoka
(2015), Eqs. 2.6–2.13, and Ametani et al. (2021), Appendix A1.1.1.
"""
@inline function insulation_impedance(
        ::Val{:Ametani1980},
        r_in::T,
        r_ex::T,
        mu_r::T,
        s::Complex{T},
        values::NamedTuple
) where {T <: Real}
    if iszero(r_in) || r_in==r_ex
        return zero(Complex{T})
    end
    μ0 = one(r_in) * 4 * (one(r_in) * π) * (one(r_in) * 10)^(-7)
    return s * μ0 * mu_r / (2 * (one(r_in) * π)) * log1p((r_ex-r_in)/r_in)
end

:Ametani1980
