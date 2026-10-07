
"""
$(TYPEDSIGNATURES)

**Identification.** Longitudinal magnetic impedance of one concentric
insulation region in Ametani's single-core cable formulation.

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
function description(::Type{<:Formula{:ametani1980}}; compact::Bool=false)
    compact ? "Ametani" : "Ametani coaxial-insulation magnetic impedance (1980)"
end

"""
$(TYPEDSIGNATURES)

Calculate the longitudinal magnetic impedance of one concentric insulation
region as used in Ametani's single-core cable assembly:

```math
z_{ab}=\\frac{s\\mu_0\\mu_i}{2\\pi}\\ln\\frac{b}{a}.
```

# Arguments

- `functor`: the Functor of the evaluation point. Its input holds:
  - `r_in`: inner insulation radius ``a`` \\[m\\].
  - `r_ex`: outer insulation radius ``b`` \\[m\\].
  - `mu_r`: relative insulation permeability ``\\mu_i`` \\[dimensionless\\].
  - `jω`: complex angular frequency ``j\\omega`` \\[rad/s\\].
  - `options`: the normalized formulation options of the formula.
- `workspace`: optional computation workspace supplying reusable numerical buffers.

# Returns

- Longitudinal insulation impedance ``z_{ab}`` \\[Ω/m\\].

# Notes

Implements Ametani (1980) as reproduced in Ametani, Ohno, and Nagaoka
(2015), Eqs. 2.6–2.13, and Ametani et al. (2021), Appendix A1.1.1.
"""
@inline function insulation_impedance(::Formula{:ametani1980}, functor, workspace)
    (; r_in, r_ex, mu_r, jω) = functor.input
    T = typeof(r_in)
    if isapprox(r_in, zero(T); atol = eps(T)) ||
       isapprox(r_in, r_ex; atol = eps(T))
        return zero(Complex{T})
    end
    μ0 = vacuum_permeability(typeof(r_in))
    return jω * μ0 * mu_r / (2 * (one(r_in) * π)) * log(r_ex / r_in)
end

formulation_options(::Expression{<:Formula{:ametani1980}, typeof(insulation_impedance)}) = FormulationOptions()

:ametani1980
