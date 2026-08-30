assumptions(::Val{:Ametani1980}) = (;)

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
@inline function ametani1980(
        r_in::T,
        r_ex::T,
        mu_r::T,
        s::Complex{T},
        values::NamedTuple
) where {T <: Real}
    if isapprox(r_in, zero(T); atol = eps(T)) ||
       isapprox(r_in, r_ex; atol = eps(T))
        return zero(Complex{T})
    end
    μ0 = one(r_in) * 4 * (one(r_in) * π) * (one(r_in) * 10)^(-7)
    return s * μ0 * mu_r / (2 * (one(r_in) * π)) * log(r_ex / r_in)
end

@inline function (::Functor{:Ametani1980})(
        r_in::T,
        r_ex::T,
        mu_r::T,
        s::Complex{T},
        values::NamedTuple
) where {T <: Real}
    return ametani1980(r_in, r_ex, mu_r, s, values)
end

:Ametani1980
