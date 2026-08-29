assumptions(::Val{:Lossless}) = (;)

description(::Formula{:Lossless}) = "Lossless insulation (ideal dielectric)"

"""
$(TYPEDSIGNATURES)

Calculate the potential coefficient of one concentric ideal dielectric layer:

```math
p = \\frac{\\ln(r_o/r_i)}{2\\pi\\varepsilon_0\\varepsilon_r}.
```

# Arguments

- `r_in`: Inner layer radius \\[m\\].
- `r_ex`: Outer layer radius \\[m\\].
- `rho`: Dielectric resistivity \\[Ω·m\\]; ignored by the lossless formula.
- `eps_r`: Relative dielectric permittivity \\[dimensionless\\].
- `s`: Complex angular frequency \\[rad/s\\]; retained for the uniform route
  contract.
- `values`: Formula assumptions.

# Returns

- Potential coefficient per unit length \\[m/F\\].
"""
@inline function lossless(
        r_in::T,
        r_ex::T,
        rho::T,
        eps_r::T,
        s::Complex{T},
        values::NamedTuple
) where {T <: Real}
    if isapprox(r_in, zero(T); atol = eps(T)) ||
       isapprox(r_in, r_ex; atol = eps(T))
        return zero(Complex{T})
    end

    eps0 = one(r_in) * 88541878128 * (one(r_in) * 10)^(-22)
    permittivity = eps0 * eps_r
    return Complex{T}(log(r_ex / r_in) / (2 * (one(r_in) * π) * permittivity))
end

@inline function (::Functor{:Lossless})(
        r_in::T,
        r_ex::T,
        rho::T,
        eps_r::T,
        s::Complex{T},
        values::NamedTuple
) where {T <: Real}
    return lossless(r_in, r_ex, rho, eps_r, s, values)
end

:Lossless
