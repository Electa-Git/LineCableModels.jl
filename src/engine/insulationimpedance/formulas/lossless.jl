assumptions(::Val{:Lossless}) = (;)

"""
$(TYPEDSIGNATURES)

**Identification.** Ideal-dielectric limit of the concentric insulation
series impedance. Dielectric conductivity is zero, but magnetic energy stored
in the annulus remains.

**Expression.**

```math
Z_{ins}=j\\omega L_{ins},\\qquad
L_{ins}=\\frac{\\mu_0\\mu_r}{2\\pi}\\ln\\frac{b}{a}.
```

**Reference.** A. Ametani, T. Ohno, and N. Nagaoka, *Cable System
Transients: Theory, Modeling and Simulation*, Wiley-IEEE Press, 2015,
Sec. 2.2.
"""
description(::Formula{:Lossless}) = "Lossless insulation (ideal dielectric)"

"""
$(TYPEDSIGNATURES)

Calculate the longitudinal series impedance of one concentric ideal dielectric
layer:

```math
z(s) = \\frac{s\\mu_0\\mu_r}{2\\pi}\\ln\\left(\\frac{r_o}{r_i}\\right).
```

# Arguments

- `r_in`: Inner layer radius \\[m\\].
- `r_ex`: Outer layer radius \\[m\\].
- `mu_r`: Relative magnetic permeability \\[dimensionless\\].
- `s`: Complex angular frequency \\[rad/s\\].
- `values`: Formula assumptions.

# Returns

- Complex series impedance per unit length \\[Ω/m\\].
"""
@inline function lossless(
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

    mu0 = one(r_in) * 4 * (one(r_in) * π) * (one(r_in) * 10)^(-7)
    permeability = mu0 * mu_r
    return Complex{T}(
        s * permeability * log(r_ex / r_in) / (2 * (one(r_in) * π))
    )
end

@inline function (::Functor{:Lossless})(
        r_in::T,
        r_ex::T,
        mu_r::T,
        s::Complex{T},
        values::NamedTuple
) where {T <: Real}
    return lossless(r_in, r_ex, mu_r, s, values)
end

:Lossless
