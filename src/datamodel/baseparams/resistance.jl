"""
$(TYPEDSIGNATURES)

Calculate the equivalent temperature coefficient of two parallel conductors:

```math
\\alpha_{eq}=\\frac{\\alpha_1R_2+\\alpha_2R_1}{R_1+R_2}.
```
"""
function equivalent_alpha(alpha1::Real, R1::Real, alpha2::Real, R2::Real)
    a1, r1, a2, r2 = promote(float(alpha1), float(R1), float(alpha2), float(R2))
    return (a1 * r2 + a2 * r1) / (r1 + r2)
end

"""
$(TYPEDSIGNATURES)

Calculate the parallel equivalent
``Z_{eq}=Z_1Z_2/(Z_1+Z_2)``. Open- and short-circuit edge cases are
handled explicitly.
"""
function parallel(Z1::Number, Z2::Number)
    z1, z2 = promote(Z1, Z2)
    isinf(z1) && return z2
    isinf(z2) && return z1
    iszero(z1) && iszero(z2) && return zero(z1)
    return z1 * z2 / (z1 + z2)
end

"""
$(TYPEDSIGNATURES)

Calculate reference-state strip resistance per unit length:

```math
R=\\frac{\\rho}{w t}.
```

# Arguments

- `thickness`: Strip thickness \\[m\\].
- `width`: Strip width \\[m\\].
- `rho`: Material resistivity at its reference temperature \\[Ω·m\\].

# Returns

- Strip resistance per unit length \\[Ω/m\\].
"""
function strip_resistance(
        thickness::Real,
        width::Real,
        rho::Real
)
    t, w, resistivity = promote(float(thickness), float(width), float(rho))
    return resistivity / (t * w)
end

"""
$(TYPEDSIGNATURES)

Calculate reference-state tubular resistance per unit length:

```math
R=\\frac{\\rho}{\\pi(r_{ex}^2-r_{in}^2)}.
```

# Arguments

- `r_in`: Inner conductor radius \\[m\\].
- `r_ex`: Outer conductor radius \\[m\\].
- `rho`: Material resistivity at its reference temperature \\[Ω·m\\].

# Returns

- Tubular resistance per unit length \\[Ω/m\\].
"""
function tubular_resistance(
        r_in::Real,
        r_ex::Real,
        rho::Real
)
    rin, rex, resistivity = promote(float(r_in), float(r_ex), float(rho))
    return resistivity / ((one(rin) * π) * (rex^2 - rin^2))
end

"""
$(TYPEDSIGNATURES)

Recover the equivalent resistivity from resistance and annular area:
``\\rho_{eq}=R\\pi(r_{ex}^2-r_{in}^2)``.
"""
function equivalent_rho(R::Real, r_ex::Real, r_in::Real)
    resistance, rex, rin = promote(float(R), float(r_ex), float(r_in))
    return resistance * (one(rin) * π) * (rex^2 - rin^2)
end
