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

Return the linear resistivity correction factor
``k(T)=1+\\alpha(T-T_0)``.

# Notes

The linear material model is accepted only for temperature differences below
150 °C.
"""
function temperature_factor(alpha::Real, temperature::Real, reference::Real)
    a, value, base = promote(float(alpha), float(temperature), float(reference))
    delta = abs(value - base)
    delta < oftype(delta, 150) || throw(DomainError(
        delta,
        "the linear resistivity model requires |T - T0| < 150 °C"
    ))
    return one(a) + a * (value - base)
end

function temperature_factor(alpha::Real, temperature::Real)
    temperature_factor(alpha, temperature, oftype(float(temperature), 20))
end

"""
$(TYPEDSIGNATURES)

Calculate strip resistance per unit length,
``R=\\rho k(T)/(w t)`` \\[Ω/m\\].
"""
function strip_resistance(
        thickness::Real,
        width::Real,
        rho::Real,
        alpha::Real,
        reference::Real,
        temperature::Real
)
    t, w, resistivity, a, base, value = promote(
        float(thickness), float(width), float(rho), float(alpha),
        float(reference), float(temperature)
    )
    return temperature_factor(a, value, base) * resistivity / (t * w)
end

"""
$(TYPEDSIGNATURES)

Calculate tubular resistance per unit length,
``R=\\rho k(T)/(\\pi(r_{ex}^2-r_{in}^2))`` \\[Ω/m\\].
"""
function tubular_resistance(
        r_in::Real,
        r_ex::Real,
        rho::Real,
        alpha::Real,
        reference::Real,
        temperature::Real
)
    rin, rex, resistivity, a, base, value = promote(
        float(r_in), float(r_ex), float(rho), float(alpha),
        float(reference), float(temperature)
    )
    return temperature_factor(a, value, base) * resistivity /
           ((one(rin) * π) * (rex^2 - rin^2))
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
