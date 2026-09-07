function routes(identifier::Val{:Knight2016})
    return (
        inner = FormulaMethod(identifier, internal_impedance, Val(:inner)),
        outer = FormulaMethod(identifier, internal_impedance, Val(:outer)),
        mutual = FormulaMethod(identifier, internal_impedance, Val(:mutual))
    )
end

assumptions(::Val{:Knight2016}) = (;)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Solid homogeneous circular cylinder of diameter ``d=2r``. |
| Calculated quantities | Continuous closed approximations for p.u.l. AC resistance and internal inductance of an isolated solid round conductor |
| Earth structure | None. |
| Model and approximation | TED-ML and PACAML are empirical modified-Lorentzian corrections of functions constrained to both exact asymptotes. Their reported errors refer to comparison with the report's Kelvin-function calculation, not independent measurements. |
| Main source | David W. Knight (version 2.08.1, 2016) |
| Citation key(s) | `:Knight2016` |
| Evidence status | Author-report page image verified |

**Expression.** The TED-ML resistance fit and PACAML inductance fit give
``Z=R_{dc}\\Xi+j\\omega\\mu\\Theta/(8\\pi)``.
The decimal coefficients are evaluated as decimal rational numbers in the
input numerical type. The zero-frequency result is the DC resistance.

**Reference.** [Knight2016](@cite), resistance fit on p. 31 and inductance
fit on p. 48 of version 2.08.1.
"""
description(::Formula{:Knight2016}) = "Knight solid-round internal impedance (2016)"

function (formula::Formula{:Knight2016})(
        r_in::T, r_ex::T, rho_c::T, mur_c::T, jω::Complex{T}
) where {T <: Real}
    iszero(r_in) || throw(DomainError(r_in, ":Knight2016 requires a solid conductor"))
    r_ex > zero(T) && rho_c > zero(T) && mur_c > zero(T) ||
        throw(DomainError((r_ex, rho_c, mur_c), "radius, resistivity, and permeability must be positive"))
    iszero(real(jω)) || throw(DomainError(jω, ":Knight2016 is a real-frequency fit"))
    μ = vacuum_permeability(r_ex) * mur_c
    state = (; r_ex, rho_c, μ, jω)
    return Functor{:Knight2016, typeof(formula.routes), typeof(state)}(formula.routes, state)
end

@inline (functor::Functor{:Knight2016})(::Val{:inner}) =
    functor.routes.inner(functor.state)
@inline (functor::Functor{:Knight2016})(::Val{:outer}) =
    functor.routes.outer(functor.state)
@inline (functor::Functor{:Knight2016})(::Val{:mutual}) =
    functor.routes.mutual(functor.state)

@inline internal_impedance(::Val{:Knight2016}, ::Val{:inner}, state) = zero(state.jω)
@inline internal_impedance(::Val{:Knight2016}, ::Val{:mutual}, state) = zero(state.jω)

function internal_impedance(::Val{:Knight2016}, ::Val{:outer}, state)
    T = typeof(state.r_ex)
    rdc = state.rho_c / ((one(T) * π) * state.r_ex^2)
    ω = abs(imag(state.jω))
    iszero(ω) && return complex(rdc, zero(T))
    c(n, d) = T(n) / T(d)
    t = state.r_ex * sqrt(ω * state.μ / (2 * state.rho_c))
    depth = -expm1(-t) / t
    z = c(62006, 100000) * t
    yr = c(189774, 1000000) /
         (1 + c(272481, 1000000) *
          (z^c(182938, 100000) - z^(-c(99457, 100000)))^2)^c(10941, 10000)
    resistance = rdc / (depth * (2 - depth) * (1 + yr))

    q = sqrt(T(2)) * t
    theta_inf = 2 / t * (
        1 + c(1209, 100000) / (q + 1) -
        c(63523, 100000) / (q^2 + 1) +
        c(16476, 100000) / (q^3 + 1)
    )
    power = c(15819, 10000)
    theta_da = theta_inf * (-expm1(-theta_inf^(-power)))^(inv(power))
    z = c(38691, 100000) * q
    yl = -c(198584, 1000000) /
         (1 + c(25741, 100000) *
          (z^c(12652, 10000) - z^(-c(39709, 100000)))^2)^c(262343, 100000)
    inductance = state.μ / (8 * (one(T) * π)) * theta_da * (1 - yl)
    return complex(resistance, imag(state.jω) * inductance)
end

:Knight2016
