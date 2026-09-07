function routes(identifier::Val{:BrandaoFaria2011})
    return (
        inner = FormulaMethod(identifier, internal_impedance, Val(:inner)),
        outer = FormulaMethod(identifier, internal_impedance, Val(:outer)),
        mutual = FormulaMethod(identifier, internal_impedance, Val(:mutual))
    )
end

assumptions(::Val{:BrandaoFaria2011}) = (exponent = 0,)

function Formula(::Val{:BrandaoFaria2011}; exponent::Real = 0, kwargs...)
    isfinite(exponent) || throw(DomainError(exponent, "the radial exponent must be finite"))
    defaults = routes(Val(:BrandaoFaria2011))
    overrides = (; kwargs...)
    unknown = setdiff(keys(overrides), keys(defaults))
    isempty(unknown) || throw(ArgumentError("unknown internal-impedance routes: $unknown"))
    return Formula(Val(:BrandaoFaria2011), merge(defaults, overrides), (; exponent))
end

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Circular tube with inner radius ``r_1`` and outer radius ``r_2``; the solid limit is ``r_1\\to0``; no insulation layer. |
| Calculated quantities | Frequency-dependent p.u.l. internal impedance of a radially inhomogeneous tubular conductor, plus its solid-cylinder limit |
| Earth structure | Not applicable. |
| Model and approximation | Not an analytical approximation after imposing the very-good-conductor diffusion model and the specific constitutive power laws. The Euler–Cauchy reduction is exact for (11)–(12); arbitrary radial profiles require numerical solution. |
| Main source | José António Brandão Faria (2011) |
| Citation key(s) | `:BrandaoFaria2011` |
| Evidence status | Original publication page images checked |

**Expression.** Equation (19) is evaluated with scaled exponentials.
The input resistivity and relative permeability are their values at the
outer radius. The keyword `exponent` specifies ``p`` in
``\\mu(r)=\\mu_2(r/r_2)^p`` and
``\\sigma(r)=\\sigma_2(r_2/r)^{2+p}``; its default is zero.
These profiles do not describe a homogeneous conductor, including when
``p=0``.

Only the outer impedance with zero current in the hollow interior is
provided. Inner and transfer surface requests are rejected for a tube.

**Reference.** [BrandaoFaria2011](@cite), equations (11)–(13) and (19).
"""
description(::Formula{:BrandaoFaria2011}) =
    "Brandão Faria power-law radial conductor impedance (2011)"

function (formula::Formula{:BrandaoFaria2011})(
        r_in::T, r_ex::T, rho_c::T, mur_c::T, jω::Complex{T}
) where {T <: Real}
    zero(T) <= r_in < r_ex && rho_c > zero(T) && mur_c > zero(T) ||
        throw(DomainError((r_in, r_ex, rho_c, mur_c),
            "require 0 ≤ inner radius < outer radius and positive material values"))
    p = T(formula.assumptions.exponent)
    μ = vacuum_permeability(r_ex) * mur_c
    state = (; r_in, r_ex, rho_c, μ, jω, p)
    return Functor{:BrandaoFaria2011, typeof(formula.routes), typeof(state)}(
        formula.routes, state
    )
end

@inline (functor::Functor{:BrandaoFaria2011})(::Val{:inner}) =
    functor.routes.inner(functor.state)
@inline (functor::Functor{:BrandaoFaria2011})(::Val{:outer}) =
    functor.routes.outer(functor.state)
@inline (functor::Functor{:BrandaoFaria2011})(::Val{:mutual}) =
    functor.routes.mutual(functor.state)

function internal_impedance(
        ::Val{:BrandaoFaria2011}, ::Union{Val{:inner}, Val{:mutual}}, state
)
    iszero(state.r_in) && return zero(state.jω)
    throw(ArgumentError(":BrandaoFaria2011 assumes zero current in the hollow interior"))
end

function internal_impedance(::Val{:BrandaoFaria2011}, ::Val{:outer}, state)
    T = typeof(state.r_ex)
    p = state.p
    factor = state.rho_c / (2 * (one(T) * π) * state.r_ex^2)
    if iszero(state.jω)
        if iszero(state.r_in)
            return complex(max(-p, zero(T)) * factor, zero(T))
        end
        ell = log(state.r_ex / state.r_in)
        dc = iszero(p) ? factor / ell : factor * p / expm1(p * ell)
        return complex(dc, zero(T))
    end
    frequency_term = state.jω * state.μ * state.r_ex^2 / state.rho_c
    root = sqrt((p / 2)^2 + frequency_term)
    # Avoid cancellation in root - p/2 for a positive exponent at low frequency.
    difference = p >= zero(T) ? frequency_term / (root + p / 2) : root - p / 2
    iszero(state.r_in) && return factor * difference
    ell = log(state.r_ex / state.r_in)
    attenuation = -2 * root * ell
    return factor * (
        difference + 2 * root * exp(attenuation) / (-expm1(attenuation))
    )
end

:BrandaoFaria2011
