function routes(::Val{:AmetaniFuse1992})
    (
        inner = ametani_fuse_zero,
        outer = ametani_fuse_outer,
        mutual = ametani_fuse_zero
    )
end

assumptions(::Val{:AmetaniFuse1992}) = (;)

description(::Formula{:AmetaniFuse1992}) = "Ametani-Fuse cross-section approximation (1992)"

"""
$(TYPEDSIGNATURES)

Construct the Ametani-Fuse approximation for the internal impedance of a
circular conductor interpreted through its cross-sectional area ``S`` and
outer perimeter ``\\ell``:

```math
Z_i\\approx R_{dc}\\sqrt{1+j\\omega\\mu_c
\\frac{S^2}{R_{dc}^2\\ell^2}}
=\\sqrt{Z_{dc}^2+Z_{hf}^2},
```

where

```math
R_{dc}=\\frac{\\rho_c}{S},
\\qquad
Z_{hf}=\\frac{\\sqrt{j\\omega\\mu_c\\rho_c}}{\\ell}.
```

For the admitted circular geometry, ``S=\\pi(b^2-a^2)`` and
``\\ell=2\\pi b``.

# Arguments

- `r_in`: Inner conductor radius ``a`` \\[m\\].
- `r_ex`: Outer conductor radius ``b`` \\[m\\].
- `rho_c`: Conductor resistivity ``\\rho_c`` \\[Ω·m\\].
- `mur_c`: Relative conductor permeability \\[dimensionless\\].
- `jω`: Complex angular frequency ``j\\omega`` \\[rad/s\\].

# Returns

- A formula functor whose outer interaction is ``Z_i`` \\[Ω/m\\]. Inner and
  mutual interactions are zero because the approximation defines only one
  longitudinal conductor surface.

# Notes

Implements Ametani and Fuse (1992) as reproduced in Ametani, Ohno, and
Nagaoka (2015), Eqs. 2.52 and 2.C.1–2.C.6. Sector and arbitrary polygonal
geometries are outside this implementation pass. The printed expanded form
in Eq. 2.C.6 is dimensionally inconsistent with its immediately preceding
definitions. Evaluation therefore uses the author-defined
``\\sqrt{Z_{dc}^2+Z_{hf}^2}`` identity and explicit ``Z_{hf}``, not that
expanded typo.
"""
function ametanifuse1992(
        formula::Formula{:AmetaniFuse1992},
        r_in::T,
        r_ex::T,
        rho_c::T,
        mur_c::T,
        jω::Complex{T}
) where {T <: Real}
    S = (one(T) * π) * (r_ex^2 - r_in^2)
    ell = 2 * (one(T) * π) * r_ex
    μ_c = vacuum_permeability(r_ex) * mur_c
    R_dc = rho_c / S
    state = (; r_in, r_ex, rho_c, mur_c, jω, μ_c, S, ell, R_dc)
    return Functor{:AmetaniFuse1992, typeof(formula.routes), typeof(state)}(
        formula.routes,
        state
    )
end

@inline function (formula::Formula{:AmetaniFuse1992})(r_in, r_ex, rho_c, mur_c, jω)
    return ametanifuse1992(formula, r_in, r_ex, rho_c, mur_c, jω)
end

@inline function (functor::Functor{:AmetaniFuse1992})(::Val{:inner})
    return functor.routes.inner(functor.state)
end

@inline function (functor::Functor{:AmetaniFuse1992})(::Val{:outer})
    return functor.routes.outer(functor.state)
end

@inline function (functor::Functor{:AmetaniFuse1992})(::Val{:mutual})
    return functor.routes.mutual(functor.state)
end

@inline ametani_fuse_zero(state) = zero(state.jω)

@inline function ametani_fuse_outer(state)
    Z_hf = sqrt(state.jω * state.μ_c * state.rho_c) / state.ell
    return sqrt(state.R_dc^2 + Z_hf^2)
end

:AmetaniFuse1992
