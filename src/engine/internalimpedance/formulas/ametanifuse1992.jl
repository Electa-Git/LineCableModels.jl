function routes(identifier::Val{:Ametani1992})
    (
        inner = FormulaMethod(identifier, internal_impedance, Val(:inner)),
        outer = FormulaMethod(identifier, internal_impedance, Val(:outer)),
        mutual = FormulaMethod(identifier, internal_impedance, Val(:mutual))
    )
end

assumptions(::Val{:Ametani1992}) = (;)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Isolated homogeneous conductor of arbitrary cross-section with area ``S`` and perimeter ``\\ell``; this includes solid and annular sectors. |
| Calculated quantities | Approximate frequency-dependent p.u.l. internal impedance and the radii of an impedance-equivalent circular annulus |
| Earth structure | Not applicable. |
| Model and approximation | ``Z_i=R_{dc}\\sqrt{1+j\\omega\\mu_cS/(R_{dc}\\ell^2)}`` interpolates between ``R_{dc}=\\rho_c/S`` and the high-frequency surface-layer limit ``\\sqrt{j\\omega\\mu_c\\rho_c}/\\ell``. The equivalent annulus has ``r_o=\\ell/(2\\pi)`` and ``r_i=\\sqrt{r_o^2-S/\\pi}``. The approximation retains skin effect through area and perimeter but does not resolve corner current crowding or proximity effect. |
| Main source | Akihiro Ametani and Ikuko Fuse (1992 English translation; 1991 Japanese original) |
| Citation key(s) | Primary English publication: `:Ametani1992`; Japanese original: `:Ametani1991`; later equation witness: `:Ametani2021` |
| Evidence status | Japanese original equations (1)–(6) and (14)–(15), sector-cable examples, and the later book's equation (A1.59) checked against page images; English publication metadata and abstract checked |

**Expression.** With cross-sectional area ``S`` and outer perimeter ``\\ell``,

```math
Z_i\\simeq R_{dc}\\sqrt{1+j\\omega\\mu_c
\\frac{S}{R_{dc}\\ell^2}}
=\\sqrt{Z_{dc}^2+Z_{hf}^2},\\qquad
R_{dc}=\\frac{\\rho_c}{S},\\quad
Z_{hf}=\\frac{\\sqrt{j\\omega\\mu_c\\rho_c}}{\\ell}.
```

The direct section call accepts physical `S` and `ell` explicitly;
the coaxial adapter retains these quantities and the original material.
Its scalar route applies to isolated homogeneous conductors, including
the supported sector shapes. It does not supply shell transfer terms,
proximity correction, or strand/composite homogenization.

The circular call uses ``S=\\pi(b^2-a^2)`` and
``\\ell=2\\pi b``. For a complete annulus this is the outer
current-carrying boundary; an open sector uses its full contour.

**Reference.** A. Ametani and I. Fuse, 1992, as reproduced in A. Ametani,
T. Ohno, and N. Nagaoka, *Cable System Transients: Theory, Modeling and
Simulation*, Wiley-IEEE Press, 2015, Eqs. 2.52 and 2.C.1–2.C.6.
"""
description(::Formula{:Ametani1992}) = "Ametani-Fuse cross-section approximation (1992)"

#=
$(TYPEDSIGNATURES)

Construct the Ametani-Fuse approximation for the internal impedance of a
circular conductor interpreted through its cross-sectional area ``S`` and
outer perimeter ``\\ell``:

```math
Z_i\\approx R_{dc}\\sqrt{1+j\\omega\\mu_c
\\frac{S}{R_{dc}\\ell^2}}
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
Nagaoka (2015), Eqs. 2.52 and 2.C.1–2.C.6. The explicit section call and
physical-section adapter retain the original area and perimeter. The printed expanded form
in Eq. 2.C.6 is dimensionally inconsistent with its immediately preceding
definitions. Evaluation therefore uses the author-defined
``\\sqrt{Z_{dc}^2+Z_{hf}^2}`` identity and explicit ``Z_{hf}``, not that
expanded typo.
=#
function (formula::Formula{:Ametani1992})(
        r_in::T,
        r_ex::T,
        rho_c::T,
        mur_c::T,
        jω::Complex{T}
) where {T <: Real}
    isfinite(r_in) && isfinite(r_ex) && 0<=r_in<r_ex ||
        throw(DomainError((r_in,r_ex),"invalid Ametani circular radii"))
    S = (one(T) * π) * (r_ex^2 - r_in^2)
    ell = 2 * (one(T) * π) * r_ex
    return formula(Val(:section),S,ell,rho_c,mur_c,jω)
end

function (formula::Formula{:Ametani1992})(::Val{:section},S::T,ell::T,
        rho_c::T,mur_c::T,jω::Complex{T}) where {T <: Real}
    all(isfinite,(S,ell,rho_c,mur_c,jω)) && min(S,ell,rho_c,mur_c)>0 ||
        throw(DomainError((S,ell,rho_c,mur_c,jω),
            "Ametani requires a finite positive area, perimeter, resistivity, and permeability"))
    r_ex=ell/(2*(one(T)*π))
    squared=r_ex^2-S/(one(T)*π)
    squared>=-64eps(T)*r_ex^2 || throw(DomainError((S,ell),
        "the area exceeds the isoperimetric bound for the supplied perimeter"))
    r_in=sqrt(max(zero(T),squared))
    μ_c = vacuum_permeability(S) * mur_c
    R_dc = rho_c / S
    state = (; r_in, r_ex, rho_c, mur_c, jω, μ_c, S, ell, R_dc)
    return Functor{:Ametani1992, typeof(formula.routes), typeof(state)}(
        formula.routes,
        state
    )
end

@inline function (functor::Functor{:Ametani1992})(::Val{:inner})
    return functor.routes.inner(functor.state)
end

@inline function (functor::Functor{:Ametani1992})(::Val{:outer})
    return functor.routes.outer(functor.state)
end

@inline function (functor::Functor{:Ametani1992})(::Val{:mutual})
    return functor.routes.mutual(functor.state)
end

@inline function internal_impedance(
        ::Val{:Ametani1992},
        ::Val{:inner},
        state
)
    return zero(state.jω)
end

@inline function internal_impedance(
        ::Val{:Ametani1992},
        ::Val{:mutual},
        state
)
    return zero(state.jω)
end

@inline function internal_impedance(
        ::Val{:Ametani1992},
        ::Val{:outer},
        state
)
    Z_hf = sqrt(state.jω * state.μ_c * state.rho_c) / state.ell
    return sqrt(state.R_dc^2 + Z_hf^2)
end

:Ametani1992
