function routes(identifier::Val{:WedepohlWilcox1973})
    (
        inner = FormulaMethod(identifier, internal_impedance, Val(:inner)),
        outer = FormulaMethod(identifier, internal_impedance, Val(:outer)),
        mutual = FormulaMethod(identifier, internal_impedance, Val(:mutual))
    )
end

assumptions(::Val{:WedepohlWilcox1973}) = (;)

"""
$(TYPEDSIGNATURES)

**Identification.** Hollow-shell approximation retaining inner, outer, and
transfer surface impedances.

**Expression.**

```math
\\begin{aligned}
Z_{iw}&=\\frac{\\rho m}{2\\pi a}\\coth[m(b-a)]
-\\frac{\\rho}{2\\pi b(a+b)},\\\\
Z_{ow}&=\\frac{\\rho m}{2\\pi b}\\coth[m(b-a)]
+\\frac{\\rho}{2\\pi b(a+b)},\\\\
Z_{mw}&=\\frac{\\rho m}{\\pi(a+b)}\\operatorname{csch}[m(b-a)].
\\end{aligned}
```

**Reference.** L. M. Wedepohl and D. J. Wilcox, 1973, as reproduced in
Ametani et al., *Electromagnetic Transients in Large HV Cable Networks*, IET,
2021, Appendix A1.4.2, Eqs. A1.56–A1.58.
"""
function description(::Formula{:WedepohlWilcox1973})
    "Wedepohl-Wilcox approximate shell impedances (1973)"
end

#=
$(TYPEDSIGNATURES)

Construct the Wedepohl-Wilcox hollow-shell approximation:

```math
Z_{iw}=\\frac{\\rho m}{2\\pi a}\\coth[m(b-a)]
-\\frac{\\rho}{2\\pi b(a+b)},
```

```math
Z_{mw}=\\frac{\\rho m}{\\pi(a+b)}\\operatorname{csch}[m(b-a)],
\\qquad
Z_{ow}=\\frac{\\rho m}{2\\pi b}\\coth[m(b-a)]
+\\frac{\\rho}{2\\pi b(a+b)},
```

where ``m=\\sqrt{j\\omega\\mu/\\rho}``.

# Arguments

- `r_in`: Shell inner radius ``a`` \\[m\\].
- `r_ex`: Shell outer radius ``b`` \\[m\\].
- `rho_c`: Shell resistivity ``\\rho`` \\[Ω·m\\].
- `mur_c`: Relative shell permeability \\[dimensionless\\].
- `jω`: Complex angular frequency ``j\\omega`` \\[rad/s\\].

# Returns

- A formula functor evaluating inner, outer, and mutual surface impedances
  \\[Ω/m\\].

# Notes

The approximation applies only to a hollow circular shell. It implements
Wedepohl and Wilcox (1973) as reproduced in Ametani et al. (2021), Appendix
A1.4.2, Eqs. A1.56–A1.58.
=#
function (formula::Formula{:WedepohlWilcox1973})(
        r_in::T,
        r_ex::T,
        rho_c::T,
        mur_c::T,
        jω::Complex{T}
) where {T <: Real}
    r_in > zero(T) || throw(DomainError(
        r_in,
        ":WedepohlWilcox1973 requires a hollow conductor"
    ))
    μ = vacuum_permeability(r_in) * mur_c
    m = sqrt(jω * μ / rho_c)
    state = (; r_in, r_ex, rho_c, mur_c, jω, μ, m)
    return Functor{
        :WedepohlWilcox1973,
        typeof(formula.routes),
        typeof(state)
    }(formula.routes, state)
end

@inline function (functor::Functor{:WedepohlWilcox1973})(::Val{:inner})
    return functor.routes.inner(functor.state)
end

@inline function (functor::Functor{:WedepohlWilcox1973})(::Val{:outer})
    return functor.routes.outer(functor.state)
end

@inline function (functor::Functor{:WedepohlWilcox1973})(::Val{:mutual})
    return functor.routes.mutual(functor.state)
end

@inline function internal_impedance(
        ::Val{:WedepohlWilcox1973},
        ::Val{:inner},
        state
)
    a = state.r_in
    b = state.r_ex
    x = state.m * (b - a)
    return state.rho_c * state.m / (2π * a) * coth(x) -
           state.rho_c / (2π * b * (a + b))
end

@inline function internal_impedance(
        ::Val{:WedepohlWilcox1973},
        ::Val{:outer},
        state
)
    a = state.r_in
    b = state.r_ex
    x = state.m * (b - a)
    return state.rho_c * state.m / (2π * b) * coth(x) +
           state.rho_c / (2π * b * (a + b))
end

@inline function internal_impedance(
        ::Val{:WedepohlWilcox1973},
        ::Val{:mutual},
        state
)
    a = state.r_in
    b = state.r_ex
    x = state.m * (b - a)
    return state.rho_c * state.m / (π * (a + b) * sinh(x))
end

:WedepohlWilcox1973
