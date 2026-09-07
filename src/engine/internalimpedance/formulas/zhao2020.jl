function routes(identifier::Val{:Zhao2020})
    (
        inner = FormulaMethod(identifier, internal_impedance, Val(:inner)),
        outer = FormulaMethod(identifier, internal_impedance, Val(:outer)),
        mutual = FormulaMethod(identifier, internal_impedance, Val(:mutual))
    )
end

assumptions(::Val{:Zhao2020}) = (;)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Homogeneous circular conducting annulus with 0<a<b. |
| Calculated quantities | Leading-asymptotic inner, outer, and transfer impedances |
| Earth structure | Not applicable. |
| Model and approximation | Large-argument Bessel approximation with geometric-mean radius in the transfer term. |
| Main source | H. Zhao, A. Ametani, A. Gole, and J. De Silva (2020) |
| Citation key(s) | Attributed source: `:Zhao2020`; equation witness: `:Ametani2021` |
| Evidence status | IET book page images checked, equations (A1.61)–(A1.65) and reference [90]. |

**Description.** High-frequency asymptotic hollow-shell impedances.

**Expression.**

```math
Z_{ia}=\\frac{\\rho m}{2\\pi a}\\coth[m(b-a)],\\qquad
Z_{oa}=\\frac{\\rho m}{2\\pi b}\\coth[m(b-a)],\\qquad
Z_{ma}=\\frac{\\rho m}{2\\pi\\sqrt{ab}}\\mathrm{csch}[m(b-a)].
```

**Reference.** Zhao et al., 2020, as reproduced in Ametani et al.,
*Electromagnetic Transients in Large HV Cable Networks*, IET, 2021,
Appendix A1.4.4.2, Eqs. A1.61–A1.65.
"""
function description(::Formula{:Zhao2020})
    "Zhao et al. high-frequency asymptotic shell impedances (2020)"
end

#=
$(TYPEDSIGNATURES)

Construct the high-frequency asymptotic shell evaluator:

```math
Z_{ia}=\\frac{\\rho m}{2\\pi a}\\coth[m(b-a)],
\\qquad
Z_{ma}=\\frac{\\rho m}{2\\pi\\sqrt{ab}}\\mathrm{csch}[m(b-a)],
\\qquad
Z_{oa}=\\frac{\\rho m}{2\\pi b}\\coth[m(b-a)],
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

Implements Zhao et al. (2020) as reproduced in Ametani et al. (2021),
Appendix A1.4.4.2, Eqs. A1.61–A1.65. The formula is an asymptotic
high-frequency approximation for hollow circular shells.
=#
function (formula::Formula{:Zhao2020})(
        r_in::T,
        r_ex::T,
        rho_c::T,
        mur_c::T,
        jω::Complex{T}
) where {T <: Real}
    r_in > zero(T) || throw(DomainError(
        r_in,
        ":Zhao2020 requires a hollow conductor"
    ))
    μ = vacuum_permeability(r_in) * mur_c
    m = sqrt(jω * μ / rho_c)
    state = (; r_in, r_ex, rho_c, mur_c, jω, μ, m)
    return Functor{:Zhao2020, typeof(formula.routes), typeof(state)}(
        formula.routes,
        state
    )
end

@inline function (functor::Functor{:Zhao2020})(::Val{:inner})
    return functor.routes.inner(functor.state)
end

@inline function (functor::Functor{:Zhao2020})(::Val{:outer})
    return functor.routes.outer(functor.state)
end

@inline function (functor::Functor{:Zhao2020})(::Val{:mutual})
    return functor.routes.mutual(functor.state)
end

@inline function internal_impedance(
        ::Val{:Zhao2020},
        ::Val{:inner},
        state
)
    x = state.m * (state.r_ex - state.r_in)
    return state.rho_c * _shell_factors(state.m,state.r_ex-state.r_in).coth / (2 * (one(state.r_in)*π) * state.r_in)
end

@inline function internal_impedance(
        ::Val{:Zhao2020},
        ::Val{:outer},
        state
)
    x = state.m * (state.r_ex - state.r_in)
    return state.rho_c * _shell_factors(state.m,state.r_ex-state.r_in).coth / (2 * (one(state.r_in)*π) * state.r_ex)
end

@inline function internal_impedance(
        ::Val{:Zhao2020},
        ::Val{:mutual},
        state
)
    x = state.m * (state.r_ex - state.r_in)
    return state.rho_c * _shell_factors(state.m,state.r_ex-state.r_in).csch /
           (2 * (one(state.r_in)*π) * sqrt(state.r_in * state.r_ex))
end

:Zhao2020
