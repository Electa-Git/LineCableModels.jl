function routes(::Val{:Schelkunoff1934})
    (
        inner = schelkunoff_inner,
        outer = schelkunoff_outer,
        mutual = schelkunoff_mutual
    )
end

assumptions(::Val{:Schelkunoff1934}) = (;)

"""
$(TYPEDSIGNATURES)

**Identification.** Exact cylindrical surface impedances for a solid or
hollow round conductor.

**Expression.**

```math
\\begin{aligned}
Z_{is}&=\\frac{\\rho m}{2\\pi aD}
[I_0(ma)K_1(mb)+K_0(ma)I_1(mb)],\\\\
Z_{os}&=\\frac{\\rho m}{2\\pi bD}
[I_0(mb)K_1(ma)+K_0(mb)I_1(ma)],\\\\
Z_{ms}&=\\frac{\\rho m}{2\\pi abD},\\\\
D&=I_1(mb)K_1(ma)-K_1(mb)I_1(ma).
\\end{aligned}
```

For ``a=0``, ``Z_{int}=\\rho mI_0(mb)/(2\\pi bI_1(mb))``.

**Reference.** S. A. Schelkunoff, “The Electromagnetic Theory of Coaxial
Transmission Lines and Cylindrical Shields,” *Bell System Technical Journal*,
13, 532–579, 1934.
"""
function description(::Formula{:Schelkunoff1934})
    "Schelkunoff exact round-conductor surface impedances (1934)"
end

"""
$(TYPEDSIGNATURES)

Construct the exact Schelkunoff surface-impedance evaluator for one solid or
hollow circular conductor:

```math
Z_{is}=\\frac{\\rho m}{2\\pi aD}
\\left[I_0(ma)K_1(mb)+K_0(ma)I_1(mb)\\right],
\\qquad
Z_{ms}=\\frac{\\rho m}{2\\pi abD},
```

```math
Z_{os}=\\frac{\\rho m}{2\\pi bD}
\\left[I_0(mb)K_1(ma)+K_0(mb)I_1(ma)\\right],
```

where

```math
D=I_1(mb)K_1(ma)-K_1(mb)I_1(ma),
\\qquad
m=\\sqrt{j\\omega\\mu/\\rho}.
```

For ``a=0``, the outer term is evaluated from the solid-cylinder limit
``Z_{int}=\\rho m I_0(mb)/(2\\pi b I_1(mb))``.

# Arguments

- `r_in`: Inner conductor radius ``a`` \\[m\\].
- `r_ex`: Outer conductor radius ``b`` \\[m\\].
- `rho_c`: Conductor resistivity ``\\rho`` \\[Ω·m\\].
- `mur_c`: Relative conductor permeability \\[dimensionless\\].
- `jω`: Complex angular frequency ``j\\omega`` \\[rad/s\\].

# Returns

- A formula functor evaluating inner, outer, and mutual surface impedances
  \\[Ω/m\\].

# Notes

Implements Schelkunoff (1934) as reproduced in Ametani et al. (2021),
Appendix A1.4.1, Eqs. A1.50–A1.55.
"""
function schelkunoff1934(
        formula::Formula{:Schelkunoff1934},
        r_in::T,
        r_ex::T,
        rho_c::T,
        mur_c::T,
        jω::Complex{T}
) where {T <: Real}
    mu_c = vacuum_permeability(r_in) * mur_c
    sigma_c = conductivity(rho_c)
    m = sqrt(jω * mu_c * sigma_c)
    w_ex = m * r_ex
    state = (;
        r_in,
        r_ex,
        rho_c,
        mur_c,
        jω,
        mu_c,
        sigma_c,
        m,
        w_ex
    )
    return Functor{:Schelkunoff1934, typeof(formula.routes), typeof(state)}(
        formula.routes,
        state
    )
end

@inline function (formula::Formula{:Schelkunoff1934})(r_in, r_ex, rho_c, mur_c, jω)
    return schelkunoff1934(formula, r_in, r_ex, rho_c, mur_c, jω)
end

@inline function (functor::Functor{:Schelkunoff1934})(::Val{:inner})
    return functor.routes.inner(functor.state)
end

@inline function (functor::Functor{:Schelkunoff1934})(::Val{:outer})
    return functor.routes.outer(functor.state)
end

@inline function (functor::Functor{:Schelkunoff1934})(::Val{:mutual})
    return functor.routes.mutual(functor.state)
end

@inline function (functor::Functor{:Schelkunoff1934})(form::Symbol)
    Base.@nospecialize form
    return form === :inner ? functor(Val(:inner)) :
           form === :outer ? functor(Val(:outer)) :
           form === :mutual ? functor(Val(:mutual)) :
           throw(ArgumentError("unknown Schelkunoff1934 interaction: $form"))
end

@inline function schelkunoff_inner(state)
    T = typeof(state.r_in)
    if isapprox(state.r_in, zero(T); atol = eps(T))
        return zero(Complex{T})
    end

    w_in = state.m * state.r_in
    sc_in = exp(abs(real(w_in)) - state.w_ex)
    sc_ex = exp(abs(real(state.w_ex)) - w_in)
    sc = sc_in / sc_ex
    numerator = special_besselkx(0, w_in) * special_besselix(1, state.w_ex) +
                sc * special_besselix(0, w_in) * special_besselkx(1, state.w_ex)
    denominator = special_besselkx(1, w_in) * special_besselix(1, state.w_ex) -
                  sc * special_besselix(1, w_in) * special_besselkx(1, state.w_ex)
    return Complex{T}(
        (state.jω * state.mu_c / 2π) * (1 / w_in) * (numerator / denominator)
    )
end

@inline function schelkunoff_outer(state)
    T = typeof(state.r_in)
    if isapprox(state.r_in, zero(T); atol = eps(T))
        numerator = special_besselix(0, state.w_ex)
        denominator = special_besselix(1, state.w_ex)
    else
        w_in = state.m * state.r_in
        sc_in = exp(abs(real(w_in)) - state.w_ex)
        sc_ex = exp(abs(real(state.w_ex)) - w_in)
        sc = sc_in / sc_ex
        numerator = special_besselix(0, state.w_ex) * special_besselkx(1, w_in) +
                    sc * special_besselkx(0, state.w_ex) * special_besselix(1, w_in)
        denominator = special_besselix(1, state.w_ex) * special_besselkx(1, w_in) -
                      sc * special_besselkx(1, state.w_ex) * special_besselix(1, w_in)
    end
    return Complex{T}(
        (state.jω * state.mu_c / 2π) * (1 / state.w_ex) *
        (numerator / denominator)
    )
end

@inline function schelkunoff_mutual(state)
    T = typeof(state.r_in)
    if isapprox(state.r_in, zero(T); atol = eps(T))
        return zero(Complex{T})
    end

    w_in = state.m * state.r_in
    sc_in = exp(abs(real(w_in)) - state.w_ex)
    sc_ex = exp(abs(real(state.w_ex)) - w_in)
    sc = sc_in / sc_ex
    numerator = one(sc_ex) / sc_ex
    denominator = special_besselix(1, state.w_ex) * special_besselkx(1, w_in) -
                  sc * special_besselix(1, w_in) * special_besselkx(1, state.w_ex)
    return Complex{T}(
        (1 / (2π * state.r_in * state.r_ex * state.sigma_c)) *
        (numerator / denominator)
    )
end

:Schelkunoff1934
