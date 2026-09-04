function routes(identifier::Val{:Ametani2004})
    (
        inner = FormulaMethod(identifier, internal_impedance, Val(:inner)),
        outer = FormulaMethod(identifier, internal_impedance, Val(:outer)),
        mutual = FormulaMethod(identifier, internal_impedance, Val(:mutual))
    )
end

assumptions(::Val{:Ametani2004}) = (;)

"""
$(TYPEDSIGNATURES)

**Identification.** Bonded two-material cylindrical conductor, including a
conductive semiconducting layer.

**Expression.** If the two materials have Schelkunoff surface terms
``Z_{1i},Z_{1o},Z_{1m}`` and ``Z_{2i},Z_{2o},Z_{2m}``,

```math
Z_{out}=Z_{2o}-\\frac{Z_{2m}^2}{Z_{1o}+Z_{2i}},\\qquad
Z_{in}=Z_{1i}-\\frac{Z_{1m}^2}{Z_{1o}+Z_{2i}},\\qquad
Z_m=\\frac{Z_{1m}Z_{2m}}{Z_{1o}+Z_{2i}}.
```

**Reference.** A. Ametani, H. Xue, T. Ohno, and H. Khalilnezhad,
*Electromagnetic Transients in Large HV Cable Networks: Modeling and
Calculations*, IET, 2021, Appendix A1.4.
"""
function description(::Formula{:Ametani2004})
    "Ametani et al. bonded two-layer conductor impedance (2004)"
end

#=
Construct the Ametani et al. two-layer conductor/semiconducting-layer
surface impedance:

```math
Z_{out}=Z_{20}-\frac{Z_{2m}^2}{Z_{10}+Z_{2i}},\qquad
Z_{in}=Z_{1i}-\frac{Z_{1m}^2}{Z_{10}+Z_{2i}},\qquad
Z_m=\frac{Z_{1m}Z_{2m}}{Z_{10}+Z_{2i}}.
```

The first cylindrical material occupies ``a\le r\le b`` and the second
occupies ``b'\le r\le c``. Every ``Z_{ki}``, ``Z_{ko}``, and ``Z_{km}`` is
the exact Schelkunoff surface term of its material. The corpus also gives

```math
Z_{out}=\frac{Z_{11}Z_{22}-Z_{12}^2}{Z_{11}+Z_{22}-2Z_{12}},
```

with ``Z_{11}=z_{10}+z_{12}+z_{2i}+z_{20}-2z_{2m}``, 
``Z_{12}=z_{20}-z_{2m}``, and ``Z_{22}=z_{20}``.

# Arguments

- `a`, `b`: Inner and outer radii of material 1 [m].
- `b_prime`, `c`: Inner and outer radii of material 2 [m].
- `rho_1`, `rho_2`: Material resistivities [Ω·m].
- `mu_r1`, `mu_r2`: Relative permeabilities [dimensionless].
- `jω`: Complex angular frequency [rad/s].

# Notes

This aggregate requires both conductive material records at once. The ordinary
five-argument single-layer call is rejected rather than silently applying the
aggregate to each flattened conductor independently.
=#
function (formula::Formula{:Ametani2004})(
        a::T,
        b::T,
        b_prime::T,
        c::T,
        rho_1::T,
        rho_2::T,
        mu_r1::T,
        mu_r2::T,
        jω::Complex{T}
) where {T <: Real}
    first = surface_impedances(Val(:Schelkunoff1934), a, b, rho_1, mu_r1, jω)
    second = surface_impedances(
        Val(:Schelkunoff1934), b_prime, c, rho_2, mu_r2, jω
    )
    denominator = first.outer + second.inner
    state = (;
        jω,
        inner = first.inner - first.mutual^2 / denominator,
        outer = second.outer - second.mutual^2 / denominator,
        mutual = first.mutual * second.mutual / denominator
    )
    return Functor{:Ametani2004, typeof(formula.routes), typeof(state)}(
        formula.routes,
        state
    )
end

function (formula::Formula{:Ametani2004})(r_in, r_ex, rho, mu_r, jω)
    throw(ArgumentError(
        ":Ametani2004 requires both bonded conductive layers in one aggregate call"
    ))
end

@inline function internal_impedance(
        ::Val{:Ametani2004},
        ::Val{:inner},
        state
)
    return state.inner
end

@inline function internal_impedance(
        ::Val{:Ametani2004},
        ::Val{:outer},
        state
)
    return state.outer
end

@inline function internal_impedance(
        ::Val{:Ametani2004},
        ::Val{:mutual},
        state
)
    return state.mutual
end

@inline function (functor::Functor{:Ametani2004})(::Val{:inner})
    return functor.routes.inner(functor.state)
end

@inline function (functor::Functor{:Ametani2004})(::Val{:outer})
    return functor.routes.outer(functor.state)
end

@inline function (functor::Functor{:Ametani2004})(::Val{:mutual})
    return functor.routes.mutual(functor.state)
end

:Ametani2004
