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

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Medium 1 occupies ``a\\le r\\le b`` and medium 2 occupies ``b'\\le r\\le c``; for bonded media ``b=b'``. The first medium may be solid (``a=0``) or hollow. Both ends are short-circuited for the outer-terminal reduction. |
| Calculated quantities | Inner-surface, outer-surface, and transfer impedance of two electrically bonded concentric conductive media; stable assembly from single-layer surface terms and direct Maxwell/Bessel representation |
| Earth structure | Not applicable. |
| Model and approximation | Appendix (A.11) is the direct Maxwell/Bessel representation with ``\\Gamma`` and intrinsic admittance retained. The main-text construction (9)–(10) is algebraically equivalent under the paper's stated good-conductor and penetration-depth assumptions: ``j\\omega\\mu\\sigma\\gg\\omega^2\\mu\\varepsilon+\\Gamma^2`` and no radial current between media. It is a model reduction, not merely a notation change. |
| Main source | A. Ametani, Y. Miyamoto, and N. Nagaoka (2004); appendix says the direct Maxwell derivation summarizes N. Amekawa's 2001 thesis |
| Citation key(s) | `:Ametani2004` |
| Evidence status | PDF page images checked |

**Numerical scope.** The aggregate reduces (5)–(8) for bonded
physical materials. Repeated application uses the same binary
interface elimination. The five-argument call is the single-material
limit; the coaxial assembler supplies the complete physical profile.
The longitudinal/displacement-dependent appendix is not substituted
for this good-conductor reduction.
The survey retains the misplaced fractions in printed (9)–(10);
the executable algebra follows the preceding complete matrix and
the later Ghosh expressions.

**Expression.** If the two materials have Schelkunoff surface terms
``Z_{1i},Z_{1o},Z_{1m}`` and ``Z_{2i},Z_{2o},Z_{2m}``,

```math
\\begin{aligned}
Z_{out}&=Z_{2o}-\\frac{Z_{2m}^2}{Z_{1o}+Z_{2i}}, \\\\
Z_{in}&=Z_{1i}-\\frac{Z_{1m}^2}{Z_{1o}+Z_{2i}}, \\\\
Z_m&=\\frac{Z_{1m}Z_{2m}}{Z_{1o}+Z_{2i}}.
\\end{aligned}
```

**Numerical evaluation.** Semiconductor/metal resistivity contrast
can make direct subtraction lose the metal contribution. The same
Bessel and interface equations are evaluated with guard precision
in a task-local scope, then rounded to the requested precision.
No semiconductor-current or skin-depth approximation is introduced.

The source's physical penetration-depth and good-conductor restrictions
remain applicable. The N-screen cable construction uses the unexpanded
effective surfaces in Ghosh (10)–(11), then the existing loop-to-terminal
matrix transformation. Attached screens are not also counted as magnetic
insulation intervals; their separate radial shunt terms remain present.

**Reference.** [Ametani2004](@cite), (9)–(10);
[Ghosh2019](@cite), (11),(16)–(21); [Ghosh2022](@cite), (9)–(13).
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

The aggregate uses the uncombined material records. The ordinary five-argument
call is the one-material limit and does not reconstruct a missing profile.
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
    layers=[(r_in=a,r_ex=b,rho=rho_1,mu_r=mu_r1),
        (r_in=b_prime,r_ex=c,rho=rho_2,mu_r=mu_r2)]
    return formula(layers,jω)
end

function (formula::Formula{:Ametani2004})(layers::AbstractVector,jω::Complex)
    state=merge((;jω),_bonded_surfaces(layers,jω))
    return Functor{:Ametani2004, typeof(formula.routes), typeof(state)}(
        formula.routes,
        state
    )
end

function (formula::Formula{:Ametani2004})(r_in, r_ex, rho, mu_r, jω)
    return formula([(r_in=r_in,r_ex=r_ex,rho=rho,mu_r=mu_r)],jω)
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
