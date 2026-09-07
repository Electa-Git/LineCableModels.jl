function routes(identifier::Val{:DiLorenzo2023})
    return (
        self=FormulaMethod(identifier,earth_potential_coefficient,Val(:self)),
        mutual=FormulaMethod(identifier,earth_potential_coefficient,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:DiLorenzo2023})=(
    air=_full,earth=_full,permeability=_material
)
propagation(::Val{:DiLorenzo2023})=Val(:zero)
media(::Formula{:DiLorenzo2023})=Val(:stratified)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Each cable has core radius ``r_1``, sheath radii ``r_2,r_3``, and outer insulation radius ``r_4``; inner and outer insulation coefficients are included in ``P'_L``. |
| Calculated quantities | Seabed-return self/mutual potential coefficients and assembled phase-domain per-unit-length admittance for submarine cables |
| Earth structure | Air over finite-depth seawater over semi-infinite seabed. |
| Model and approximation | The parent electromagnetic problem is reduced by the quasi-TEM assumption. Equations (33)–(35) remain infinite spectral integrals and require numerical quadrature; no fit or truncation order is specified. Equation (39) is the ideal coaxial lossless-insulation coefficient used by the source. |
| Main source | G. Di Lorenzo, E. Stracqualursi, M. Marzinotto, J. Brandao Faria, and R. Araneo (2023) |
| Citation key(s) | `:DiLorenzo2023` |
| Evidence status | Original publication page images checked; DOI retained |

**Expression.** The source's bottom-layer coefficient is evaluated as

```math
P_{ij}=\\frac{j\\omega}{2\\pi\\kappa_2}
\\left[
K_0(\\gamma_2d_{ij})+
\\int_0^\\infty\\left(\\frac{R_m}{\\alpha_2}+G_b^0\\right)
e^{-\\alpha_2(h_i+h_j-2h_s)}\\cos(\\lambda q_{ij})\\,d\\lambda
\\right].
```

Here ``R_m=-(s_{10}d_{21}-d_{10}s_{21}E)/D_m``,
``E=e^{-2\\alpha_1h_s}``, and
``D_m=s_{10}s_{21}-d_{10}d_{21}E``.
The electric correction is the numerator printed in (34), without its
burial exponential, divided by ``D_mD_e``, where

```math
\\begin{aligned}
D_e&=A_{10}A_{21}+\\Delta_{10}\\Delta_{21}E,\\
A_{21}&=\\alpha_1\\gamma_2^2\\mu_1+
         \\alpha_2\\gamma_1^2\\mu_2,\\
\\Delta_{21}&=\\alpha_1\\gamma_2^2\\mu_1-
         \\alpha_2\\gamma_1^2\\mu_2.
\\end{aligned}
```

This denominator follows from the eight Hertz-potential boundary
conditions (27); the printed numerator alone does not have the spectral
units of ``F_3``. Independent boundary-system solutions, zero-water-depth
and homogeneous-layer limits test the normalization. Permeability and
bulk-wave factors are scaled in the numerical quotients to avoid
underflow without altering the kernels.

**Assembly and limits.** Both source and observation must lie in the
bottom half-space, below a finite water layer. Independent positive
layer permeabilities and displacement current are retained; the common
axial input is zero as in the source qTEM reduction. Self terms sample
at the outer cable radius; mutual terms use the actual transverse
distance. Coaxial insulation remains in the existing local potential
assembly. The full phase-domain matrix is inverted once, reproducing
(36)–(40), not elementwise reciprocals of the external coefficients.

**Reference.** [DiLorenzo2023](@cite), equations (20), (25)–(40).
The series kernel is already covered by Tsiamitros2008 under the
restrictions recorded in the survey; the electric correction is not
inferred from the series impedance alone.
"""
description(::Formula{:DiLorenzo2023}) =
    "Di Lorenzo et al. three-medium seabed potential coefficient (2023)"

function propagation_constant(::Val{:DiLorenzo2023},s,mu,epsilon)
    return (Γ=zero(s),squared=zero(s))
end

function (formula::Formula{:DiLorenzo2023})(rho,epsilon,mu,s,Γ,segments,thickness)
    _check(rho,epsilon,mu)
    length(rho)==length(thickness)==3 ||
        throw(DimensionMismatch(":DiLorenzo2023 requires air, one finite layer, and a bottom half-space"))
    isinf(rho[1]) && all(x->isfinite(x)&&x>0,rho[2:3]) ||
        throw(DomainError(rho,":DiLorenzo2023 requires lossless air and two conducting earth media"))
    all(x->isfinite(x)&&x>0,epsilon) && all(x->isfinite(x)&&x>0,mu) ||
        throw(DomainError((epsilon,mu),":DiLorenzo2023 requires positive finite material constants"))
    isinf(thickness[1]) && isinf(thickness[3]) &&
        isfinite(thickness[2]) && thickness[2]>=0 ||
        throw(DomainError(thickness,":DiLorenzo2023 requires one finite nonnegative water depth"))
    isfinite(s) && !iszero(s) && iszero(real(s)) ||
        throw(DomainError(s,":DiLorenzo2023 requires nonzero real frequency"))
    return _stratified_functor(
        Val(:DiLorenzo2023),formula,rho,epsilon,mu,s,Γ,segments,thickness)
end

function earth_potential_coefficient(::Val{:DiLorenzo2023},::Val{:kernel},lambda,state)
    roots=map(g->spectral_root(lambda^2+g,state.jω),state.gamma_medium_squared)
    scale_a=maximum(abs,roots)
    a0,a1,a2=roots./scale_a
    mu0,mu1,mu2=state.mu./maximum(state.mu)
    g0,g1,g2=state.gamma_medium_squared./maximum(abs,state.gamma_medium_squared)
    decay=exp(-2roots[2]*state.thickness[2])
    s10=mu0*a1+mu1*a0; d10=mu0*a1-mu1*a0
    s21=mu2*a1+mu1*a2; d21=mu2*a1-mu1*a2
    A10=a0*g1*mu0+a1*g0*mu1; D10=a0*g1*mu0-a1*g0*mu1
    A21=a1*g2*mu1+a2*g1*mu2; D21=a1*g2*mu1-a2*g1*mu2
    magnetic=s10*s21-d10*d21*decay
    electric=A10*A21+D10*D21*decay
    numerator=(g2-g1)*(s10*A10-d10*D10*decay^2)-
        2mu0*mu1*decay*(a0^2*g1*(g2-g1)+a1^2*g0*(g2+g1)-2a1^2*g1*g2)
    G=2mu1*mu2*a2*numerator/(magnetic*electric)/scale_a
    R=-(s10*d21-d10*s21*decay)/magnetic
    return (;R,G,alpha=roots[3])
end

function earth_potential_coefficient(::Val{:DiLorenzo2023},::Val{:mutual},functor,pair)
    all(==(3),pair.layers) ||
        throw(ArgumentError(":DiLorenzo2023 requires both conductors in the bottom half-space"))
    state=functor.state
    di=-pair.heights[1]-state.thickness[2]
    dj=-pair.heights[2]-state.thickness[2]
    x=abs(pair.separation)
    d=hypot(x,di-dj); H=di+dj
    all(isfinite,(di,dj,x,d)) && di>0 && dj>0 && d>0 ||
        throw(DomainError((pair.heights,pair.separation),
            ":DiLorenzo2023 requires positive burial below the water layer and nonzero distance"))
    correction=_quadrature(state) do t
        iszero(exp(-t)) && return zero(state.jω)
        lambda=t/H
        kernel=earth_potential_coefficient(Val(:DiLorenzo2023),Val(:kernel),lambda,state)
        (kernel.R/kernel.alpha+kernel.G)*exp(-kernel.alpha*H)*cos(lambda*x)/H
    end
    direct=special_besselk(0,state.gamma[3]*d)
    kappa=state.sigma[3]+state.jω*state.epsilon[3]
    return _complex_result(state.jω,state.jω*(direct+correction)/
        (2*(one(H)*π)*kappa))
end

:DiLorenzo2023
