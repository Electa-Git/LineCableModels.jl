function routes(identifier::Val{:Mariscotti2019})
    return (
        self=FormulaMethod(identifier,earth_potential_coefficient,Val(:self)),
        mutual=FormulaMethod(identifier,earth_potential_coefficient,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant))
end
assumptions(::Val{:Mariscotti2019})=(air=_full,earth=_full,permeability=_material)
propagation(::Val{:Mariscotti2019})=Val(:explicit)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Scalar two-conductor setup; self/mutual positions. Filamentary external field; no coating. |
| Calculated quantities | Self/mutual generalized-potential coefficient and source shunt-admittance extraction |
| Earth structure | Air over homogeneous lossy earth. |
| Model and approximation | Generalized potential within the source's thin-wire/modal model; no extra closed-form approximation is applied to the displayed integral. |
| Main source | Andrea Mariscotti (2019), building on D’Amore–Sarto |
| Citation key(s) | `:Mariscotti2019` |
| Evidence status | Original publication page images checked |

**Expression.** Charge continuity gives
``P=V/Q=j\\omega V/(\\gamma_{source}I)``, not the reciprocal of
the source's ``V/I_t`` entry. For the medium ``m`` containing both wires,

```math
P_{ij}=\\frac{j\\omega}{2\\pi\\kappa_m}
\\left[K_0(\\chi_m d_{ij})-K_0(\\chi_m D_{ij})
+2\\int_0^\\infty
\\frac{\\gamma_m^2e^{-Hq_m}\\cos(x_{ij}u)}
{\\gamma_1^2q_2+\\gamma_2^2q_1}\\,du\\right],
```

where ``\\kappa_m=\\sigma_m+j\\omega\\varepsilon_m``,
``\\gamma_m^2=j\\omega\\mu_m\\kappa_m``,
``\\chi_m^2=\\gamma_m^2+k_x^2``, and
``q_m=\\sqrt{u^2+\\chi_m^2}``.
The engine wavenumber is ``k_x=j\\gamma_{source}``.

**Numerical interpretation.** The image uses ``K_0``, consistent
with the two-dimensional potential Green function; the printed
``K_1`` remains in the survey transcription. The default uses the
lossless-air reference from ``\\gamma_{source}\\simeq jk_1`` and
the source definition ``k_1^2=\\omega^2\\mu_1\\varepsilon_1``.
The extra imaginary factor in the radical printed in (17) is not
used to change that definition. An explicit engine wavenumber
retains the prescribed-propagation potential. No modal root is solved.

Both same-medium placements are available. Equations (14)–(15)
do not provide a cross-interface source/receiver coefficient, so
mixed pairs are rejected. The complete potential matrix, including
separate radial insulation contributions, is inverted only after
assembly. This source normalization is distinct from Xue's buried
potential kernel and must not be substituted for it silently.

**Reference.** [Mariscotti2019](@cite), (8),(12)–(21).
"""
description(::Formula{:Mariscotti2019})="Mariscotti same-medium generalized scalar potential (2019)"
propagation_constant(::Val{:Mariscotti2019},s,mu,epsilon)=
    (Γ=sqrt(-s*s*mu*epsilon),squared=-s*s*mu*epsilon)
function (formula::Formula{:Mariscotti2019})(rho,epsilon,mu,s,Γ,segments=nothing)
    return _scalar_potential_functor(Val(:Mariscotti2019),formula,rho,epsilon,mu,s,Γ,segments)
end
function earth_potential_coefficient(::Val{:Mariscotti2019},::Val{:self},functor,pair)
    pair.row==pair.column && pair.heights[1]==pair.heights[2] ||
        throw(ArgumentError("Mariscotti self coefficient requires repeated wire geometry"))
    0<pair.separation<abs(pair.heights[1]) ||
        throw(DomainError(pair,"the complete wire section must lie in one half-space"))
    return earth_potential_coefficient(Val(:Mariscotti2019),Val(:coefficient),functor,pair)
end
function earth_potential_coefficient(::Val{:Mariscotti2019},::Val{:mutual},functor,pair)
    pair.row!=pair.column || throw(ArgumentError("mutual coefficient requires distinct wires"))
    return earth_potential_coefficient(Val(:Mariscotti2019),Val(:coefficient),functor,pair)
end
function earth_potential_coefficient(::Val{:Mariscotti2019},::Val{:coefficient},functor,pair)
    _placement(pair)===Val(:mixed) &&
        throw(ArgumentError("Mariscotti (14)–(15) do not supply a mixed source/receiver kernel"))
    all(h->!iszero(h),pair.heights) || throw(DomainError(pair,"wire centres must lie strictly in one half-space"))
    geometry=_geometry(pair)
    return _scalar_potential_coefficient(functor.state,pair.layers[1]==1 ? 1 : 2,
        geometry.d_ij,geometry.D_ij,geometry.H,geometry.y_ij)
end

:Mariscotti2019
