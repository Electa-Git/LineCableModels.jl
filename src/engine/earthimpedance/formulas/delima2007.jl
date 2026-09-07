function routes(identifier::Val{:DeLima2007})
    return (
        self=FormulaMethod(identifier,earth_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,earth_impedance,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:DeLima2007}) = (
    air=_lossless,earth=_full,permeability=vacuum_permeability
)
propagation(::Val{:DeLima2007}) = Val(:zero)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel overhead conductors represented by heights and separation. No conductor radius or insulation parameter enters this ground correction. |
| Calculated quantities | Self and mutual overhead ground-return corrections per unit length |
| Earth structure | Linear homogeneous isotropic half-space with a planar air interface. |
| Model and approximation | Integral representation within the declared TEM/quasi-TEM model, with a separate empirical constitutive approximation (7)–(8). The authors evaluate the integrals by Gauss–Kronrod quadrature; the numerical evaluation is not a new physical kernel or an analytical asymptote. In the conduction-only limit the source sets ``\\varepsilon_s=0`` and ``\\sigma_s=\\sigma_0``, hence ``\\eta_s=\\sqrt{j\\omega\\mu_0\\sigma_0}``. |
| Main source | Antonio Carlos Siqueira de Lima and Carlos Portela (2007), extending Carson-type expressions to complex frequency-dependent soil parameters; the constitutive soil model is attributed to earlier Portela work |
| Citation key(s) | `:DeLima2007` |
| Evidence status | Original PDF equations checked against page images; longitudinal-reduction prescription unresolved |

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Cable external radius ``r`` for self evaluation; burial depths ``h_i,h_j``; horizontal separation ``d_{ij}``; direct and reflected distances ``d,D,D_c``. |
| Calculated quantities | Buried-cable self and mutual ground-return impedance per length |
| Earth structure | Homogeneous isotropic linear half-space under air. |
| Model and approximation | Integral representation within the source's reduced field model, combined with an empirical joint soil law. No analytical approximation of these integrals is introduced; Gauss–Kronrod is their numerical evaluation method. The stated conduction-only limit sets ``\\varepsilon_s=0`` and ``\\sigma_s=\\sigma_0``, so ``\\eta_s=\\sqrt{j\\omega\\mu_0\\sigma_0}``. |
| Main source | Antonio Carlos Siqueira de Lima and Carlos Portela (2007), extending the Pollaczek-type earth representation to complex frequency-dependent soil |
| Citation key(s) | `:DeLima2007` |
| Evidence status | Original PDF equations checked against page images; main-text/appendix Bessel-order disagreement and appendix dependencies unresolved |

**Expression.** The final kernels use
`η² = jωμ₀(σ(f) + jωε(f))` with no independent air wave term.
For overhead pairs, the finite-earth correction is

```math
z_g=\\frac{j\\omega\\mu_0}{\\pi}
\\int_0^\\infty
\\frac{e^{-H\\lambda}\\cos(x\\lambda)}
{\\lambda+\\sqrt{\\lambda^2+\\eta^2}}\\,d\\lambda.
```

The engine adds the separate ideal-ground logarithm from source equation
(11), since this family returns the complete external coefficient.
The overhead self correction uses `x=0`, with the radius only in the
ideal-ground self logarithm.

For buried pairs, the Bessel/image integral uses the source appendix's
`K₀` image term, equation (30). The main-text `K₁` image would not
recover the stated conductive Pollaczek limit. This numerical selection is
explicit; the source Markdown retains both printed witnesses.

**Soil relation.** Material values are supplied by the selected soil law
before this evaluator runs; they are not fitted a second time here.
The source joint power law is already represented by `Portela1999` with

```math
\\alpha=\\text{exponent},\\qquad
\\beta=\\frac{\\Delta_i\\cot(\\pi\\alpha/2)}
{10^{-6}(2\\pi\\,10^6)^\\alpha}.
```

A mixed overhead/buried route is not inferred from these two source records.
The zero longitudinal value describes their final transverse kernels, not
a claimed prescription for the appendix's unknown modal constant.

**Reference.** [DeLima2007](@cite), equations (1)–(4), (7)–(8), (11),
and appendix equation (30).
"""
description(::Formula{:DeLima2007}) =
    "De Lima–Portela impedance with complex frequency-dependent soil (2007)"

function propagation_constant(::Val{:DeLima2007},jω,permeability,permittivity)
    return (Γ=zero(jω),squared=zero(jω))
end

function (formula::Formula{:DeLima2007})(rho,epsilon,mu,jω,Γ,segments=nothing)
    return _homogeneous_functor(
        Val(:DeLima2007),formula,rho,epsilon,mu,jω,Γ,segments
    )
end

function earth_impedance(::Val{:DeLima2007},::Val{:self},functor,pair)
    placement=_placement(pair)
    placement===Val(:underground) && return earth_impedance(
        Val(:Pollaczek1926),Val(:underground),functor,pair
    )
    _require(pair,Val(:overhead))
    return _homogeneous_overhead_coefficient(functor.state,pair)
end

function earth_impedance(::Val{:DeLima2007},::Val{:mutual},functor,pair)
    placement=_placement(pair)
    placement===Val(:overhead) && return earth_impedance(
        Val(:Carson1926),Val(:mutual),functor,pair
    )
    placement===Val(:underground) && return earth_impedance(
        Val(:Pollaczek1926),Val(:underground),functor,pair
    )
    throw(ArgumentError(":DeLima2007 does not supply a mixed overhead/buried kernel"))
end

:DeLima2007
