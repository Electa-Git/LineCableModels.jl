function routes(identifier::Val{:DeLima2018})
    return (
        self=FormulaMethod(identifier,earth_potential_coefficient,Val(:self)),
        mutual=FormulaMethod(identifier,earth_potential_coefficient,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:DeLima2018})=(
    air=_full,earth=_full,permeability=vacuum_permeability,approximation=:image
)
propagation(::Val{:DeLima2018})=Val(:explicit)

function Formula(::Val{:DeLima2018};approximation::Symbol=:image,kwargs...)
    approximation in (:image,:quasi_full_wave) ||
        throw(ArgumentError(":DeLima2018 approximation must be :image or :quasi_full_wave"))
    identifier=Val(:DeLima2018); defaults=routes(identifier); overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) ||
        throw(ArgumentError("unknown routes for earth-admittance formula :DeLima2018"))
    selected=merge(defaults,overrides)
    values=merge(assumptions(identifier),(;approximation))
    return Formula{:DeLima2018,typeof(selected),typeof(values)}(selected,values)
end


"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite circular conductor of radius ``r`` at height/depth ``h``; no finite insulation layer in the displayed formula. |
| Calculated quantities | Full-wave parent and quasi-full-wave per-unit-length admittance of one overhead or bare buried conductor, normalized as the scalar line admittance ``Y=\\gamma/Z_c`` |
| Earth structure | Two homogeneous half-spaces. |
| Model and approximation | qFW replaces the unknown longitudinal root by a predefined image-derived value inside the full-wave spectral quantities. It does not apply scalar reciprocals to potential coefficients, and it does not replace the remaining integral by a closed form. |
| Main source | A. C. S. de Lima, A. P. C. Magalhães, P. E. D. Rocha, R. A. Meyberg, and M. T. C. de Barros (2018) |
| Citation key(s) | `:DeLima2018` |
| Evidence status | Original publication page images and Appendix D image expressions checked |

**Expression.** Voltage definition (10)–(11) gives
``Z_c=\\gamma_p(\\Lambda_1-S_4)/(2\\pi\\kappa_1)``, hence
``Y=2\\pi\\kappa_1/(\\Lambda_1-S_4)`` and the returned scalar
coefficient is ``P=j\\omega(\\Lambda_1-S_4)/(2\\pi\\kappa_1)``.
The inverse follows independently by integrating the scalar and vector
potentials in (3). It is also explicit in Appendix (23). The physical
material factor is ``\\kappa_1=\\sigma_1+j\\omega\\varepsilon_1``.

For the image selection, the same Appendix gives
``P_{image}=j\\omega[\\ln(2h/r)-\\bar S_4]/(2\\pi\\kappa_1)``.
The shared coefficient evaluator belongs to EarthImpedance.DeLima2018;
there is no second copy of its image or spectral kernels here.

**Image expression.** The default `approximation=:image` evaluates
Appendix D, equations (23)–(24), which are printed in the same paper.
The separate conductor impedance is supplied by the internal family.
For the medium containing the wire, ``\\kappa_1=\\sigma_1+j\\omega\\varepsilon_1``,
``N=\\gamma_2^2/\\gamma_1^2``, ``\\beta=\\sqrt{\\gamma_2^2-\\gamma_1^2}``,
``D=\\sqrt{4h^2+r^2}``, and ``A=(N+1)/(\\beta D)``:

```math
\\begin{aligned}
\\bar S_1&=\\ln\\left(1+\\frac{2}{\\beta D}\\right),\\\\
\\bar S_2&=\\frac{2}{N+1}\\ln(1+A),\\\\
\\bar S_4&=\\frac{2\\ln2}{N+1}
+\\frac{2N}{N+1}\\ln\\left(1+\\frac{1}{1+2A}\\right).
\\end{aligned}
```

The last expression is algebraically (24), arranged to retain the
cancellation at large earth-to-air material ratios.

**Propagation input.** Select `approximation=:quasi_full_wave`
with an explicit backend input `k_x=j*gamma_image`. The image estimate is
``\\gamma_{image}=\\sqrt{(z_i+Z_{image})j\\omega/P_{image}}``.
The source's internal impedance ``z_i`` must be supplied explicitly
when forming that estimate; it is not silently set to zero. The helper
`EarthImpedance.propagation_constant(Val(:DeLima2018), Val(:image),
parameters, z_i, jomega)` returns the backend wavenumber and its square.
Here `parameters` contains the image `impedance` and `potential`
coefficients returned by the shared `Val(:parameters)` route.

The qFW selection evaluates the integral coefficients at that prescribed
wavenumber. It does not solve the full-wave modal equation, select between
fast and TL modes, or claim that the image estimate satisfies the exact
modal equation. A missing qFW propagation input is rejected.


**Scope.** This is a one-exterior-unit scalar line coefficient with the
source's interface-referenced voltage, not a mutual Maxwell-potential
formula. Normal assembly still uses ``Y=j\\omega P^{-1}``.
A second external wire is rejected. The same bare-wire, nonmagnetic
two-half-space and conductor-loss assumptions as the impedance source
apply. No automatic full-wave modal root selection is introduced.

**Reference.** [DeLima2018](@cite), equations (3), (10)–(16), and
Appendix D, equations (23)–(24).
"""
description(::Formula{:DeLima2018}) =
    "De Lima et al. single-wire image and prescribed qFW potential (2018)"

function propagation_constant(::Val{:DeLima2018},s,mu,epsilon)
    return (Γ=zero(s),squared=zero(s))
end

function (formula::Formula{:DeLima2018})(rho,epsilon,mu,s,Γ,segments=nothing)
    EarthImpedance.earth_impedance(
        Val(:DeLima2018),Val(:validate),formula,rho,epsilon,mu,s,Γ)
    return _homogeneous_functor(Val(:DeLima2018),formula,rho,epsilon,mu,s,Γ,segments)
end

function earth_potential_coefficient(::Val{:DeLima2018},::Val{:self},functor,pair)
    return EarthImpedance.earth_impedance(
        Val(:DeLima2018),Val(:parameters),functor,pair).potential
end

function earth_potential_coefficient(::Val{:DeLima2018},::Val{:mutual},functor,pair)
    throw(ArgumentError(":DeLima2018 supplies no mutual potential coefficient"))
end

:DeLima2018
