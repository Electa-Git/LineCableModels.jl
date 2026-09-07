function routes(identifier::Val{:Wedepohl1973})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Wedepohl1973})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Wedepohl1973}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Buried filamentary cable pair; self uses cable radius and mutual terms use conductor separation and burial depths. |
| Calculated quantities | Low-order closed self and mutual earth-return impedance approximations |
| Earth structure | Homogeneous earth below air. |
| Model and approximation | The Wedepohl–Wilcox reduction truncates the Pollaczek decomposition for small complex propagation-distance products. |
| Main source | L. M. Wedepohl and D. J. Wilcox (1973) |
| Citation key(s) | Primary: `:Wedepohl1973`; comparative equation witness: `:Guneri2018` |
| Evidence status | Original publication and comparative PDF equations checked |

**Expression.**

```math
Z_{e,ii}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[-\\ln
\\left(\\frac{e_c\\gamma_1r_i}{2}\\right)+\\frac12-
\\frac43\\gamma_1h_i\\right],
```

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[-\\ln
\\left(\\frac{e_c\\gamma_1d_{ij}}{2}\\right)+\\frac12-
\\frac23\\gamma_1H\\right],\\qquad e_c=e^{\\gamma_E}.
```

**Reference.** L. M. Wedepohl and D. J. Wilcox, “Transient Analysis of
Underground Power-Transmission Systems: System-Model and Wave-Propagation
Characteristics,” *Proceedings of the IEE*, 120, 253–260, 1973.
"""
function description(::Formula{:Wedepohl1973})
    "Wedepohl-Wilcox low-frequency underground approximation (1973)"
end

function propagation_constant(
        ::Val{:Wedepohl1973}, jω, permeability, permittivity
)
    (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Wedepohl1973})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Wedepohl1973), formula,
        rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Wedepohl-Wilcox low-frequency underground terms:

```math
Z_{e,ii}=\frac{j\omega\mu_0}{2\pi}
\left[-\ln\left(\frac{e_c\gamma_1r_i}{2}\right)+\frac12
-\frac43\gamma_1h_i\right],
```

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}
\left[-\ln\left(\frac{e_c\gamma_1d_{ij}}{2}\right)+\frac12
-\frac23\gamma_1(h_i+h_j)\right],\qquad e_c=e^{\gamma_E}.
```
"""
function earth_impedance(
        ::Val{:Wedepohl1973}, ::Val{:self}, functor, pair
)
    _require(pair, Val(:underground))
    geometry = _geometry(pair)
    return earth_impedance(Val(:Wedepohl1973),Val(:kernel),functor,
        geometry.y_ij,2geometry.h_i)
end

function earth_impedance(
        ::Val{:Wedepohl1973}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    geometry = _geometry(pair)
    return earth_impedance(Val(:Wedepohl1973),Val(:kernel),functor,
        geometry.d_ij,geometry.H)
end

function earth_impedance(::Val{:Wedepohl1973},::Val{:kernel},functor,d,H)
    state=functor.state
    e_c=exp(one(H)*Base.MathConstants.eulergamma)
    bracket=-log(e_c*state.gamma[2]*d/2)+one(H)/2-
        (2one(H)/3)*state.gamma[2]*H
    return _complex_result(state.jω,state.jω*state.mu[1]/(2*(one(H)*π))*bracket)
end

:Wedepohl1973
