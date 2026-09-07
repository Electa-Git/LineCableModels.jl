function routes(identifier::Val{:Pettersson1994})
    return (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Pettersson1994})
    (
        air = _full,
        earth = _full,
        permeability = vacuum_permeability,
        lossless_air = false
    )
end

propagation(::Val{:Pettersson1994}) = Val(:zero)

function Formula(::Val{:Pettersson1994};lossless_air::Bool=false,kwargs...)
    defaults=routes(Val(:Pettersson1994))
    overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) ||
        throw(ArgumentError("unknown routes for earth-impedance formula :Pettersson1994"))
    selected=merge(defaults,overrides)
    values=merge(assumptions(Val(:Pettersson1994)),
        (;air=lossless_air ? _lossless : _full,lossless_air))
    return Formula{:Pettersson1994,typeof(selected),typeof(values)}(selected,values)
end
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filamentary line source representing a sufficiently thin round wire of radius ``a``; no coating layer. |
| Calculated quantities | Generalized p.u.l. self/mutual series impedance of a thin wire above, at, or below a lossy planar interface |
| Earth structure | Two homogeneous half-spaces separated by a plane. |
| Model and approximation | The exact spectral integrand in (6) is replaced by the asymptotically matched expression (9), which is then integrated using identity (7). The self case ``x=0,y=h+a`` gives ``P\\simeq\\ln[1+1/(\\beta h)]`` in (11). Separate on-interface distances are printed in (14)–(15). |
| Main source | Pär Pettersson (1994 publication; 1993 conference manuscript) |
| Citation key(s) | `:Pettersson1994` |
| Evidence status | Original publication page images checked |

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+\\ln\\frac{
\\sqrt{(H+2/\\beta_\\gamma)^2+y_{ij}^2}}{D_{ij}}\\right],
\\qquad \\beta_\\gamma=\\sqrt{\\gamma_g^2-\\gamma_0^2}.
```

**Equivalent overhead images.** [Maaouni2001](@cite), equation (12), and
[Ametani2014](@cite), equations (17)–(18), give the same series image with
lossless air. The former `Ametani2014` selector now selects this formula
with `lossless_air=true`; no second numerical image evaluator is retained.
The image depth uses the air-referenced difference of bulk propagation
constants, not the conduction-only depth. Self correction uses zero
horizontal separation and the radius only in the direct-distance denominator.
Equations (10)–(11) also apply to buried wires after interchanging the
source and other half-spaces. On-interface wires use (14)–(15). These closed
images do not supply mixed or interface/off-interface mutual entries, which
require a common modal prescription. The coupled mode equation is not solved.

**Reference.** P. Pettersson, “Image Representation of Wave Propagation on
Wires Above, On and Under Ground,” *IEEE Transactions on Power Delivery*, 9,
1049–1055, 1994. DOI: 10.1109/61.296290.
"""
description(::Formula{:Pettersson1994}) =
    "Pettersson two-half-space and interface image approximation (1994)"

function propagation_constant(::Val{:Pettersson1994}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Pettersson1994})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    if formula.assumptions.lossless_air
        isinf(first(rho)) || throw(DomainError(first(rho), "this image selection requires lossless air"))
        isfinite(rho[2]) && rho[2]>zero(rho[2]) ||
            throw(DomainError(rho[2], "this image selection requires conducting earth"))
        iszero(jω) && throw(DomainError(jω,"this image selection requires nonzero frequency"))
    end
    return _homogeneous_functor(
        Val(:Pettersson1994), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate Pettersson's wideband image approximation for overhead earth-return
impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[
\ln\frac{D_{ij}}{d_{ij}}+M_{e,ij}\right],
\qquad
M_{e,ij}=\ln\frac{
\sqrt{(H+2/\beta_\gamma)^2+y_{ij}^2}}{D_{ij}},
```

```math
\beta_\gamma=\sqrt{\gamma_g^2-\gamma_0^2},\qquad
\gamma_m^2=j\omega\mu_m(\sigma_m+j\omega\varepsilon_m).
```

# Reference

P. Pettersson, "Image representation of wave propagation on wires above,
on and under ground," *IEEE Transactions on Power Delivery*, vol. 9,
pp. 1049-1055, 1994. DOI: 10.1109/61.296290.
"""
function earth_impedance(
        ::Val{:Pettersson1994}, ::Val{:mutual}, functor, pair
)
    state=functor.state
    images=pettersson_images(state,pair)
    return _complex_result(state.jω,state.jω*state.mu[1]/
        (2*(one(pair.separation)*π))*images.magnetic)
end

function earth_impedance(::Val{:Pettersson1994},::Val{:self},functor,pair)
    return earth_impedance(Val(:Pettersson1994),Val(:mutual),functor,pair)
end

:Pettersson1994
