function routes(identifier::Val{:Sunde1949})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        homogeneous = FormulaMethod(
            identifier, earth_impedance, Val(:homogeneous)
        ),
        two_layer = FormulaMethod(
            identifier, earth_impedance, Val(:two_layer)
        ),
        multilayer = FormulaMethod(
            identifier, earth_impedance, Val(:multilayer)
        ),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Sunde1949})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Sunde1949}) = Val(:zero)
media(::Formula{:Sunde1949}) = Val(:stratified)

function Formula(::Val{:Sunde1949}; approximation::Symbol=:integral,
        displacement_current::Bool=false, kwargs...)
    approximation in (:integral,:large_spacing) || throw(ArgumentError(
        ":Sunde1949 approximation must be :integral or :large_spacing"
    ))
    identifier=Val(:Sunde1949)
    defaults=routes(identifier)
    if approximation === :large_spacing
        defaults=merge(defaults,(
            mutual=FormulaMethod(identifier,earth_impedance,Val(:large_spacing)),
        ))
    end
    overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) || throw(ArgumentError(
        "unknown routes for earth-impedance formula :Sunde1949"
    ))
    values=merge(assumptions(identifier),(
        earth=displacement_current ? _full : _conductive,
    ))
    return Formula(identifier,merge(defaults,overrides),values)
end
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filamentary external wire positions; heights ``h_1,h_2`` and separation ``y``; finite conductor radius and insulation do not enter. |
| Calculated quantities | Mutual external inductance/impedance kernel for infinite overhead wires above two-layer earth; large-separation two-term expansion |
| Earth structure | Finite upper layer of depth ``d`` over a semi-infinite lower layer. |
| Model and approximation | Equation (4.55) already neglects the direct-distance logarithmic term for large separation and replaces horizontal by radial separation within the stated practical range. Equation (4.56) then expands ``F(u)`` as in the uniform-earth case and retains only its first two terms. |
| Main source | E. D. Sunde; formula appears in a 1968 corrected Dover republication of a work first published in 1949, but first-edition token identity was not inspected |
| Citation key(s) | `:Sunde1949` |
| Evidence status | Corrected-republication page image verified; original-1949 priority and token identity unresolved |

Overhead impedance for homogeneous or horizontally stratified earth;
the number of earth layers selects the homogeneous, two-layer, or recursive
kernel. The source's ground correction is combined with the ideal-ground
logarithm to supply the complete external matrix coefficient.

**Expression.** The two-layer kernel is

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+2\\int_0^\\infty F_{ij}^{S}(\\lambda)
\\cos(y_{ij}\\lambda)d\\lambda\\right],
```

```math
F_{ij}^{S}=\\frac{a_1+a_2+(a_1-a_2)e^{-2a_1d}}
{(a_1+a_2)(\\lambda+a_1)+(a_1-a_2)(\\lambda-a_1)e^{-2a_1d}}
e^{-\\lambda H}.
```

The same two-interface relation is applied recursively for more layers:

```math
b_N=a_N,\\qquad
b_m=a_m\\frac{b_{m+1}+a_m\\tanh(a_md_m)}
{a_m+b_{m+1}\\tanh(a_md_m)},\\qquad
F(\\lambda)=\\frac{1}{\\lambda+b_1}.
```

This is a boundary-matching extension of the two-layer coefficient.
The default uses the conductive source model. Select
`displacement_current=true` for the earlier engine's bulk full-admittivity
extension. Neither selection subtracts the air propagation constant.

Select `approximation=:large_spacing` for the two-term mutual correction
(4.56). The ideal-ground logarithm is still assembled separately; the
self route retains the integral, not the large-spacing expansion.

```math
Z_{g,ij}\\simeq\\frac{j\\omega\\mu_0}{\\pi}
\\left[F(0)\\frac{H}{H^2+x^2}
-F(0)^2\\frac{H^2-x^2}{(H^2+x^2)^2}\\right].
```

The decaying square-root branch and the SI coefficient are fixed by
the MKS witness below. Numerical quadrature evaluates the parent directly;
the source's graphical normalization and crossed operator are not required.

## Two-layer self witness

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Overhead line-current self geometry specified by height; conductor skin, insulation, and free-space geometric terms are excluded. |
| Calculated quantities | Two-layer self earth-return impedance, nonconducting-lower-region specialization, normalization, and graphical evaluator. |
| Earth structure | An upper layer of resistivity ``\\rho_1`` and thickness ``D`` above a lower half-space of resistivity ``\\rho_2``; ``\\rho_2=\\infty`` gives a nonconducting lower region. |
| Model and approximation | Conduction-only, fixed-permeability integral attributed to Sunde. The chart uses the approximate ``r_e\\sim0.1\\sqrt{f/R_e}``; it does not replace printed equation (21). |
| Main source | K. Iwamoto (August 1958), citing Sunde (1949). |
| Citation key(s) | Primary: `:Iwamoto1958b`; auxiliary evaluator: `:Iwamoto1958a`; attributed sources: `:Sunde1949`, `:Carson1926` |
| Evidence status | All 12 August pages checked. The root branch, unit scaling, and the January operator's step argument remain unresolved. |

**Reference.** E. D. Sunde, *Earth Conduction Effects in Transmission
Systems*, Dover, 1968.
"""
function description(::Formula{:Sunde1949})
    "Sunde conductive-earth overhead integral and large-spacing approximation (1949)"
end

function propagation_constant(::Val{:Sunde1949}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Sunde1949})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    length(rho) == 2 || throw(DimensionMismatch(
        ":Sunde1949 requires layer thicknesses for a stratified earth"
    ))
    return _homogeneous_functor(
        Val(:Sunde1949), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

function (formula::Formula{:Sunde1949})(
        rho, epsilon, mu, jω, Γ, segments, thickness
)
    length(rho) >= 2 || throw(DimensionMismatch(
        ":Sunde1949 requires air and at least one earth layer"
    ))
    return _stratified_functor(
        Val(:Sunde1949), formula,
        rho, epsilon, mu, jω, Γ, segments, thickness
    )
end

"""
$(TYPEDSIGNATURES)

Select Sunde's earth-return kernel from the number of physical earth layers.

One earth layer uses the homogeneous formula. Two and more layers apply
the same boundary relation. The two-layer and multilayer routes remain
individually replaceable.
"""
function earth_impedance(
        ::Val{:Sunde1949}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    count = length(functor.state.rho) - 1
    count == 1 && return functor.routes.homogeneous(functor, pair)
    count == 2 && return functor.routes.two_layer(functor, pair)
    return functor.routes.multilayer(functor, pair)
end

raw"""
Evaluate Sunde's homogeneous-earth overhead impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[\ln\frac{D_{ij}}{d_{ij}}+
2\int_0^\infty\frac{e^{-(h_i+h_j)\lambda}\cos(y_{ij}\lambda)}
{\lambda+\sqrt{\lambda^2+\gamma_1^2}}\,d\lambda\right],
```

with ``\gamma_1^2=j\omega\mu_0(\sigma_1+j\omega\varepsilon_1)``.
"""
function earth_impedance(
        ::Val{:Sunde1949}, ::Val{:homogeneous}, functor, pair
)
    return _homogeneous_overhead_coefficient(functor.state,pair)
end

raw"""
Evaluate Sunde's two-layer overhead earth-return impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[\ln\frac{D_{ij}}{d_{ij}}+
2\int_0^\infty F_{ij}^{S}(\lambda)\cos(y_{ij}\lambda)d\lambda\right],
```

```math
F_{ij}^{S}=\frac{a_1+a_2+(a_1-a_2)e^{-2a_1d}}
{(a_1+a_2)(\lambda+a_1)+(a_1-a_2)(\lambda-a_1)e^{-2a_1d}}
e^{-\lambda(h_i+h_j)},\qquad
a_m=\sqrt{\lambda^2+\gamma_m^2}.
```
"""
function earth_impedance(
        identifier::Val{:Sunde1949}, ::Val{:two_layer}, functor, pair
)
    return earth_impedance(identifier,Val(:multilayer),functor,pair)
end

"""
$(TYPEDSIGNATURES)

Apply the nonmagnetic two-layer boundary relation recursively to the
remaining soil stack. The kernel has length dimensions and reduces to
the homogeneous and two-layer coefficients without a separate unit factor.

The root in this Sunde selection is the bulk soil root, without the
air-reference subtraction used by the later Lee2014 and Xue2021 models.
"""
function spectral_kernel(::Val{:Sunde1949},lambda,state)
    return _layered_overhead_kernel(lambda,state)
end

function earth_impedance(::Val{:Sunde1949},::Val{:multilayer},functor,pair)
    return _layered_overhead_coefficient(functor.state,pair)
end

function earth_impedance(
        identifier::Val{:Sunde1949},::Val{:large_spacing},functor,pair
)
    _require(pair,Val(:overhead))
    pair.row != pair.column || throw(ArgumentError(
        "Sunde's large-spacing expansion is a mutual correction"
    ))
    length(functor.state.rho)==3 || throw(DimensionMismatch(
        "Sunde's surveyed large-spacing expansion requires two earth layers"
    ))
    state=functor.state; geometry=_geometry(pair)
    H=geometry.H; x=geometry.y_ij; distance2=H^2+x^2
    F0=spectral_kernel(identifier,zero(H),state)
    correction=F0*H/distance2-F0^2*(H^2-x^2)/distance2^2
    ideal=log(geometry.D_ij/geometry.d_ij)
    return state.jω*state.mu[1]/(2*(one(H)*π))*(ideal+2correction)
end

:Sunde1949
