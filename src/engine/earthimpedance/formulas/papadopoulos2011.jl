function routes(identifier::Val{:Papadopoulos2011})
    (
        self = formula_method(identifier, earth_impedance, Val(:self)),
        mutual = formula_method(identifier, earth_impedance, Val(:mutual)),
        Γ = formula_method(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Papadopoulos2011})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Papadopoulos2011}) = Val(:explicit)
media(::Formula{:Papadopoulos2011}) = Val(:stratified)
"""
$(TYPEDSIGNATURES)

**Identification.** Two-layer underground mutual impedance for conductors in
the finite upper soil layer.

**Expression.**

```math
Z'_{e,ij}=\\frac{j\\omega\\mu_1}{2\\pi}\\int_0^\\infty
F_{ij}^{strat}(\\lambda)\\cos(y_{ij}\\lambda)d\\lambda,
```

```math
F_{ij}^{strat}=\\frac{
s_{10}s_{21}e^{-\\alpha_1|h_i-h_j|}+
s_{10}d_{21}e^{-\\alpha_1(2d-h_i-h_j)}-
d_{10}s_{21}e^{-\\alpha_1H}-
d_{10}d_{21}e^{-\\alpha_1(2d-|h_i-h_j|)}}
{\\alpha_1(s_{10}s_{21}+d_{10}d_{21}e^{-2\\alpha_1d})}.
```

**Reference.** T. A. Papadopoulos, D. A. Tsiamitros, and G. K. Papagiannis,
“Earth Return Admittances and Impedances of Underground Cables in
Non-Homogeneous Earth,” *IET Generation, Transmission & Distribution*, 5(2),
161–171, 2011.
"""
function description(::Formula{:Papadopoulos2011})
    "Papadopoulos et al. two-layer underground mutual impedance (2011)"
end

function propagation_constant(
        ::Val{:Papadopoulos2011}, jω, permeability, permittivity
)
    squared = oftype(jω, (-jω^2) * permeability * permittivity)
    return (Γ = sqrt(squared), squared)
end

function (formula::Formula{:Papadopoulos2011})(
        rho, epsilon, mu, jω, Γ, segments, thickness
)
    length(rho) == 3 || throw(DimensionMismatch(
        ":Papadopoulos2011 requires air and exactly two earth layers"
    ))
    return _stratified_functor(
        Val(:Papadopoulos2011), formula,
        rho, epsilon, mu, jω, Γ, segments, thickness
    )
end

raw"""
Evaluate the Papadopoulos et al. two-layer underground impedance:

```math
Z'_{e,ij}=\frac{j\omega\mu_1}{2\pi}\int_0^\infty
F_{ij}^{strat}(\lambda)\cos(y_{ij}\lambda)d\lambda,
```

```math
F_{ij}^{strat}=\frac{
s_{10}s_{21}e^{-\alpha_1|h_i-h_j|}+
s_{10}d_{21}e^{-\alpha_1(2d-h_i-h_j)}-
d_{10}s_{21}e^{-\alpha_1(h_i+h_j)}-
d_{10}d_{21}e^{-\alpha_1(2d-|h_i-h_j|)}}
{\alpha_1(s_{10}s_{21}+d_{10}d_{21}e^{-2\alpha_1d})}.
```

Here ``s_{mn}=\mu_n\alpha_m+\mu_m\alpha_n``, 
``d_{mn}=\mu_m\alpha_n-\mu_n\alpha_m``, and
``\alpha_m=\sqrt{\lambda^2+\gamma_m^2+k_x^2}``. Conductors must lie in the
finite top earth layer.
"""
function earth_impedance(
        ::Val{:Papadopoulos2011}, ::Val{:mutual}, functor, pair
)
    pair.layers == (2, 2) || throw(ArgumentError(
        ":Papadopoulos2011 requires both conductors in the top earth layer"
    ))
    state = functor.state
    geometry = _geometry(pair)
    d = state.thickness[2]
    difference = abs(geometry.h_i - geometry.h_j)
    radial = sqrt(
        state.gamma_medium_squared[2] + state.gamma_squared
    )
    direct = special_besselk(0, radial * geometry.d_ij)
    integral = _quadrature(state) do lambda
        alpha0 = sqrt(lambda^2 + state.gamma_medium_squared[1] + state.gamma_squared)
        alpha1 = sqrt(lambda^2 + state.gamma_medium_squared[2] + state.gamma_squared)
        alpha2 = sqrt(lambda^2 + state.gamma_medium_squared[3] + state.gamma_squared)
        s10 = state.mu[1] * alpha1 + state.mu[2] * alpha0
        d10 = state.mu[2] * alpha0 - state.mu[1] * alpha1
        s21 = state.mu[2] * alpha2 + state.mu[3] * alpha1
        d21 = state.mu[3] * alpha1 - state.mu[2] * alpha2
        numerator = s10 * s21 * exp(-alpha1 * difference) +
                    s10 * d21 * exp(-alpha1 * (2d - geometry.H)) -
                    d10 * s21 * exp(-alpha1 * geometry.H) -
                    d10 * d21 * exp(-alpha1 * (2d - difference))
        denominator = alpha1 * (
            s10 * s21 + d10 * d21 * exp(-2alpha1 * d)
        )
        F = numerator / denominator
        direct_spectrum = exp(-alpha1 * difference) / alpha1
        (F - direct_spectrum) * cos(geometry.y_ij * lambda)
    end
    return state.jω * state.mu[2] / (2π) * (direct + integral)
end

:Papadopoulos2011
