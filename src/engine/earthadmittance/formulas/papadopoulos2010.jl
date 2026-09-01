function routes(identifier::Val{:Papadopoulos2010})
    (
        self = formula_method(identifier, earth_potential_coefficient, Val(:self)),
        mutual = formula_method(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = formula_method(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Papadopoulos2010})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Papadopoulos2010}) = Val(:explicit)
"""
$(TYPEDSIGNATURES)

**Identification.** Homogeneous-earth underground potential coefficient with
explicit longitudinal propagation.

**Expression.**

```math
P_{e,ij}=\\frac{j\\omega}{2\\pi(\\sigma_1+j\\omega\\varepsilon_1)}
(\\Delta_1+2S_P),
```

```math
\\Delta_1=\\int_0^\\infty
\\frac{e^{-|h_i-h_j|\\alpha_1}-e^{-H\\alpha_1}}{\\alpha_1}
\\cos(y\\lambda)d\\lambda,\\quad
S_P=\\int_0^\\infty\\frac{e^{-H\\alpha_1}}
{\\alpha_0+(\\gamma_0^2/\\gamma_1^2)\\alpha_1}
\\cos(y\\lambda)d\\lambda.
```

**Reference.** T. A. Papadopoulos, D. A. Tsiamitros, and G. K. Papagiannis,
“Impedances and Admittances of Underground Cables for the Homogeneous Earth
Case,” *IEEE Transactions on Power Delivery*, 25(2), 961–969, 2010.
DOI: 10.1109/TPWRD.2009.2034797.
"""
function description(::Formula{:Papadopoulos2010})
    "Papadopoulos et al. homogeneous-earth underground potential coefficient (2010)"
end

function propagation_constant(
        ::Val{:Papadopoulos2010}, jω, permeability, permittivity
)
    squared = oftype(jω, (-jω^2) * permeability * permittivity)
    return (Γ = sqrt(squared), squared)
end

function (formula::Formula{:Papadopoulos2010})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Papadopoulos2010), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Papadopoulos et al. underground potential coefficient:

```math
P_{e,ij}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
(\Delta_1+2S_P),
```

```math
\Delta_1=\int_0^\infty\frac{e^{-|h_i-h_j|\alpha_1}-e^{-H\alpha_1}}
{\alpha_1}\cos(y\lambda)d\lambda,\qquad
S_P=\int_0^\infty\frac{e^{-H\alpha_1}}
{\alpha_0+(\gamma_0^2/\gamma_1^2)\alpha_1}
\cos(y\lambda)d\lambda,
```

where ``\alpha_m=\sqrt{\lambda^2+\gamma_m^2+k_x^2}``.

The JSON corpus rendered the denominator of ``S_P`` as the product
``\alpha_0\alpha_1``. Ametani et al. (2021), Eq. 2.69, has the sum shown
above; the product is logarithmically singular at ``\lambda=0`` under the
stated ``k_x`` and therefore cannot represent the cited author formula.
"""
function earth_potential_coefficient(
        ::Val{:Papadopoulos2010}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    radial_squared = (
        state.gamma_medium_squared[1] + state.gamma_squared,
        state.gamma_medium_squared[2] + state.gamma_squared
    )
    radial = sqrt(radial_squared[2])
    delta_1 = special_besselk(0, radial * geometry.d_ij) -
              special_besselk(0, radial * geometry.D_ij)
    S_P = _quadrature(state) do lambda
        alpha_0 = sqrt(lambda^2 + radial_squared[1])
        alpha_1 = sqrt(lambda^2 + radial_squared[2])
        denominator = alpha_0 +
                      state.gamma_medium_squared[1] /
                      state.gamma_medium_squared[2] * alpha_1
        exp(-geometry.H * alpha_1) / denominator *
        cos(geometry.y_ij * lambda)
    end
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    return state.jω / (2π * kappa) * (delta_1 + 2 * S_P)
end

:Papadopoulos2010
