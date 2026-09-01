function routes(identifier::Val{:Theethayi2007})
    (
        self = formula_method(identifier, earth_impedance, Val(:self)),
        mutual = formula_method(identifier, earth_impedance, Val(:mutual)),
        Γ = formula_method(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Theethayi2007})
    (
        air = _lossless,
        earth = _full,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Theethayi2007}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Logarithmic-exponential underground approximation with
earth displacement current retained.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\left(\\frac{1+\\gamma_1R_{ab}}{\\gamma_1R_{ab}}\\right)+
\\frac{2e^{-H|\\gamma_1|}}{4+\\gamma_1^2R_{ab}^2}\\right].
```

**Reference.** N. Theethayi, R. Thottappillil, M. Paolone, C. A. Nucci, and
F. Rachidi, “External Impedance and Admittance of Buried Horizontal Wires for
Transient Studies Using Transmission Line Analysis,” *IEEE Transactions on
Dielectrics and Electrical Insulation*, 14(3), 751–761, 2007.
DOI: 10.1109/TDEI.2007.369540.
"""
function description(::Formula{:Theethayi2007})
    "Theethayi logarithmic-exponential underground approximation (2007)"
end

function propagation_constant(::Val{:Theethayi2007}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Theethayi2007})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Theethayi2007), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Theethayi et al. underground approximation:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[
\ln\left(\frac{1+\gamma_1R_{ab}}{\gamma_1R_{ab}}\right)+
\frac{2e^{-(h_i+h_j)|\gamma_1|}}{4+\gamma_1^2R_{ab}^2}\right],
```

where ``\gamma_1^2=j\omega\mu_0(\sigma_1+j\omega\varepsilon_1)``.
"""
function earth_impedance(
        ::Val{:Theethayi2007}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    pair.row == pair.column || _require_horizontal_separation(pair)
    state = functor.state
    geometry = _geometry(pair)
    gamma = state.gamma[2]
    argument = gamma * geometry.y_ij
    correction = 2exp(-geometry.H * abs(gamma)) / (4 + argument^2)
    return state.jω * state.mu[1] / (2π) *
           (log((1 + argument) / argument) + correction)
end

:Theethayi2007
