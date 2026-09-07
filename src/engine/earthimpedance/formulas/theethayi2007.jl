function routes(identifier::Val{:Theethayi2007})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
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

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Circular bare conductor radius ``a`` or concentric insulated wire with outer radius ``b``; ``R_{ab}=a`` for bare, ``R_{ab}=b`` for insulated. Burial depth is ``d``. |
| Calculated quantities | Buried bare/insulated-wire earth impedance; source-prescribed parallel-wire mutual substitution |
| Earth structure | Homogeneous conducting dielectric half-space below air; empirical correction represents burial-depth/interface effects. |
| Model and approximation | Empirical modification of the infinite-earth logarithmic approximation attributed to Petrache et al. (2005) in (8). A depth-dependent exponential term is added in (9); the inspected paper does not derive an expansion parameter, retained order, or discarded remainder for it. Its Sunde/Wait comparisons are presented as numerical evidence within the stated geometries. |
| Main source | Nelson Theethayi's 2005 doctoral thesis, explicitly cited as reference [14] by the inspected 2007 paper; the record year identifies this later published witness |
| Citation key(s) | Primary publication: `:Theethayi2007`; thesis source: `:Theethayi2005` |
| Evidence status | Author's later explicit restatement and original thesis equation (7.10) checked against page images; the relevant thesis sections contain no additional distinct formula in the requested families. |

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
    return state.jω * state.mu[1] / (2*(one(geometry.H)*π)) *
           (log((1 + argument) / argument) + correction)
end

:Theethayi2007
