function routes(identifier::Val{:Theethayi2007})
    (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        impedance = FormulaMethod(identifier, earth_impedance, Val(:support)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Theethayi2007})
    (
        air = _vacuum,
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
| Family | External admittance |
| Geometry | Buried self or horizontally separated mutual conductor pair. |
| Calculated quantities | Potential coefficients from the selected logarithmic-exponential impedance |
| Earth structure | Homogeneous conducting nonmagnetic earth below air. |
| Model and approximation | Vance scalar conversion with the Theethayi impedance approximation; matrix extension follows the IET book. |
| Main source | N. Theethayi, R. Thottappillil, M. Paolone, C. A. Nucci, and F. Rachidi (2007) |
| Citation key(s) | `:Theethayi2007`; scalar attribution: `:Vance1978`; matrix extension: `:Ametani2021` |
| Evidence status | Scalar source equations (9)–(10) and IET book Section 2.6.3 checked. |

**Description.** Underground earth-admittance relation derived from the
matching logarithmic-exponential impedance model.

**Expression.** At matrix level,

```math
\\mathbf Y_e=\\gamma_1^2\\mathbf Z_e^{-1},\\qquad
P_{e,ij}=\\frac{j\\omega Z_{e,ij}}{\\gamma_1^2}.
```

The inverse is taken only after assembling ``\\mathbf P_e``; the package does
not divide admittance element by impedance element.

**Reference.** N. Theethayi, R. Thottappillil, M. Paolone, C. A. Nucci, and
F. Rachidi, “External Impedance and Admittance of Buried Horizontal Wires for
Transient Studies Using Transmission Line Analysis,” *IEEE Transactions on
Dielectrics and Electrical Insulation*, 14(3), 751–761, 2007.
DOI: 10.1109/TDEI.2007.369540.
"""
function description(::Formula{:Theethayi2007})
    "Theethayi et al. earth admittance from impedance (2007)"
end

function propagation_constant(
        ::Val{:Theethayi2007}, jω, permeability, permittivity
)
    (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Theethayi2007})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Theethayi2007), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Theethayi et al. earth admittance relation

```math
Y_e=\frac{\gamma_1^2}{Z_e}.
```

The default `impedance` leaf is the Theethayi logarithmic-exponential
approximation from the same corpus. Because the engine assembles the full
potential-coefficient matrix before inversion, each route returns
``P_{e,ij}=j\omega Z_{e,ij}/\gamma_1^2``. Consequently the assembled matrices
satisfy ``\mathbf Y_e=\gamma_1^2\mathbf Z_e^{-1}``; no elementwise
admittance quotient is taken. Self evaluation uses the same pair formula with
horizontal separation set to the conductor radius.
"""
function earth_potential_coefficient(
        ::Val{:Theethayi2007}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    pair.row == pair.column || _require_horizontal_separation(pair)
    state = functor.state
    Z_e = functor.routes.impedance(state, pair)
    return impedance_potential_coefficient(state, Z_e)
end

function earth_impedance(
        ::Val{:Theethayi2007}, ::Val{:support}, state, pair
)
    proxy=(;state)
    return EarthImpedance.earth_impedance(Val(:Theethayi2007),Val(:mutual),proxy,pair)
end

:Theethayi2007
