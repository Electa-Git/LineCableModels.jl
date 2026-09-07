function routes(identifier::Val{:Vance1978})
    return (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        impedance = FormulaMethod(identifier, earth_impedance, Val(:support)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

assumptions(::Val{:Vance1978}) = (
    air = _vacuum, earth = _full, permeability = vacuum_permeability
)
propagation(::Val{:Vance1978}) = Val(:zero)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Bare circular radius ``a`` or concentric insulated wire with outer radius ``b`` at depth ``d``. Ground impedance selects ``R_{ab}=a`` or ``b`` respectively. |
| Calculated quantities | Scalar ground admittance of a bare wire and of an insulated wire; source-supplied insulation series assembly |
| Earth structure | Homogeneous conducting dielectric half-space below air. |
| Model and approximation | A source-attributed scalar relation in the selected TL description. When ``Z_g`` is approximated by (9), its empirical approximation carries into ``Y_g``. |
| Main source | E. F. Vance (1978), as reproduced by later cable-parameter publications |
| Citation key(s) | Attributed source: `:Vance1978`; equation witnesses: `:Theethayi2007`, `:Zhang2017`, `:Guneri2018` |
| Evidence status | The relation is independently printed in three accessible publications |

**Expression.** The source supplies the scalar relation

```math
Y_g=\\frac{\\gamma_g^2}{Z_g},\\qquad
\\gamma_g^2=j\\omega\\mu_0(\\sigma_g+j\\omega\\varepsilon_g).
```

For the engine's single-exterior-unit matrix, the returned coefficient is
``P_g=j\\omega/Y_g=j\\omega Z_g/\\gamma_g^2``. The default impedance
dependency is Vance1978's cylindrical exterior impedance. The
`impedance` route may be replaced by a compatible scalar earth-impedance
evaluator; it receives `(functor, pair)`.

Only underground self evaluation is supported. A second external unit requires
mutual potential coefficients not established by this scalar source and is
rejected. Internal coaxial dielectric coefficients are assembled separately;
neither conductor internal impedance nor insulation series impedance belongs
in the input ``Z_g``.

**Reference.** [Vance1978](@cite), as explicitly reproduced in
[Theethayi2007](@cite), equations (10a)–(10b), and the other witnesses listed
in the source table.
"""
description(::Formula{:Vance1978}) = "Vance scalar buried-wire ground admittance (1978)"

function propagation_constant(::Val{:Vance1978}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Vance1978})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Vance1978), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

function earth_potential_coefficient(::Val{:Vance1978}, ::Val{:self}, functor, pair)
    _require(pair, Val(:underground))
    pair.row == pair.column || throw(ArgumentError(
        ":Vance1978 supplies a scalar ground admittance, not mutual potential coefficients"
    ))
    state = functor.state
    impedance = functor.routes.impedance(functor, pair)
    return impedance_potential_coefficient(state, impedance)
end

function earth_potential_coefficient(::Val{:Vance1978}, ::Val{:mutual}, functor, pair)
    throw(ArgumentError(
        ":Vance1978 supplies a scalar ground admittance; multiple external units require a matrix-capable formula"
    ))
end

function earth_impedance(::Val{:Vance1978}, ::Val{:support}, functor, pair)
    return EarthImpedance.earth_impedance(Val(:Vance1978), Val(:self), functor, pair)
end

:Vance1978
