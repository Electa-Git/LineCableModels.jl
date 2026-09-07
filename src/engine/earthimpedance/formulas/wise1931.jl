function routes(identifier::Val{:Wise1931})
    return (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

assumptions(::Val{:Wise1931}) = (
    air = _conductive, earth = _conductive, permeability = _material
)
propagation(::Val{:Wise1931}) = Val(:zero)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Thin round overhead wires; self expression uses conductor radius ``a`` and mutual expression uses geometric distances ``\\rho'`` and ``\\rho''``. |
| Calculated quantities | P.u.l. self and mutual series impedances of parallel overhead ground-return circuits |
| Earth structure | Homogeneous conducting half-space below air. |
| Model and approximation | None in the displayed integral beyond the inherited Carson model. Equations (23)–(24) retain the integral; Wise's separate asymptotic and convergent evaluator series are excluded. |
| Main source | W. Howard Wise (1931), extending Carson's 1926 ground-return derivation to non-unit earth permeability |
| Citation key(s) | `:Wise1931` |
| Evidence status | Original publication page images checked |

**Expression.** The SI evaluator retains the permeability-weighted denominator
of the 1931 integral, with both displacement-current terms omitted.
It shares the numerical integral with Wise1934 under those constitutive
restrictions. The nonmagnetic case reduces to Carson1926.

**Reference.** [Wise1931](@cite), equations (23)–(24).
"""
description(::Formula{:Wise1931}) = "Wise permeable-earth overhead impedance (1931)"

function propagation_constant(::Val{:Wise1931}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Wise1931})(rho, epsilon, mu, jω, Γ, segments = nothing)
    isinf(first(rho)) || throw(ArgumentError(":Wise1931 requires nonconducting air"))
    return _homogeneous_functor(Val(:Wise1931), formula, rho, epsilon, mu, jω, Γ, segments)
end

function earth_impedance(::Val{:Wise1931}, ::Val{:mutual}, functor, pair)
    return _complex_result(functor.state.jω,
        earth_impedance(Val(:Wise1934), Val(:mutual), functor, pair))
end

:Wise1931
