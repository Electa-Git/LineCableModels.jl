function routes(identifier::Val{:Wait1972a})
    return (
        self=FormulaMethod(identifier,earth_potential_coefficient,Val(:self)),
        mutual=FormulaMethod(identifier,earth_potential_coefficient,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:Wait1972a})=(air=_full,earth=_full,permeability=vacuum_permeability)
propagation(::Val{:Wait1972a})=Val(:zero)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite thin circular conductor of radius ``a`` at height ``h`` parallel to a planar interface. |
| Calculated quantities | Generalized full-wave shunt admittance and its Carson qTEM reduction |
| Earth structure | Two homogeneous half-spaces. |
| Model and approximation | Equation (25) retains the solved modal propagation constant and the interface spectral term; the qTEM reduction follows from equations (30)–(36). |
| Main source | James R. Wait (1972) |
| Citation key(s) | Primary: `:Wait1972a`; earlier result: `:Kikuchi1956`; validity analysis: `:Pogorzelski1977`; excitation witness: `:Kuester1978` |
| Evidence status | Original PDF equations (23)–(36) checked; the shunt-admittance definition is explicit |

**Implemented scope.** Only the qTEM self reduction (32) is selected.
The full-wave admittance and its implicit modal root problem are not
represented by this frequency-only coefficient.

**Expression.** Return the complete external potential coefficient,

```math
P_{ii}=\\frac{\\ln(2h/a)}{2\\pi\\varepsilon_1},\\qquad
Y_{ii}=\\frac{j\\omega}{P_{ii}}.
```

The earth family returns ``P``, not the scalar reciprocal ``Y``.
The engine therefore applies its usual full potential-matrix inversion.
This single-wire expression has no independently supplied mutual route;
use a multiconductor potential formula for more than one external wire.
Unlike the zero-correction IdealGround selector, this includes the
external geometric logarithm.

**Reference.** [Wait1972a](@cite), equation (32), p. 678.
"""
description(::Formula{:Wait1972a})="Wait overhead-wire qTEM self potential (1972)"
propagation_constant(::Val{:Wait1972a},jω,permeability,permittivity)=(Γ=zero(jω),squared=zero(jω))
function (formula::Formula{:Wait1972a})(rho,epsilon,mu,jω,Γ,segments=nothing)
    return _homogeneous_functor(Val(:Wait1972a),formula,rho,epsilon,mu,jω,Γ,segments)
end
function earth_potential_coefficient(::Val{:Wait1972a},::Val{:mutual},functor,pair)
    throw(ArgumentError("Wait1972a qTEM reduction supplies only an overhead self coefficient"))
end
function earth_potential_coefficient(::Val{:Wait1972a},::Val{:self},functor,pair)
    _require(pair,Val(:overhead))
    pair.row==pair.column || throw(ArgumentError("Wait1972a requires self geometry"))
    state=functor.state; geometry=_geometry(pair)
    return complex(log(geometry.H/geometry.y_ij)/(2*(one(geometry.H)*π)*state.epsilon[1]))
end

:Wait1972a
