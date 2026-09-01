function routes(identifier::Val{:IdealGround})
    return (
        self = formula_method(identifier, earth_potential_coefficient, Val(:self)),
        mutual = formula_method(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = formula_method(identifier, propagation_constant)
    )
end

assumptions(::Val{:IdealGround}) = (;)
propagation(::Val{:IdealGround}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Limiting reference in which the ground is an ideal
equipotential conductor and contributes no earth potential coefficient.

**Expression.**

```math
P_{e,ij}=0,\\qquad \\Gamma=0.
```

Only the non-earth electrostatic or insulation terms remain in the assembled
potential-coefficient matrix.

**Reference.** Ideal-conductor boundary condition; no empirical literature
fit is introduced by this reference case.
"""
description(::Formula{:IdealGround}) = "Ideal ground reference"

function propagation_constant(
        ::Val{:IdealGround}, jω, permeability, permittivity
)
    value = zero(jω)
    return (Γ = value, squared = value)
end

function (formula::Formula{:IdealGround})(
        resistivity::AbstractVector{T},
        permittivity::AbstractVector{T},
        permeability::AbstractVector{T},
        jω::Complex{T},
        Γ,
        segments = nothing
) where {T <: Real}
    longitudinal = _longitudinal(
        formula,
        Γ,
        jω,
        first(permeability),
        first(permittivity)
    )
    state = (; formula, jω, Γ = longitudinal.Γ)
    return Functor{:IdealGround, typeof(formula.routes), typeof(state)}(
        formula.routes,
        state
    )
end

function earth_potential_coefficient(
        ::Val{:IdealGround}, ::Val{:mutual}, functor, pair
)
    return zero(functor.state.jω)
end

@inline function (functor::Functor{:IdealGround})(::Val{:self}, pair)
    return functor.routes.self(functor, pair)
end

@inline function (functor::Functor{:IdealGround})(::Val{:mutual}, pair)
    return functor.routes.mutual(functor, pair)
end

:IdealGround
