routes(::Val{:IdealGround}) = (
    self = ideal_zero,
    mutual = ideal_zero,
    Γ = ideal_gamma
)

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

function ideal_gamma(jω, permeability, permittivity)
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

ideal_zero(functor, pair) = zero(functor.state.jω)

@inline function (functor::Functor{:IdealGround})(::Val{:self}, pair)
    return functor.routes.self(functor, pair)
end

@inline function (functor::Functor{:IdealGround})(::Val{:mutual}, pair)
    return functor.routes.mutual(functor, pair)
end

:IdealGround
