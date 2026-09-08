function assumptions(::Val{:IdealGround})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

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

Γ(::Val{:IdealGround}, jω, materials, layers) = zero(jω)

function earth_potential_coefficient(
        ::Val{:IdealGround}, ::Union{Val{:self}, Val{:mutual}}, ::Val{1}, ::Val{1},
        functor, pair, workspace)
    return zero(functor.state.jω)
end

function earth_potential_coefficient(
        ::Val{:IdealGround}, ::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace)
    return zero(functor.state.jω)
end

function earth_potential_coefficient(
        ::Val{:IdealGround}, ::Val{:mutual}, ::Val{1}, ::Val{2},
        functor, pair, workspace)
    return zero(functor.state.jω)
end

function earth_potential_coefficient(
        ::Val{:IdealGround}, ::Val{:mutual}, ::Val{2}, ::Val{1},
        functor, pair, workspace)
    return zero(functor.state.jω)
end

Formulation(::LineCableModelsCoaxial, selected::Formula{:IdealGround}) = selected

function hooks(::FormulaMethod{:IdealGround, typeof(earth_potential_coefficient)})
    return (configurable = (:Γ, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:IdealGround), Γ),
            air = FormulaMethod(Val(:full), propagation),
            earth = FormulaMethod(Val(:full), propagation),
            permeability = identity,
            contribution = nothing))
end

function computation_options(::FormulaMethod{
        :IdealGround, typeof(earth_potential_coefficient)})
    (;)
end

function validate(
        binding::FormulaMethod{:IdealGround, typeof(earth_potential_coefficient)},
        ::EquivalentHomogeneous.Formula{:default})
    binding
end

:IdealGround
