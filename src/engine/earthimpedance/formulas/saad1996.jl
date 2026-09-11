function assumptions(::Val{:Saad1996})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Closed-form underground approximation combining the
direct cylindrical term and an interface correction.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
K_0(\\gamma_1R_{ab})+
\\frac{2e^{-H\\gamma_1}}{4+\\gamma_1^2R_{ab}^2}\\right].
```

**Reference.** O. Saad, G. Gaba, and M. Giroux, “A Closed-Form Approximation
for Ground Return Impedance of Underground Cables,” *IEEE Transactions on
Power Delivery*, 11(3), 1536–1545, 1996.
"""
description(::Formula{:Saad1996}) = "Saad underground closed form (1996)"

function Γ(::Val{:Saad1996}, jω, materials, layers)
    return zero(jω)
end

raw"""
Evaluate the Saad et al. underground approximation:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}
\left[K_0(\gamma_1R_{ab})+
\frac{2e^{-(h_i+h_j)\gamma_1}}{4+\gamma_1^2R_{ab}^2}\right].
```

For a self term ``R_{ab}=r_i``; for a mutual term it is the horizontal
center separation supplied by `pair`.

# Reference

O. Saad, G. Gaba, and M. Giroux, "A closed-form approximation for ground
return impedance of underground cables," *IEEE Transactions on Power
Delivery*, vol. 11, no. 3, pp. 1536-1545, 1996.
"""
function earth_impedance(
        ::Val{:Saad1996}, ::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    gamma = state.gamma[2]
    radius = geometry.y_ij
    correction = 2exp(-geometry.H * gamma) / (4 + gamma^2 * radius^2)
    direct = oftype(
        state.jω, special_besselk(0, gamma * radius)
    )
    πT = one(radius) * π
    return state.jω * state.mu[1] / (2πT) * (direct + correction)
end

function validate(pair::EarthPair, ::FormulaMethod{:Saad1996, typeof(earth_impedance)})
    pair.row != pair.column && iszero(pair.separation) &&
        throw(DomainError(
            pair.separation, ":Saad1996 mutual closed form requires nonzero horizontal cable separation"))
    return pair
end

Formulation(::LineCableModelsCoaxial, selected::Formula{:Saad1996}) = selected

function hooks(::FormulaMethod{:Saad1996, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{2}, Val{2}}}
    return (configurable = (:Γ, :earth, :permeability, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:Saad1996), Γ),
            air = FormulaMethod(Val(:lossless), propagation),
            earth = FormulaMethod(Val(:conductive), propagation),
            permeability = vacuum_permeability,
            contribution = nothing))
end

function computation_options(::FormulaMethod{:Saad1996, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{2}, Val{2}}}
    (;)
end

function validate(binding::FormulaMethod{:Saad1996, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:default})
    binding
end

:Saad1996
