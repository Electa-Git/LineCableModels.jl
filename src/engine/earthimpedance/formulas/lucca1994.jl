function assumptions(::Val{:Lucca1994})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Homogeneous-earth mixed-pair model with a corrected
complex-depth approximation for mixed overhead-underground coupling.

**Expression.** Its distinctive mixed term is

```math
Z_{e,ij}^{01}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{S}{D}-\\frac23\\left(\\frac{h_e}{S^2}\\right)^3
H(H^2-3y_{ij}^2)\\right],
```

```math
h_e=(j\\omega\\mu_0\\sigma_g)^{-1/2},\\quad
H=h_a+h_g+2h_e,\\quad S=\\sqrt{H^2+y_{ij}^2},\\quad
D=\\sqrt{(h_a+h_g)^2+y_{ij}^2}.
```

**Reference.** G. Lucca, “Mutual Impedance Between an Overhead and a Buried
Line with Earth Return,” *9th International Conference on Electromagnetic
Compatibility*, 1994. DOI: 10.1049/cp:19940679.
"""
description(::Formula{:Lucca1994}) = "Lucca mixed-pair homogeneous-earth impedance (1994)"

function Γ(::Val{:Lucca1994}, jω, materials, layers)
    return zero(jω)
end

raw"""
Evaluate Lucca's approximation for the mutual impedance between one overhead
and one buried conductor:

```math
Z_{e,ij}^{01}=\frac{j\omega\mu_0}{2\pi}\left[
\ln\frac{S}{D}-\frac23\left(\frac{h_e}{S^2}\right)^3
H(H^2-3y_{ij}^2)\right],
```

```math
h_e=\frac1{\sqrt{j\omega\mu_0\sigma_g}},\qquad
H=h_a+h_g+2h_e,\qquad
S=\sqrt{H^2+y_{ij}^2},\qquad
D=\sqrt{(h_a+h_g)^2+y_{ij}^2}.
```

Only the published mixed interaction is registered.

# Reference

G. Lucca, "Mutual impedance between an overhead and a buried line with earth
return," *9th International Conference on Electromagnetic Compatibility*, 1994.
DOI: 10.1049/cp:19940679.
"""
function earth_impedance(
        ::Val{:Lucca1994}, ::Val{:mutual}, ::Val{1}, ::Val{2},
        functor, pair, workspace
)
    state = functor.state
    air = pair.layers[1] == 1 ? 1 : 2
    earth = air == 1 ? 2 : 1
    h_a = abs(pair.heights[air])
    h_g = abs(pair.heights[earth])
    h_e = inv(state.gamma[2])
    H = h_a + h_g + 2h_e
    S_squared = H^2 + pair.separation^2
    S = sqrt(S_squared)
    D = hypot(pair.separation, h_a + h_g)
    correction = (2 * one(h_a) / 3) * (h_e / S_squared)^3 *
                 H * (H^2 - 3 * pair.separation^2)
    πT = one(h_a) * π
    return state.jω * state.mu[1] / (2πT) * (log(S / D) - correction)
end

function earth_impedance(
        ::Val{:Lucca1994}, ::Val{:mutual}, ::Val{2}, ::Val{1},
        functor, pair, workspace
)
    state = functor.state
    air = pair.layers[1] == 1 ? 1 : 2
    earth = air == 1 ? 2 : 1
    h_a = abs(pair.heights[air])
    h_g = abs(pair.heights[earth])
    h_e = inv(state.gamma[2])
    H = h_a + h_g + 2h_e
    S_squared = H^2 + pair.separation^2
    S = sqrt(S_squared)
    D = hypot(pair.separation, h_a + h_g)
    correction = (2 * one(h_a) / 3) * (h_e / S_squared)^3 *
                 H * (H^2 - 3 * pair.separation^2)
    πT = one(h_a) * π
    return state.jω * state.mu[1] / (2πT) * (log(S / D) - correction)
end

Formulation(::LineCableModelsCoaxial, selected::Formula{:Lucca1994}) = selected

function hooks(::FormulaMethod{:Lucca1994, typeof(earth_impedance),
        A}) where {A <: Tuple{Val{:mutual}, Val{1}, Val{2}}}
    return (configurable = (:Γ, :earth, :permeability, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:Lucca1994), Γ),
            air = FormulaMethod(Val(:lossless), propagation),
            earth = FormulaMethod(Val(:conductive), propagation),
            permeability = vacuum_permeability,
            contribution = nothing))
end

function computation_options(::FormulaMethod{:Lucca1994, typeof(earth_impedance),
        A}) where {A <: Tuple{Val{:mutual}, Val{1}, Val{2}}}
    (;)
end

function hooks(::FormulaMethod{:Lucca1994, typeof(earth_impedance),
        A}) where {A <: Tuple{Val{:mutual}, Val{2}, Val{1}}}
    return (configurable = (:Γ, :earth, :permeability, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:Lucca1994), Γ),
            air = FormulaMethod(Val(:lossless), propagation),
            earth = FormulaMethod(Val(:conductive), propagation),
            permeability = vacuum_permeability,
            contribution = nothing))
end

function computation_options(::FormulaMethod{:Lucca1994, typeof(earth_impedance),
        A}) where {A <: Tuple{Val{:mutual}, Val{2}, Val{1}}}
    (;)
end

function validate(binding::FormulaMethod{:Lucca1994, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:default})
    binding
end

:Lucca1994
