function assumptions(::Val{:gary1976})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Complex-depth logarithmic approximation for overhead
conductors.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\ln\\frac{S_{ij}}{d_{ij}},\\qquad
S_{ij}=\\sqrt{(H+2h_e)^2+y_{ij}^2},\\qquad
h_e=(j\\omega\\mu_0\\sigma_1)^{-1/2}.
```

**Reference.** C. Gary, “Approche complète de la propagation multifilaire en
haute fréquence par utilisation des matrices complexes,” *EDF Bulletin de la
Direction des Études et Recherches*, série B, 1976; formula as reproduced in
Ametani et al., IET, 2021.
"""
description(::Type{<:Formula{:gary1976}}; compact::Bool=false) = compact ? "Gary" : "Gary complex-depth approximation (1976)"


raw"""
Evaluate Gary's complex-depth overhead earth-return impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\ln\frac{S_{ij}}{d_{ij}},\qquad
S_{ij}=\sqrt{(h_i+h_j+2h_e)^2+y_{ij}^2},\qquad
h_e=\frac{1}{\sqrt{j\omega\mu_0\sigma_1}}.
```

The direct distance is
``d_{ij}=\sqrt{(h_i-h_j)^2+y_{ij}^2}``. For self terms, `pair`
supplies the outer radius as ``y_{ii}=r_i``.
"""
function earth_impedance(
        ::Formula{:gary1976}, ::Union{Val{:self}, Val{:mutual}}, ::Val{1}, ::Val{1},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    h_e = inv(state.gamma[2])
    S_ij = sqrt((geometry.H + 2h_e)^2 + geometry.y_ij^2)
    return state.jω * state.mu[1] / (2π) * log(S_ij / geometry.d_ij)
end



function formulation_options(::FormulaMethod{<:Formula{:gary1976}, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{1}, Val{1}}}
    return FormulationOptions((;))
end

function validate(binding::FormulaMethod{<:Formula{:gary1976}, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:bottommost})
    binding
end

"""
$(TYPEDSIGNATURES)

Evaluate this formulation's medium state: absolute permeability \\[H/m\\]
and transverse propagation constant \\[1/m\\]. Air retains its prescribed
permeability; soil follows the selected source's magnetic approximation.
"""
function constitutive(::Formula{:gary1976}, ::Val{:air}, jω, μ, σ, ε)
    return (mu=μ, gamma=propagation(Val(:lossless), jω, μ, σ, ε))
end

function constitutive(::Formula{:gary1976}, ::Val{:earth}, jω, μ, σ, ε)
    permeability = vacuum_permeability(μ)
    return (mu=permeability, gamma=propagation(Val(:conductive), jω, permeability, σ, ε))
end

:gary1976
