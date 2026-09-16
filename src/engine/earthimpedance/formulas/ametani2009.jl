function assumptions(::Val{:ametani2009})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Homogeneous-earth approximation for mixed overhead-underground pairs.

**Expression.** Its distinctive mixed term is

```math
Z_{e,ij}^{01}=\\frac{j\\omega\\mu_0}{2\\pi}e^{-h_g/h_e}\\ln\\frac{S}{D},
\\quad h_e=(j\\omega\\mu_0\\sigma_g)^{-1/2},
```

```math
S=\\sqrt{(h_a+h_g+2h_e)^2+y_{ij}^2},\\qquad
D=\\sqrt{(h_a+h_g)^2+y_{ij}^2}.
```

**Reference.** A. Ametani, T. Yoneda, Y. Baba, and N. Nagaoka, “An
Investigation of Earth-Return Impedance Between Overhead and Underground
Conductors and Its Approximation,” *IEEE Transactions on Electromagnetic
Compatibility*, 51, 860–867, 2009.
DOI: 10.1109/TEMC.2009.2019953.
"""
function description(::Type{<:Formula{:ametani2009}}; compact::Bool=false)
    compact ? "Ametani" : "Ametani mixed-pair homogeneous-earth impedance (2009)"
end


raw"""
Evaluate Ametani's approximation for the mutual impedance between one
overhead and one buried conductor:

```math
Z_{e,ij}^{01}=\frac{j\omega\mu_0}{2\pi}
e^{-h_g/h_e}\ln\frac{S}{D},
```

```math
h_e=\frac1{\sqrt{j\omega\mu_0\sigma_g}},\qquad
S=\sqrt{(h_a+h_g+2h_e)^2+y_{ij}^2},\qquad
D=\sqrt{(h_a+h_g)^2+y_{ij}^2}.
```

Here ``h_a`` and ``h_g`` are positive height and burial-depth magnitudes.
Only the published mixed interaction is registered.

# Reference

A. Ametani, "An investigation of earth-return impedance between overhead
and underground conductors and its approximation," *IEEE Transactions on
Electromagnetic Compatibility*, vol. 51, pp. 860-867, 2009.
DOI: 10.1109/TEMC.2009.2019953.
"""
function earth_impedance(
        ::Formula{:ametani2009}, ::Val{:mutual}, ::Val{1}, ::Val{2},
        functor, pair, workspace
)
    state = functor.state
    air = pair.layers[1] == 1 ? 1 : 2
    earth = air == 1 ? 2 : 1
    h_a = abs(pair.heights[air])
    h_g = abs(pair.heights[earth])
    h_e = inv(state.gamma[2])
    D = hypot(pair.separation, h_a + h_g)
    S = sqrt((h_a + h_g + 2h_e)^2 + pair.separation^2)
    πT = one(h_a) * π
    return state.jω * state.mu[1] / (2πT) *
           exp(-h_g / h_e) * log(S / D)
end

function earth_impedance(
        ::Formula{:ametani2009}, ::Val{:mutual}, ::Val{2}, ::Val{1},
        functor, pair, workspace
)
    state = functor.state
    air = pair.layers[1] == 1 ? 1 : 2
    earth = air == 1 ? 2 : 1
    h_a = abs(pair.heights[air])
    h_g = abs(pair.heights[earth])
    h_e = inv(state.gamma[2])
    D = hypot(pair.separation, h_a + h_g)
    S = sqrt((h_a + h_g + 2h_e)^2 + pair.separation^2)
    πT = one(h_a) * π
    return state.jω * state.mu[1] / (2πT) *
           exp(-h_g / h_e) * log(S / D)
end



function formulation_options(::FormulaMethod{<:Formula{:ametani2009}, typeof(earth_impedance),
        A}) where {A <: Tuple{Val{:mutual}, Val{1}, Val{2}}}
    return FormulationOptions((;))
end


function formulation_options(::FormulaMethod{<:Formula{:ametani2009}, typeof(earth_impedance),
        A}) where {A <: Tuple{Val{:mutual}, Val{2}, Val{1}}}
    return FormulationOptions((;))
end

function validate(binding::FormulaMethod{<:Formula{:ametani2009}, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:bottommost})
    binding
end

"""
$(TYPEDSIGNATURES)

Evaluate this formulation's medium state: absolute permeability \\[H/m\\]
and transverse propagation constant \\[1/m\\]. Air retains its prescribed
permeability; soil follows the selected source's magnetic approximation.
"""
function constitutive(::Formula{:ametani2009}, ::Val{:air}, jω, μ, σ, ε)
    return (mu=μ, gamma=propagation(Val(:lossless), jω, μ, σ, ε))
end

function constitutive(::Formula{:ametani2009}, ::Val{:earth}, jω, μ, σ, ε)
    permeability = vacuum_permeability(μ)
    return (mu=permeability, gamma=propagation(Val(:conductive), jω, permeability, σ, ε))
end

:ametani2009
