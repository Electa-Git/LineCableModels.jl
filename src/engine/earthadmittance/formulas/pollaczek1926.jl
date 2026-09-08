function assumptions(::Val{:Pollaczek1926})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Classical homogeneous-earth underground potential coefficient.

**Expression.**

```math
P_{e,ij}^{11}=\\frac{j\\omega}{2\\pi(\\sigma_1+j\\omega\\varepsilon_1)}
[K_0(\\gamma_0d_{ij})-K_0(\\gamma_0D_{ij})].
```

**Reference.** F. Pollaczek, “Über das Feld einer unendlich langen
wechselstromdurchflossenen Einfachleitung,” *Elektrische Nachrichtentechnik*,
3, 339–360, 1926; potential-coefficient transcription follows Ametani et al.,
IET, 2021.
"""
function description(::Formula{:Pollaczek1926})
    "Pollaczek underground potential coefficients (1926)"
end

function Γ(
        ::Val{:Pollaczek1926}, jω, materials, layers
)
    zero(jω)
end

raw"""
Evaluate Pollaczek's classical homogeneous-earth underground potential
coefficient:

```math
P_{e,ij}^{11}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
\left[K_0(\gamma_0d_{ij})-K_0(\gamma_0D_{ij})\right],
\qquad \gamma_0=j\omega\sqrt{\mu_0\varepsilon_0}.
```

Earth conductivity remains in the potential-coefficient prefactor.
"""
function earth_potential_coefficient(
        ::Val{:Pollaczek1926}, ::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    gamma_0 = state.gamma[2]
    direct = special_besselk(0, gamma_0 * geometry.d_ij) -
             special_besselk(0, gamma_0 * geometry.D_ij)
    kappa_1 = state.sigma[2] + state.jω * state.epsilon[2]
    return state.jω / (2π * kappa_1) * direct
end

Formulation(::LineCableModelsCoaxial, selected::Formula{:Pollaczek1926}) = selected

function hooks(::FormulaMethod{:Pollaczek1926, typeof(earth_potential_coefficient),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{2}, Val{2}}}
    return (configurable = (:Γ, :earth, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:Pollaczek1926), Γ),
            air = FormulaMethod(Val(:vacuum), propagation),
            earth = FormulaMethod(Val(:vacuum), propagation),
            permeability = vacuum_permeability,
            contribution = nothing))
end

function computation_options(::FormulaMethod{
        :Pollaczek1926, typeof(earth_potential_coefficient),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{2}, Val{2}}}
    (;)
end

function validate(
        binding::FormulaMethod{:Pollaczek1926, typeof(earth_potential_coefficient)},
        ::EquivalentHomogeneous.Formula{:default})
    binding
end

:Pollaczek1926
