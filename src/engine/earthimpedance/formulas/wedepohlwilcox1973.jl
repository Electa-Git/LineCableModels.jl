function assumptions(::Val{:WedepohlWilcox1973})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Low-frequency underground expansion for conductive,
nonmagnetic earth.

**Expression.**

```math
Z_{e,ii}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[-\\ln
\\left(\\frac{e_c\\gamma_1r_i}{2}\\right)+\\frac12-
\\frac43\\gamma_1h_i\\right],
```

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[-\\ln
\\left(\\frac{e_c\\gamma_1d_{ij}}{2}\\right)+\\frac12-
\\frac23\\gamma_1H\\right],\\qquad e_c=1.7811.
```

Here ``d_{ij}=\\sqrt{x_{ij}^2+(h_i-h_j)^2}`` is the distance between
cable axes and ``H=h_i+h_j`` is the sum of burial depths, all in metres.

**Reference.** L. M. Wedepohl and D. J. Wilcox, “Transient Analysis of
Underground Power-Transmission Systems: System-Model and Wave-Propagation
Characteristics,” *Proceedings of the IEE*, 120, 253–260, 1973, Eqs. (7)–(8).
DOI: 10.1049/piee.1973.0056.
"""
function description(::Formula{:WedepohlWilcox1973})
    "Wedepohl-Wilcox low-frequency underground approximation (1973)"
end

function Γ(
        ::Val{:WedepohlWilcox1973}, jω, materials, layers
)
    zero(jω)
end

raw"""
Evaluate the Wedepohl-Wilcox low-frequency underground terms:

```math
Z_{e,ii}=\frac{j\omega\mu_0}{2\pi}
\left[-\ln\left(\frac{e_c\gamma_1r_i}{2}\right)+\frac12
-\frac43\gamma_1h_i\right],
```

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}
\left[-\ln\left(\frac{e_c\gamma_1d_{ij}}{2}\right)+\frac12
-\frac23\gamma_1(h_i+h_j)\right],\qquad e_c=1.7811.
```

The mutual distance is ``d_{ij}=\sqrt{x_{ij}^2+(h_i-h_j)^2}``, not
the horizontal projection ``x_{ij}``. Distinct cable axes can therefore
be vertically aligned. The shared pair validator rejects coincident axes.
"""
function earth_impedance(
        ::Val{:WedepohlWilcox1973}, ::Val{:self}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    e_c = oftype(geometry.h_i, 1.7811)
    bracket = -log(e_c * state.gamma[2] * geometry.y_ij / 2) +
              one(geometry.h_i) / 2 -
              (4 * one(geometry.h_i) / 3) * state.gamma[2] * geometry.h_i
    return state.jω * state.mu[1] / (2π) * bracket
end

function earth_impedance(
        ::Val{:WedepohlWilcox1973}, ::Val{:mutual}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    e_c = oftype(geometry.H, 1.7811)
    bracket = -log(e_c * state.gamma[2] * geometry.d_ij / 2) +
              one(geometry.H) / 2 -
              (2 * one(geometry.H) / 3) * state.gamma[2] * geometry.H
    return state.jω * state.mu[1] / (2π) * bracket
end

Formulation(::LineCableModelsCoaxial, selected::Formula{:WedepohlWilcox1973}) = selected

function hooks(::FormulaMethod{:WedepohlWilcox1973, typeof(earth_impedance),
        A}) where {A <: Tuple{Val{:self}, Val{2}, Val{2}}}
    return (configurable = (:Γ, :earth, :permeability, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:WedepohlWilcox1973), Γ),
            air = FormulaMethod(Val(:lossless), propagation),
            earth = FormulaMethod(Val(:conductive), propagation),
            permeability = vacuum_permeability,
            contribution = nothing))
end

function computation_options(::FormulaMethod{:WedepohlWilcox1973, typeof(earth_impedance),
        A}) where {A <: Tuple{Val{:self}, Val{2}, Val{2}}}
    (;)
end

function hooks(::FormulaMethod{:WedepohlWilcox1973, typeof(earth_impedance),
        A}) where {A <: Tuple{Val{:mutual}, Val{2}, Val{2}}}
    return (configurable = (:Γ, :earth, :permeability, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:WedepohlWilcox1973), Γ),
            air = FormulaMethod(Val(:lossless), propagation),
            earth = FormulaMethod(Val(:conductive), propagation),
            permeability = vacuum_permeability,
            contribution = nothing))
end

function computation_options(::FormulaMethod{:WedepohlWilcox1973, typeof(earth_impedance),
        A}) where {A <: Tuple{Val{:mutual}, Val{2}, Val{2}}}
    (;)
end

function validate(binding::FormulaMethod{:WedepohlWilcox1973, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:default})
    binding
end

:WedepohlWilcox1973
