function routes(::Val{:Ametani2021})
    (
        self = ametani_space2021,
        mutual = ametani_space2021,
        Γ = ametani_space2021_gamma
    )
end

function assumptions(::Val{:Ametani2021})
    (
        air = _vacuum,
        earth = _vacuum,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Ametani2021}) = Val(:zero)
function description(::Formula{:Ametani2021})
    "Ametani et al. classical overhead space potential coefficient (2021)"
end

ametani_space2021_gamma(jω, permeability, permittivity) = (Γ = zero(jω), squared = zero(jω))

function (formula::Formula{:Ametani2021})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Ametani2021), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the classical overhead space potential coefficient:

```math
P_{0,ij}=\frac{1}{2\pi\varepsilon_0}\ln\frac{D_{ij}}{d_{ij}},
\qquad \mathbf Y_0=j\omega\mathbf P_0^{-1}.
```

For a self term this reduces to
``P_{0,ii}=(2\pi\varepsilon_0)^{-1}\ln(2h_i/r_i)`` under the usual
thin-conductor approximation.
"""
function ametani_space2021(functor, pair)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    return log(geometry.D_ij / geometry.d_ij) / (2π * state.epsilon[1])
end

:Ametani2021
