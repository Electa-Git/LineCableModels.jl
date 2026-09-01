function routes(identifier::Val{:Ametani2021})
    (
        self = formula_method(identifier, earth_potential_coefficient, Val(:self)),
        mutual = formula_method(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = formula_method(identifier, propagation_constant)
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
"""
$(TYPEDSIGNATURES)

**Identification.** Classical electrostatic image coefficient for overhead
conductors above an equipotential plane.

**Expression.**

```math
P_{0,ij}=\\frac{1}{2\\pi\\varepsilon_0}\\ln\\frac{D_{ij}}{d_{ij}},
\\qquad \\mathbf Y_0=j\\omega\\mathbf P_0^{-1}.
```

For a self term, ``P_{0,ii}=(2\\pi\\varepsilon_0)^{-1}\\ln(2h_i/r_i)``.

**Reference.** A. Ametani, H. Xue, T. Ohno, and H. Khalilnezhad,
*Electromagnetic Transients in Large HV Cable Networks: Modeling and
Calculations*, IET, 2021.
"""
function description(::Formula{:Ametani2021})
    "Ametani et al. classical overhead space potential coefficient (2021)"
end

function propagation_constant(
        ::Val{:Ametani2021}, jω, permeability, permittivity
)
    return (Γ = zero(jω), squared = zero(jω))
end

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
function earth_potential_coefficient(
        ::Val{:Ametani2021}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    return log(geometry.D_ij / geometry.d_ij) / (2π * state.epsilon[1])
end

:Ametani2021
