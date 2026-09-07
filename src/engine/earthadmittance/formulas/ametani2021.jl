function routes(identifier::Val{:Ametani2021})
    (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
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

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Thin overhead conductors or buried coaxial units with separately assembled radial insulation. |
| Calculated quantities | Constant overhead space potential; zero exterior potential for buried classical TL units |
| Earth structure | Equipotential reference plane; no finite-conductivity potential correction. |
| Model and approximation | Classical TL reference, not a lossy-earth potential solution. |
| Main source | A. Ametani, H. Xue, T. Ohno, and H. Khalilnezhad (2021) |
| Citation key(s) | `:Ametani2021` |
| Evidence status | IET book page image checked, Table 2.1 and equations (2.76)–(2.78), p. 24. |

**Description.** Classical TL potential reference: constant overhead space potential and zero
buried exterior correction, with radial insulation assembled separately.

**Expression.**

```math
P_{0,ij}=\\frac{1}{2\\pi\\varepsilon_0}\\ln\\frac{D_{ij}}{d_{ij}},
\\qquad \\mathbf Y_0=j\\omega\\mathbf P_0^{-1}.
```

For a self term, ``P_{0,ii}=(2\\pi\\varepsilon_0)^{-1}\\ln(2h_i/r_i)``.

**Reference.** [Ametani2021](@cite). A. Ametani, H. Xue, T. Ohno, and H. Khalilnezhad,
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
    _placement(pair) === Val(:overhead) || return zero(functor.state.jω)
    state = functor.state
    geometry = _geometry(pair)
    image_distance = pair.row == pair.column ? geometry.H : geometry.D_ij
    return log(image_distance / geometry.d_ij) / (2 * (one(geometry.H)*π) * state.epsilon[1])
end

:Ametani2021
