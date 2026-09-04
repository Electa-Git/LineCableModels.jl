function routes(identifier::Val{:Theodoulidis2015})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Theodoulidis2015})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Theodoulidis2015}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Closed form for Carson's overhead correction using
Struve and Bessel functions.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}
\\left[\\ln\\frac{D_{ij}}{d_{ij}}+S_{ij}\\right],
```

```math
S_{ij}=\\sum_{q\\in\\{1,2\\}}\\left[
\\frac{\\pi}{2u_q}(\\mathbf H_1(u_q)-Y_1(u_q))-\\frac1{u_q^2}\\right],
\\quad u_{1,2}=\\gamma_1(H\\mp jy_{ij}).
```

**Reference.** T. P. Theodoulidis, “On the Closed-Form Expression of Carson's
Integral,” *Periodica Polytechnica Electrical Engineering and Computer
Science*, 59(1), 26–29, 2015.
"""
function description(::Formula{:Theodoulidis2015})
    "Theodoulidis closed-form Carson correction (2015)"
end

function propagation_constant(
        ::Val{:Theodoulidis2015}, jω, permeability, permittivity
)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Theodoulidis2015})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Theodoulidis2015), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate Theodoulidis' closed form for Carson's correction:

```math
S_{ij}=\frac{\pi}{2u_1}[\mathbf H_1(u_1)-Y_1(u_1)]-\frac1{u_1^2}
+\frac{\pi}{2u_2}[\mathbf H_1(u_2)-Y_1(u_2)]-\frac1{u_2^2},
```

```math
u_1=\gamma_1(H-jy_{ij}),\qquad
u_2=\gamma_1(H+jy_{ij}),\qquad H=h_i+h_j.
```

The contribution is inserted exactly as
``Z_{e,ij}=(j\omega\mu_0/(2\pi))[\ln(D_{ij}/d_{ij})+S_{ij}]``.
The complex Struve function ``\mathbf H_1`` is evaluated by its convergent
power series because the sanctioned `SpecialFunctions` dependency does not
export Struve functions.
"""
function struve_h1(::Val{:Theodoulidis2015}, z::Complex{T}) where {T <: Real}
    term = 2z^2 / (3 * (one(T) * π))
    value = term
    tolerance = max(eps(T), T(1e-15))
    for k in 0:9999
        term *= -(z / 2)^2 / ((k + T(3) / 2) * (k + T(5) / 2))
        updated = value + term
        abs(term) <= tolerance * max(one(T), abs(updated)) && return updated
        value = updated
    end
    throw(ErrorException("complex Struve H1 series did not converge"))
end

@inline function closed_form_term(identifier::Val{:Theodoulidis2015}, u)
    return π / (2u) * (struve_h1(identifier, u) - bessely(1, u)) - inv(u^2)
end

function earth_impedance(
        identifier::Val{:Theodoulidis2015}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    u_1 = state.gamma[2] * (geometry.H - im * geometry.y_ij)
    u_2 = state.gamma[2] * (geometry.H + im * geometry.y_ij)
    correction = closed_form_term(identifier, u_1) +
                 closed_form_term(identifier, u_2)
    return state.jω * state.mu[1] / (2π) *
           (log(geometry.D_ij / geometry.d_ij) + correction)
end

:Theodoulidis2015
