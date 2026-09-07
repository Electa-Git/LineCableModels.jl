function routes(identifier::Val{:Xue2021})
    (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        overhead = FormulaMethod(identifier, earth_potential_coefficient, Val(:overhead)),
        underground = FormulaMethod(identifier, earth_potential_coefficient, Val(:underground)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Xue2021})
    (
        air = _full,
        earth = _full,
        permeability = vacuum_permeability,
        evaluation = :ehem
    )
end

propagation(::Val{:Xue2021}) = Val(:zero)
propagation(formula::Formula{:Xue2021}) =
    formula.assumptions.evaluation===:exact ? Val(:explicit) : Val(:zero)
media(formula::Formula{:Xue2021}) =
    formula.assumptions.evaluation===:exact ? Val(:stratified) : Val(:homogeneous)

function Formula(::Val{:Xue2021};evaluation::Symbol=:ehem,kwargs...)
    evaluation in (:ehem,:exact) ||
        throw(ArgumentError(":Xue2021 potential evaluation must be :ehem or :exact"))
    identifier=Val(:Xue2021); defaults=routes(identifier); overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) ||
        throw(ArgumentError("unknown routes for earth-admittance formula :Xue2021"))
    selected=merge(defaults,overrides)
    values=merge(assumptions(identifier),(;evaluation,
        permeability=evaluation===:exact ? _material : vacuum_permeability))
    return Formula{:Xue2021,typeof(selected),typeof(values)}(selected,values)
end
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite parallel round overhead conductors. |
| Calculated quantities | Approximate N-layer overhead potential coefficients and assembled earth-return admittance via equivalent propagation |
| Earth structure | Arbitrary horizontal ``N``-layer earth. |
| Model and approximation | The layered electric interface response is replaced by one equivalent constant. The potential matrix is nevertheless assembled before inversion; no entrywise reciprocal is permitted. |
| Main source | H. Xue, J. Mahseredjian, A. Ametani, J. Morales, and I. Kocar (2021) |
| Citation key(s) | `:Xue2021` |
| Evidence status | Accepted-publication page image verified |

**Overhead EHEM expression.** The selected equivalent material reconstructs
the bulk soil constant. Recover the source's transverse constant before
evaluating equation (8):

```math
q^2=\\gamma_{soil,eq}^2-\\gamma_0^2,\\qquad
P_{e,ij}=\\frac{1}{2\\pi\\varepsilon_0}
\\left[\\ln(D_{ij}/d_{ij})+
2\\int_0^\\infty\\frac{e^{-H\\lambda}\\cos(x\\lambda)}
{(q^2/\\gamma_0^2)\\lambda+\\sqrt{\\lambda^2+q^2}}\\,d\\lambda\\right].
```

Compose this route with the Xue2021 equivalent-earth rule. The denominator
uses q², not the reconstructed bulk constant q²+γ₀². It is therefore not
identified with the Wise1948 material substitution.

**Exact four-layer selection.** With `evaluation=:exact`, the
following distinct source record applies.

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite parallel round overhead conductors; ``D_1,D_2`` from (2). |
| Calculated quantities | Four-layer overhead earth-return potential coefficients and assembled admittance matrix |
| Earth structure | Four horizontal earth layers below air. |
| Model and approximation | No fitted approximation is applied to (5); it is the source's exact four-layer potential entry within the selected longitudinal/TL model. |
| Main source | H. Xue, J. Mahseredjian, A. Ametani, J. Morales, and I. Kocar (2021) |
| Citation key(s) | `:Xue2021` |
| Evidence status | Accepted-publication page image verified |

The coefficient in source (5) uses the exact ``G`` of (6), not
the equivalent-earth value in (8). Magnetic and driven electric
interface equations are eliminated with decaying layer transmissions.
The electric forcing is retained and the two corrections are combined
before subtraction. Appendix B of the final publication is numbered
(25)–(44). Its complete expressions are independently compared with
the reduced numerical evaluation.

The exact selection requires air and four earth layers, lossless air,
overhead conductors, and the source air-reference longitudinal
wavenumber. Layer permeabilities remain independent. It is composed
with the actual layered earth, not an equivalent homogeneous earth.
Only the completed potential matrix is inverted.

**Reference.** [Xue2021](@cite), equations (4)–(6), (10)–(44).

**Retained underground selection.** The earlier underground closed form
below remains available for compatibility. It is the book-attributed Xue
approximation, not the overhead EHEM paper's equation (8).



**Expression.**

```math
P_{e,ij}=\\frac{j\\omega}{2\\pi(\\sigma_1+j\\omega\\varepsilon_1)}
\\left[\\ln\\left(\\frac{1+\\gamma_1R_{ab}}{\\gamma_1R_{ab}}\\right)+
\\frac{2e^{-H\\gamma_1}}{4+\\gamma_1^2R_{ab}^2}\\right].
```

**Reference.** A. Ametani, H. Xue, T. Ohno, and H. Khalilnezhad,
*Electromagnetic Transients in Large HV Cable Networks: Modeling and
Calculations*, IET, 2021, Section 2.6.3, equation (2.73). The book presents
this closed form as Xue's approximation based on the Saad and Petrache
expressions.
"""
function description(::Formula{:Xue2021})
    "Xue exact four-layer/EHEM overhead potential and retained underground approximation (2021)"
end

function propagation_constant(::Val{:Xue2021}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Xue2021})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    formula.assumptions.evaluation===:exact &&
        throw(DimensionMismatch(":Xue2021 exact potential requires four-layer thicknesses"))
    return _homogeneous_functor(
        Val(:Xue2021), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

function (formula::Formula{:Xue2021})(rho,epsilon,mu,s,Γ,segments,thickness)
    formula.assumptions.evaluation===:exact ||
        return formula(rho,epsilon,mu,s,Γ,segments)
    reference=EarthImpedance.earth_impedance(
        Val(:Xue2021),Val(:validate),rho,epsilon,mu,s,Γ,thickness)
    functor=_stratified_functor(Val(:Xue2021),formula,rho,epsilon,mu,s,reference.Γ,segments,thickness)
    state=merge(functor.state,(gamma_squared=-functor.state.gamma_medium_squared[1],))
    return Functor{:Xue2021,typeof(formula.routes),typeof(state)}(formula.routes,state)
end

raw"""
Evaluate Xue's closed-form underground potential coefficient:

```math
P_{e,ij}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
\left[\ln\left(\frac{1+\gamma_1R_{ab}}{\gamma_1R_{ab}}\right)+
\frac{2e^{-(h_i+h_j)\gamma_1}}{4+\gamma_1^2R_{ab}^2}\right].
```
"""
function earth_potential_coefficient(
        ::Val{:Xue2021}, ::Val{:mutual}, functor, pair
)
    functor.state.formula.assumptions.evaluation===:exact &&
        return earth_potential_coefficient(Val(:Xue2021),Val(:exact),functor,pair)
    placement=_placement(pair)
    placement===Val(:overhead) && return functor.routes.overhead(functor,pair)
    placement===Val(:underground) && return functor.routes.underground(functor,pair)
    throw(ArgumentError("Xue2021 potential formulas do not support mixed pairs"))
end

# Stable elimination of the exact four-layer Hertz-potential boundary
# equations. X is the magnetic potential at each interface. The electric
# response includes its interface forcing, not only a TM reflection factor.
function earth_potential_coefficient(::Val{:Xue2021},::Val{:kernel},lambda,state)
    count=length(state.rho)
    g=collect(state.gamma_medium_squared)
    q=[spectral_root(lambda^2+v-g[1],state.jω) for v in g]
    q[1]=complex(lambda)
    g./=maximum(abs,g)
    mu=state.mu./state.mu[1]
    input=q[end]/mu[end]; transmission=similar(q)
    for n in count-1:-1:2
        step=EarthImpedance._layered_input_step(q[n],inv(mu[n]),state.thickness[n],input)
        input=step.input; transmission[n]=step.transmission
    end
    F=2/(lambda+input)
    X=similar(q); X[1]=lambda*F
    for n in 2:count-1
        X[n]=transmission[n]*X[n-1]
    end
    electric=mu[end]*q[end]/g[end]; forcing=zero(state.jω)
    for n in count-1:-1:2
        step=EarthImpedance._layered_input_step(q[n],mu[n]/g[n],state.thickness[n],electric)
        forcing=(forcing+(inv(g[n])-inv(g[n+1]))*X[n])*step.transmission
        electric=step.input
    end
    # Substitute the top-interface forcing before subtraction. Computing
    # F plus the electric correction directly loses the small potential.
    G=(F*(electric+lambda/g[2])-forcing)/(lambda/g[1]+electric)
    return (;F,G)
end

function earth_potential_coefficient(::Val{:Xue2021},::Val{:exact},functor,pair)
    geometry=earth_potential_coefficient(Val(:Wise1948),Val(:geometry),pair)
    state=functor.state
    integral=_quadrature(state) do t
        lambda=t/geometry.H
        kernel=earth_potential_coefficient(Val(:Xue2021),Val(:kernel),lambda,state)
        kernel.G*exp(-t)*cos(geometry.x*lambda)/geometry.H
    end
    return _complex_result(state.jω,(log(geometry.D/geometry.d)+integral)/
        (2*(one(geometry.H)*π)*state.epsilon[1]))
end

function earth_potential_coefficient(
        ::Val{:Xue2021}, ::Val{:overhead}, functor, pair
)
    _require(pair,Val(:overhead))
    state=functor.state; geometry=_geometry(pair)
    gamma0_squared=state.gamma_medium_squared[1]
    equivalent_squared=state.gamma_medium_squared[2]-gamma0_squared
    iszero(gamma0_squared) && throw(DomainError(gamma0_squared,
        "Xue2021 overhead potential requires a nonzero air propagation constant"))
    integral=_quadrature(state) do λ
        denominator=equivalent_squared/gamma0_squared*λ+sqrt(λ^2+equivalent_squared)
        exp(-geometry.H*λ)*cos(geometry.y_ij*λ)/denominator
    end
    return (log(geometry.D_ij/geometry.d_ij)+2integral)/
        (2*(one(geometry.H)*π)*state.epsilon[1])
end

function earth_potential_coefficient(
        ::Val{:Xue2021}, ::Val{:underground}, functor, pair
)
    _require(pair, Val(:underground))
    pair.row == pair.column || _require_horizontal_separation(pair)
    state = functor.state
    geometry = _geometry(pair)
    gamma = state.gamma[2]
    argument = gamma * geometry.y_ij
    bracket = log((1 + argument) / argument) +
              2exp(-geometry.H * gamma) / (4 + argument^2)
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    return state.jω / (2π * kappa) * bracket
end

:Xue2021
