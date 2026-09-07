function routes(identifier::Val{:Wait1978})
    return (
        self=FormulaMethod(identifier,earth_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,earth_impedance,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:Wait1978})=(air=_lossless,earth=_full,permeability=vacuum_permeability)
propagation(::Val{:Wait1978})=Val(:zero)
function Formula(::Val{:Wait1978};approximation::Symbol=:bessel,kwargs...)
    approximation in (:bessel,:small_argument) || throw(ArgumentError(
        ":Wait1978 approximation must be :bessel or :small_argument"))
    identifier=Val(:Wait1978); defaults=routes(identifier)
    if approximation === :small_argument
        defaults=merge(defaults,(
            self=FormulaMethod(identifier,earth_impedance,Val(:small_argument)),))
    end
    overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) || throw(ArgumentError(
        "unknown routes for earth-impedance formula :Wait1978"))
    return Formula(identifier,merge(defaults,overrides),assumptions(identifier))
end

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Buried insulated thin cable with core radius ``c_0``, outer insulation radius ``c``, and burial depth ``h``. |
| Calculated quantities | Full-wave spectral external impedance ``Z_e(\\beta)`` and low-frequency earth-return reduction |
| Earth structure | Homogeneous earth below air. |
| Model and approximation | The parent expression retains longitudinal propagation and displacement current in both media. The low-frequency form uses the source's stated good-ground and electrically thin limits. |
| Main source | James R. Wait (1978) |
| Citation key(s) | `:Wait1978` |
| Evidence status | Original PDF equations (5), (6), and (10)–(18) checked |

**Implemented scope.** The low-frequency self reduction (12), with the
alternative small-radius expression (13) selected by
`approximation=:small_argument`, supplies an external coefficient
independent of the unknown modal propagation constant. The full-wave
dispersion problem (5)–(6) is not implemented by this selection.
No mutual expression is inferred from the single-wire approximation.

**Expression.** For ``g=jk=\\sqrt{j\\omega\\mu_0(\\sigma+j\\omega\\varepsilon)}``
and ``a=2gh``,

```math
Z_e=\\frac{j\\omega\\mu_0}{2\\pi}
\\left[K_0(gc)+K_2(a)-\\frac{2(1+a)e^{-a}}{a^2}\\right].
```

The recurrence ``K_2(a)=K_0(a)+2K_1(a)/a`` gives the source's
equation (12). Equation (13) replaces the direct term with
``-\\ln(0.89gc)``; the printed decimal remains part of that approximation.
The combined image term is evaluated without subtracting divergent
small-argument contributions. This numerical helper is shared with the
later Sunde-image approximations.

**Reference.** [Wait1978](@cite), equations (12)–(14), pp. 876–877.
"""
description(::Formula{:Wait1978})="Wait buried-wire low-frequency self impedance (1978)"
propagation_constant(::Val{:Wait1978},jω,permeability,permittivity)=(Γ=zero(jω),squared=zero(jω))
function (formula::Formula{:Wait1978})(rho,epsilon,mu,jω,Γ,segments=nothing)
    return _homogeneous_functor(Val(:Wait1978),formula,rho,epsilon,mu,jω,Γ,segments)
end
function earth_impedance(::Val{:Wait1978},::Val{:mutual},functor,pair)
    throw(ArgumentError("Wait1978 low-frequency reduction supplies only a buried self coefficient"))
end
function earth_impedance(identifier::Val{:Wait1978},::Val{:self},functor,pair)
    return earth_impedance(identifier,Val(:reduced),functor,pair,false)
end
function earth_impedance(identifier::Val{:Wait1978},::Val{:small_argument},functor,pair)
    return earth_impedance(identifier,Val(:reduced),functor,pair,true)
end
function earth_impedance(::Val{:Wait1978},::Val{:reduced},functor,pair,small)
    _require(pair,Val(:underground))
    pair.row==pair.column || throw(ArgumentError("Wait1978 requires self geometry"))
    state=functor.state; geometry=_geometry(pair); g=state.gamma[2]
    direct=small ? -log(oftype(geometry.y_ij,0.89)*g*geometry.y_ij) :
        special_besselk(0,g*geometry.y_ij)
    correction=_sunde_image_correction(g,geometry.H,zero(geometry.H))
    return state.jω*state.mu[2]/(2*(one(geometry.H)*π))*(direct+correction)
end

:Wait1978
