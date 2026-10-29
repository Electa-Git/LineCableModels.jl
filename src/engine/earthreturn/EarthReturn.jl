"""
    LineCableModels.Engine.EarthReturn

Compute the earth-return impedances and potential coefficients of cables. A formula gives
the impedance, the potential coefficients, or both, by the methods it defines.

# Formulas

## `Formula{:unified}`

Earth-return impedances and potential coefficients of conductors in two homogeneous
half-spaces, the air, medium 1 with ``y > 0``, and the earth, medium 2 with ``y < 0``, for
fields with the dependence ``e^{-\\Gamma \\ell}`` along the line, where ``\\ell`` is the axial
coordinate. Each cable enters as a solid conductor bounded by its exterior circumference,
which lies wholly in one medium. The field of each conductor is that of one axial line
source, and each conductor takes the circumferential mean of the field on its circumference, the
zeroth azimuthal harmonic. The formula contains the solutions of the underground case
[5437464](@cite) and of the mixed overhead and underground case [MartinsBritto2024](@cite).

For each medium ``t``, the complex conductivity ``\\widehat\\sigma_t`` \\[S/m\\], the
permeability ``\\mu_t`` \\[H/m\\], the transverse propagation parameter ``\\kappa_t`` \\[1/m\\]
and the vertical decay parameter ``a_t`` \\[1/m\\] are

```math
\\widehat\\sigma_t = \\frac{1}{\\rho_t} + s \\varepsilon_0 \\varepsilon_{r,t}, \\qquad
\\mu_t = \\mu_0 \\mu_{r,t}, \\qquad
\\kappa_t^2 = s \\mu_t \\widehat\\sigma_t - \\Gamma^2, \\qquad
a_t = \\sqrt{\\lambda^2 + \\kappa_t^2},
```

where ``s = \\mathrm{j} 2 \\pi f`` is the complex frequency \\[1/s\\], ``\\rho_t`` is the
resistivity \\[Ω·m\\], ``\\varepsilon_{r,t}`` and ``\\mu_{r,t}`` are the relative permittivity
and permeability \\[dimensionless\\], ``\\Gamma`` is the longitudinal wavenumber \\[1/m\\], and
``\\lambda`` is the transverse Fourier variable \\[1/m\\]. The constants ``\\varepsilon_0`` and
``\\mu_0`` are [`vacuum_permittivity`](@ref) and [`vacuum_permeability`](@ref). The square
roots take the decaying branch. The spectral denominators are

```math
\\Delta_\\mu = \\mu_2 a_1 + \\mu_1 a_2, \\qquad
\\Delta_\\sigma = \\widehat\\sigma_2 a_1 + \\widehat\\sigma_1 a_2.
```

Conductor ``i`` has its center at ``(x_i, y_i)``, the distance ``h_i = |y_i|`` to the
interface, the exterior radius ``r_i`` \\[m\\] and the medium ``t = m(i)``. In a coefficient
``ij``, conductor ``i`` is the target and conductor ``j`` the source, in the medium
``s = m(j)``, and ``\\bar s`` is the other medium. The distance between the centers is
``D_{ij}``, the distance from conductor ``i`` to the image of conductor ``j`` is
``D'_{ij} = \\sqrt{x_{ij}^2 + (h_i + h_j)^2}``, with ``D'_{ii} = 2 h_i``, and
``x_{ij} = x_i - x_j`` \\[m\\].

For conductors in the same medium, ``t = s``, the means of the direct and image fields and
the interface corrections are

```math
\\begin{aligned}
\\Lambda_{ii} &= K_0(\\kappa_s r_i) - I_0(\\kappa_s r_i) K_0(2 \\kappa_s h_i), \\qquad
\\Lambda_{ij} = I_0(\\kappa_s r_i)
    \\left[K_0(\\kappa_s D_{ij}) - K_0(\\kappa_s D'_{ij})\\right], \\\\
J_{ij} &= \\int_0^\\infty \\frac{\\mu_{\\bar s}}{\\Delta_\\mu}
    e^{-a_s (h_i + h_j)} \\cos(\\lambda x_{ij}) \\, \\mathrm{d}\\lambda, \\qquad
M_{ij} = \\int_0^\\infty
    \\frac{\\widehat\\sigma_s (\\mu_s a_s + \\mu_{\\bar s} a_{\\bar s})}{\\Delta_\\mu \\Delta_\\sigma}
    e^{-a_s (h_i + h_j)} \\cos(\\lambda x_{ij}) \\, \\mathrm{d}\\lambda,
\\end{aligned}
```

where ``I_\\nu`` and ``K_\\nu`` are the modified Bessel functions of order ``\\nu``. The axial
vector-potential kernel ``\\widetilde L_{ij}`` \\[H/m\\] and the scalar-potential kernel
``\\widetilde\\Phi_{ij}`` \\[m/F\\] are

```math
\\widetilde L_{ij} = \\frac{\\mu_s}{2 \\pi}
    \\left[\\Lambda_{ij} + 2 I_0(\\kappa_s r_i) J_{ij}\\right], \\qquad
\\widetilde\\Phi_{ij} = \\frac{s}{2 \\pi \\widehat\\sigma_s}
    \\left[\\Lambda_{ij} + 2 I_0(\\kappa_s r_i) M_{ij}\\right]
```

for ``t = s``, and, from the transmitted fields,

```math
\\begin{aligned}
\\widetilde L_{ij} &= \\frac{\\mu_t \\mu_s}{\\pi} I_0(\\kappa_t r_i) \\int_0^\\infty
    \\frac{e^{-a_t h_i - a_s h_j}}{\\Delta_\\mu} \\cos(\\lambda x_{ij}) \\, \\mathrm{d}\\lambda, \\\\
\\widetilde\\Phi_{ij} &= \\frac{s}{\\pi} I_0(\\kappa_t r_i) \\int_0^\\infty
    \\frac{\\mu_1 a_1 + \\mu_2 a_2}{\\Delta_\\mu \\Delta_\\sigma}
    e^{-a_t h_i - a_s h_j} \\cos(\\lambda x_{ij}) \\, \\mathrm{d}\\lambda
\\end{aligned}
```

for ``t \\neq s``. The kernel
``\\widetilde K_{ij} = s \\widetilde L_{ij} - (\\Gamma^2 / s) \\widetilde\\Phi_{ij}`` \\[Ω/m\\] is
the negative mean axial electric field at conductor ``i`` per source amplitude at conductor
``j``.

The voltage of a conductor is the circumferential mean of the integral of the vertical
electric field from a reference level. A conductor in the earth is referred to infinite
depth, and a conductor in the air to the interface below it. The source-normalized voltage
coefficients ``\\widetilde P_{ij}`` \\[m/F\\] are

```math
\\begin{aligned}
\\widetilde P_{ij} &= \\delta_{ts} \\frac{s}{2 \\pi \\widehat\\sigma_2} \\Lambda_{ij}
    + \\frac{s}{\\pi} I_0(\\kappa_2 r_i) \\int_0^\\infty \\frac{a_1}{a_2 \\Delta_\\sigma}
    e^{-a_2 h_i - a_s h_j} \\cos(\\lambda x_{ij}) \\, \\mathrm{d}\\lambda,
    & t &= 2, \\\\
\\widetilde P_{ij} &= \\delta_{ts} \\frac{s}{2 \\pi \\widehat\\sigma_1} \\Lambda_{ij}
    + \\frac{s}{\\pi} \\int_0^\\infty \\frac{a_2}{a_1 \\Delta_\\sigma}
    \\left[I_0(\\kappa_1 r_i) e^{-a_1 h_i} - J_0(\\lambda r_i)\\right]
    e^{-a_s h_j} \\cos(\\lambda x_{ij}) \\, \\mathrm{d}\\lambda,
    & t &= 1,
\\end{aligned}
```

where ``\\delta`` is the Kronecker delta and ``J_0`` is the Bessel function of order zero,
which averages the interface endpoint over the horizontal extent of the circumference.

The enclosed currents follow from the complete field on the circumferences, with the reflected and mutual
fields, through the dimensionless current matrix ``\\mathbf T``, held as `L` in the state of
the formula, and the diagonal matrices ``\\mathbf W`` \\[dimensionless\\] and ``\\mathbf F``
\\[m/Ω\\]. Elimination of the source amplitudes gives the earth-return potential coefficients
``\\mathbf P_g`` \\[m/F\\] and impedances ``\\mathbf Z_g`` \\[Ω/m\\]:

```math
\\mathbf T = \\mathbf W^{-1} - \\mathbf F \\widetilde{\\mathbf K}, \\qquad
W_{ii} = I_0(\\kappa_t r_i), \\qquad
F_{ii} = \\frac{2 \\pi \\widehat\\sigma_t r_i I_1(\\kappa_t r_i)}{\\kappa_t I_0(\\kappa_t r_i)},
```

```math
\\mathbf P_g = \\widetilde{\\mathbf P} \\mathbf T^{-1}, \\qquad
\\mathbf Z_g = \\left(\\widetilde{\\mathbf K}
    + \\frac{\\Gamma^2}{s} \\widetilde{\\mathbf P}\\right) \\mathbf T^{-1}.
```

The internal impedances of the conductors are not part of ``\\mathbf Z_g``: the impedance
inside each cable, `Zin`, gives them.

The integrals run along the complex ray that [`ray!`](@ref) selects, and
[`integrate`](@ref) evaluates them. The kernel ``\\widetilde\\Phi_{ij}`` enters only for
``\\Gamma \\neq 0``. With ``\\xi = \\kappa_t r_i``, the ratio ``I_1(\\xi) / (\\xi I_0(\\xi))`` in
``F_{ii}`` is the series ``1/2 - \\xi^2/16 + \\xi^4/96 - 11 \\xi^6/6144`` for a nominal
``|\\xi| < 10^{-3}``. The direct and image terms are evaluated together for a nominal
``|\\kappa_s D'_{ij}| < 0.5``, and their difference is ``\\ln(D'_{ij} / D_{ij})`` for
``\\kappa_s = 0``, with ``D_{ii} = r_i``. The implementation scales the column of each
source ``j`` by ``e^{|\\operatorname{Re} \\kappa_s r_j|}`` with exponentially scaled Bessel
functions, and the scaling cancels in ``\\mathbf P_g`` and ``\\mathbf Z_g``.
"""
module EarthReturn

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults, validate
using ...Commons: FormulationOptions
using ...DataModel: support
using ..Engine: AbstractEarthReturn, CoaxialInput, LineParametersProblem, SpectralIntegral,
                integrate, nominal, ray!, vacuum_permittivity, vacuum_permeability
import ..Engine: impedance!, potential!, prepare!, initialize_buffers
using LinearAlgebra: ldiv!, lu!
using SpecialFunctions: besselix, besselkx, besseljx, besselk

export Formula, FORMULAS

"Identifiers of the earth-return formulas."
const FORMULAS = (:unified,)

"""
$(TYPEDEF)

Represent an earth-return impedance and potential-coefficient formula.
Its options form a record of type `FormulationOptions{Formula{ID}}`, with `ID` the formula
identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, O <: FormulationOptions} <: AbstractEarthReturn
    "Options of the formula."
    options::O
    function Formula{ID}(options::O) where {ID, O <: FormulationOptions{Formula{ID}}}
        return new{ID, O}(options)
    end
end

"""
$(TYPEDSIGNATURES)

Build the formula named by `ID` with its declared defaults and supplied options.

# Arguments

- `Val{ID}`: formula identifier from [`FORMULAS`](@ref).

# Keywords

- `kw...`: options declared by [`defaults`](@ref) for `Formula{ID}`.

# Returns

- A [`Formula`](@ref) with a typed option record.

# Errors

- Throws `ArgumentError` for an unknown identifier or option.
- Throws `ArgumentError` for invalid integration controls.
"""
Formula(::Val{ID}; kw...) where {ID} = Formula{ID}(FormulationOptions(Formula{ID}; kw...))

# An identifier without a formula has no defaults.
function defaults(::Type{Formula{ID}}) where {ID}
    throw(ArgumentError(
        "unknown earth-return formula :$ID; the formulas are $FORMULAS"))
end

# A formula gives the quantities it has a method for, and a slot refuses a formula that does
# not give its quantity.
function validate(formula::Formula{ID}, verb) where {ID}
    throw(DomainError(formula, "the earth-return formula :$ID gives no $(nameof(verb))"))
end

# Formulas of one identifier are equal when their options are equal.
Base.:(==)(a::Formula{ID}, b::Formula{ID}) where {ID} = a.options == b.options

include("unified.jl")

end
