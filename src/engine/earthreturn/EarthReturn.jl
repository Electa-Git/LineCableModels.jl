"""
    LineCableModels.Engine.EarthReturn

Compute the earth-return impedances and potential coefficients of cables. A formula gives
the impedance, the potential coefficients, or both, by the methods it defines. It treats the
pairs of conductors it declares, with `Val(1)` for the air and `Val(2)` for the earth.

# Formulas

## `Formula{:martins2026}`

Earth-return impedances and potential coefficients of conductors in two homogeneous
half-spaces, the air, medium 1 with ``y > 0``, and the earth, medium 2 with ``y < 0``, for
fields with the dependence ``e^{-\\Gamma \\ell}`` along the line, where ``\\ell`` is the axial
coordinate. Each cable enters as a solid conductor bounded by its exterior circumference,
which lies wholly in one medium. The field of each conductor is that of one axial line
source, and each conductor takes the circumferential mean of the field on its circumference, the
zeroth azimuthal harmonic. The formula contains the solutions of the underground case
[5437464](@cite) and of the mixed overhead and underground case [MartinsBritto2024](@cite).
Its form follows [MartinsBritto2026](@cite).

For each medium ``t``, the complex conductivity ``\\widetilde\\sigma_t`` \\[S/m\\], the
permeability ``\\mu_t`` \\[H/m\\], the transverse propagation parameter ``\\kappa_t`` \\[1/m\\]
and the vertical decay parameter ``a_t`` \\[1/m\\] are

```math
\\widetilde\\sigma_t = \\frac{1}{\\rho_t} + s \\varepsilon_0 \\varepsilon_{r,t}, \\qquad
\\mu_t = \\mu_0 \\mu_{r,t}, \\qquad
\\kappa_t^2 = s \\mu_t \\widetilde\\sigma_t - \\Gamma^2, \\qquad
a_t = \\sqrt{\\lambda^2 + \\kappa_t^2},
```

where ``s = \\mathrm{j} 2 \\pi f`` is the complex frequency \\[1/s\\], ``\\rho_t`` is the
resistivity \\[Ω·m\\], ``\\varepsilon_{r,t}`` and ``\\mu_{r,t}`` are the relative permittivity
and permeability \\[dimensionless\\], ``\\Gamma`` is the longitudinal wavenumber \\[1/m\\], and
``\\lambda`` is the transverse Fourier variable \\[1/m\\]. The constants ``\\varepsilon_0`` and
``\\mu_0`` are [`ε₀`](@ref) and [`μ₀`](@ref). The square
roots take the decaying branch. The spectral denominators are

```math
\\Delta_\\mu = \\mu_2 a_1 + \\mu_1 a_2, \\qquad
\\Delta_\\sigma = \\widetilde\\sigma_2 a_1 + \\widetilde\\sigma_1 a_2.
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
    \\frac{\\widetilde\\sigma_s (\\mu_s a_s + \\mu_{\\bar s} a_{\\bar s})}{\\Delta_\\mu \\Delta_\\sigma}
    e^{-a_s (h_i + h_j)} \\cos(\\lambda x_{ij}) \\, \\mathrm{d}\\lambda,
\\end{aligned}
```

where ``I_\\nu`` and ``K_\\nu`` are the modified Bessel functions of order ``\\nu``. The axial
vector-potential kernel ``\\widetilde L_{ij}`` \\[H/m\\] and the scalar-potential kernel
``\\widetilde\\Phi_{ij}`` \\[m/F\\] are

```math
\\widetilde L_{ij} = \\frac{\\mu_s}{2 \\pi}
    \\left[\\Lambda_{ij} + 2 I_0(\\kappa_s r_i) J_{ij}\\right], \\qquad
\\widetilde\\Phi_{ij} = \\frac{s}{2 \\pi \\widetilde\\sigma_s}
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
\\widetilde P_{ij} &= \\delta_{ts} \\frac{s}{2 \\pi \\widetilde\\sigma_2} \\Lambda_{ij}
    + \\frac{s}{\\pi} I_0(\\kappa_2 r_i) \\int_0^\\infty \\frac{a_1}{a_2 \\Delta_\\sigma}
    e^{-a_2 h_i - a_s h_j} \\cos(\\lambda x_{ij}) \\, \\mathrm{d}\\lambda,
    & t &= 2, \\\\
\\widetilde P_{ij} &= \\delta_{ts} \\frac{s}{2 \\pi \\widetilde\\sigma_1} \\Lambda_{ij}
    + \\frac{s}{\\pi} \\int_0^\\infty \\frac{a_2}{a_1 \\Delta_\\sigma}
    \\left[I_0(\\kappa_1 r_i) e^{-a_1 h_i} - J_0(\\lambda r_i)\\right]
    e^{-a_s h_j} \\cos(\\lambda x_{ij}) \\, \\mathrm{d}\\lambda,
    & t &= 1,
\\end{aligned}
```

where ``\\delta`` is the Kronecker delta and ``J_0`` is the Bessel function of order zero,
which averages the interface endpoint over the horizontal extent of the circumference. For
``s = 1``, the integral of the line ``t = 1`` is ``-\\Theta_{ij} / \\widetilde\\sigma_1``, with
``\\Theta_{ij}`` the averaged interface correction of [MartinsBritto2026](@cite). The kernels
evaluate it as one integral for a nominal ``|\\operatorname{Re} \\kappa_1 r_i| < 300`` and as its
two endpoint terms above.

The enclosed currents follow from the complete field on the circumferences, with the reflected and mutual
fields, through the dimensionless current matrix ``\\mathbf T``, the key `T` of the formula's
block, and the diagonal matrices ``\\mathbf W`` \\[dimensionless\\] and ``\\mathbf F``
\\[m/Ω\\]. Elimination of the source amplitudes gives the earth-return potential coefficients
``\\mathbf P_g`` \\[m/F\\] and impedances ``\\mathbf Z_g`` \\[Ω/m\\]:

```math
\\mathbf T = \\mathbf W^{-1} - \\mathbf F \\widetilde{\\mathbf K}, \\qquad
W_{ii} = I_0(\\kappa_t r_i), \\qquad
F_{ii} = \\frac{2 \\pi \\widetilde\\sigma_t r_i I_1(\\kappa_t r_i)}{\\kappa_t I_0(\\kappa_t r_i)},
```

```math
\\mathbf P_g = \\widetilde{\\mathbf P} \\mathbf T^{-1}, \\qquad
\\mathbf Z_g = \\left(\\widetilde{\\mathbf K}
    + \\frac{\\Gamma^2}{s} \\widetilde{\\mathbf P}\\right) \\mathbf T^{-1}.
```

The internal impedances of the conductors are not part of ``\\mathbf Z_g``: the impedance
inside each cable, `Zin`, gives them.

The integrals run along the complex path that each [`SpectralIntegral`](@ref) chooses, and
[`integrate`](@ref) evaluates them. The kernel ``\\widetilde\\Phi_{ij}`` enters only for
``\\Gamma \\neq 0``. With ``\\xi = \\kappa_t r_i``, the ratio ``I_1(\\xi) / (\\xi I_0(\\xi))`` in
``F_{ii}`` is the series ``1/2 - \\xi^2/16 + \\xi^4/96 - 11 \\xi^6/6144`` for a nominal
``|\\xi| < 10^{-3}``. The direct and image terms are evaluated together for a nominal
``|\\kappa_s D'_{ij}| < 0.5``, and their difference is ``\\ln(D'_{ij} / D_{ij})`` for
``\\kappa_s = 0``, with ``D_{ii} = r_i``. The implementation scales the column of each
source ``j`` by ``e^{|\\operatorname{Re} \\kappa_s r_j|}`` with exponentially scaled Bessel
functions, and the scaling cancels in ``\\mathbf P_g`` and ``\\mathbf Z_g``.

### Callables

The kernels are the integrands of the integrals of a pair in ``\\lambda``, with the states
`medium` and `pair` that [`to_buffer!`](@ref) builds:

- `formula(Val(:J), λ, medium, pair)`: the magnetic coupling.
- `formula(Val(:M), λ, medium, pair)`: the scalar-potential correction.
- `formula(Val(:P̃), λ, medium, pair)`: the voltage correction.
- `formula(Val(:Θ), λ, medium, pair)`: the averaged interface correction of a target in the
  air.
- `formula(Val(:Θ), Val(:interface), λ, medium, pair)`: the interface endpoint of that
  correction.

The entries read the block of the formula, with `i` the target and `j` the source:

- `formula(Val(:Z), Val(:self), i, buffers, input, media, s, k)` and
  `formula(Val(:Z), Val(:mutual), i, j, buffers, input, media, s, k)`: an entry of
  ``\\mathbf Z_g`` \\[Ω/m\\].
- `formula(Val(:P), Val(:self), i, buffers, input, media, s, k)` and
  `formula(Val(:P), Val(:mutual), i, j, buffers, input, media, s, k)`: an entry of
  ``\\mathbf P_g`` \\[m/F\\].

## `Formula{:ideal}`

Potential coefficients of conductors over an earth of infinite conductivity, from the images of
the sources, for conductors in the air:

```math
P_{g,ii} = \\frac{\\ln(2 h_i / r_i)}{2 \\pi \\varepsilon_1}, \\qquad
P_{g,ij} = \\frac{\\ln(D'_{ij} / D_{ij})}{2 \\pi \\varepsilon_1},
```

where ``\\varepsilon_1 = \\varepsilon_0 \\varepsilon_{r,1}`` is the permittivity of the air
\\[F/m\\]. These are ``\\Lambda_{ii}`` and ``\\Lambda_{ij}`` at ``\\kappa = 0``. The earth shields
a buried conductor, so each coefficient of a buried or mixed pair is zero. The formula gives
the potential coefficients alone.

### Callables

`formula(Val(:P), Val(:self), i, buffers, input, media, s, k)` and
`formula(Val(:P), Val(:mutual), i, j, buffers, input, media, s, k)` return an entry of
``\\mathbf P_g`` \\[m/F\\], with `i` the target and `j` the source.

## `Formula{:gary1976}`

Impedance of conductors in the air over a nonmagnetic earth, from an image plane at the complex
depth ``p`` \\[m\\] [Gary1976, Deri1981](@cite):

```math
\\begin{aligned}
p &= \\frac{1}{\\sqrt{s \\mu_0 \\sigma_2}}, \\\\
Z_{g,ii} &= \\frac{s \\mu_0}{2 \\pi} \\left[\\ln\\frac{2 h_i}{r_i} + \\ln\\frac{h_i + p}{h_i}\\right], \\\\
Z_{g,ij} &= \\frac{s \\mu_0}{2 \\pi} \\left[\\ln\\frac{D'_{ij}}{D_{ij}}
    + \\ln\\frac{\\sqrt{y_{ij}^2 + (H_{ij} + 2 p)^2}}{D'_{ij}}\\right],
\\end{aligned}
```

where ``\\sigma_2 = 1 / \\rho_2`` \\[S/m\\] is the conductivity of the earth. The lengths \\[m\\] are the
height ``h_i``, the exterior radius ``r_i``, the horizontal distance ``y_{ij}`` and
``H_{ij} = h_i + h_j``. The first logarithm of each entry is ``\\Lambda_{ii}`` or ``\\Lambda_{ij}`` at ``\\kappa = 0``. The formula
leaves out the displacement current of the earth, and it gives the impedance alone. A computation refuses a
magnetic earth layer.

### Callables

`formula(Val(:Z), Val(:self), i, buffers, input, media, s, k)` and
`formula(Val(:Z), Val(:mutual), i, j, buffers, input, media, s, k)` return an entry of
``\\mathbf Z_g`` \\[Ω/m\\], with `i` the target and `j` the source.

## `Formula{:wedepohl1973}`

Impedance of buried conductors by the low-order expansion of Wedepohl and Wilcox
[Wedepohl1973](@cite). The earth has the permeability ``\\mu_2 = \\mu_0 \\mu_{r,2}``, with
``m = \\sqrt{s \\mu_2 / \\rho_2}`` \\[1/m\\], and ``e_c`` is the exponential of Euler's constant:

```math
\\begin{aligned}
Z_{g,ii} &= \\frac{s \\mu_2}{2 \\pi} \\left[-\\ln\\frac{e_c m r_i}{2} + \\frac{1}{2} - \\frac{4 m h_i}{3}\\right], \\\\
Z_{g,ij} &= \\frac{s \\mu_2}{2 \\pi} \\left[-\\ln\\frac{e_c m D_{ij}}{2} + \\frac{1}{2} - \\frac{2 m H_{ij}}{3}\\right],
\\end{aligned}
```

where ``h_i`` is the depth \\[m\\]. The formula leaves out the displacement current of the earth, and it
gives the impedance alone. The source bounds the entries by ``|m r_i| < 0.25`` and
``|m D_{ij}| < 0.25``. A computation warns when the highest frequency of its scan leaves that bound, with
``m`` from the static resistivity and permeability of the deepest layer.

### Callables

`formula(Val(:Z), Val(:self), i, buffers, input, media, s, k)` and
`formula(Val(:Z), Val(:mutual), i, j, buffers, input, media, s, k)` return an entry of
``\\mathbf Z_g`` \\[Ω/m\\], with `i` the target and `j` the source.

## `Formula{:saad1996}`

Impedance of buried conductors in a nonmagnetic earth by the closed form of Saad, Gaba and Giroux
[Saad1996](@cite), with ``m = \\sqrt{s \\mu_0 / \\rho_2}`` \\[1/m\\]:

```math
\\begin{aligned}
Z_{g,ii} &= \\frac{s \\mu_0}{2 \\pi} \\left[K_0(m r_i) + \\frac{2 e^{-2 m h_i}}{4 + m^2 r_i^2}\\right], \\\\
Z_{g,ij} &= \\frac{s \\mu_0}{2 \\pi} \\left[K_0(m D_{ij}) + \\frac{2 e^{-m H_{ij}}}{4 + m^2 y_{ij}^2}\\right],
\\end{aligned}
```

where ``K_0`` is the modified Bessel function of the second kind and order zero. The source proves
the form for ``y_{ij} / H_{ij} < 1``. The formula leaves out the displacement current of the earth, and it
gives the impedance alone. A computation refuses a magnetic earth layer.

### Callables

`formula(Val(:Z), Val(:self), i, buffers, input, media, s, k)` and
`formula(Val(:Z), Val(:mutual), i, j, buffers, input, media, s, k)` return an entry of
``\\mathbf Z_g`` \\[Ω/m\\], with `i` the target and `j` the source.

## `Formula{:lucca1994}`

Mutual impedance between a conductor in the air and a conductor in a nonmagnetic earth, by the
approximation of Lucca [Lucca1994](@cite). Equation (28) of [Ametani2009](@cite) prints the expression.
With the complex depth ``p = 1 / \\sqrt{s \\mu_0 \\sigma_2}`` \\[m\\]:

```math
\\begin{aligned}
V &= h_a + h_g, \\qquad H_c = V + 2 p, \\qquad
S = \\sqrt{H_c^2 + y_{ij}^2}, \\qquad D = \\sqrt{V^2 + y_{ij}^2}, \\\\
Z_{g,ij} &= \\frac{s \\mu_0}{2 \\pi} \\left[\\ln\\frac{S}{D}
    - \\frac{2}{3} \\left(\\frac{p}{S^2}\\right)^3 H_c \\left(H_c^2 - 3 y_{ij}^2\\right)\\right],
\\end{aligned}
```

where ``h_a`` is the height of the conductor in the air and ``h_g`` the depth of the conductor in
the earth \\[m\\]. The horizontal distance is ``y_{ij}`` \\[m\\], and ``\\sigma_2 = 1 / \\rho_2`` \\[S/m\\].
The formula leaves out the displacement current of the earth. It gives the entries of a mixed pair
alone, so other formulas fill the air and earth places. A computation refuses a magnetic earth layer.

### Callables

`formula(Val(:Z), Val(:mutual), i, j, buffers, input, media, s, k)` returns an entry of
``\\mathbf Z_g`` \\[Ω/m\\], with `i` the target and `j` the source, in either order of a mixed pair.

## `Formula{:ametani2009}`

Mutual impedance between a conductor in the air and a conductor in a nonmagnetic earth, by
equation (27) of Ametani, Yoneda, Baba and Nagaoka [Ametani2009](@cite). With
``m = \\sqrt{s \\mu_0 / \\rho_2}`` \\[1/m\\] and the complex depth ``p = 1 / m`` \\[m\\]:

```math
\\begin{aligned}
S &= \\sqrt{(h_a + h_g + 2 p)^2 + y_{ij}^2}, \\qquad
D = \\sqrt{(h_a + h_g)^2 + y_{ij}^2}, \\\\
Z_{g,ij} &= \\frac{s \\mu_0}{2 \\pi} e^{-h_g / p} \\ln\\frac{S}{D},
\\end{aligned}
```

where ``h_a``, ``h_g`` and ``y_{ij}`` are the lengths of `Formula{:lucca1994}`. The exponential takes
``h_g`` as a positive depth. The formula leaves out the displacement current of the earth. It gives
the entries of a mixed pair alone, so other formulas fill the air and earth places. A computation
refuses a magnetic earth layer.

### Callables

`formula(Val(:Z), Val(:mutual), i, j, buffers, input, media, s, k)` returns an entry of
``\\mathbf Z_g`` \\[Ω/m\\], with `i` the target and `j` the source, in either order of a mixed pair.

## `Formula{:carson1926}`

Impedance of conductors in the air over a nonmagnetic earth, by the correction integral of Carson
[Carson1926](@cite), with ``\\gamma^2 = s \\mu_0 \\sigma_2`` \\[1/m²\\]:

```math
\\begin{aligned}
J_{ij} &= \\int_0^\\infty \\frac{e^{-H_{ij} \\lambda} \\cos(y_{ij} \\lambda)}
    {\\lambda + \\sqrt{\\lambda^2 + \\gamma^2}} \\, \\mathrm{d}\\lambda, \\\\
Z_{g,ii} &= \\frac{s \\mu_0}{2 \\pi} \\left[\\ln\\frac{2 h_i}{r_i} + 2 J_{ii}\\right], \\\\
Z_{g,ij} &= \\frac{s \\mu_0}{2 \\pi} \\left[\\ln\\frac{D'_{ij}}{D_{ij}} + 2 J_{ij}\\right],
\\end{aligned}
```

where ``\\sigma_2 = 1 / \\rho_2`` \\[S/m\\], ``J_{ij}`` is in \\[dimensionless\\], and the lengths are
those of `Formula{:gary1976}`. The first logarithm of each entry is ``\\Lambda_{ii}`` or
``\\Lambda_{ij}`` at ``\\kappa = 0``. The integral of the source is ``\\mathrm{j} J_{ij} = P + \\mathrm{j} Q``,
at ``p = H_{ij} \\sqrt{\\omega \\mu_0 \\sigma_2}`` and ``q = y_{ij} \\sqrt{\\omega \\mu_0 \\sigma_2}``. Its series
evaluate that integral, and this formula evaluates it on the path of [`SpectralIntegral`](@ref), under the
options `integral` of [`integrate`](@ref). The formula leaves out the displacement currents of the air
and the earth, and it gives the impedance alone. A computation refuses a magnetic earth layer.

### Callables

`formula(Val(:J), λ, γ², H, y)` returns the integrand of ``J_{ij}`` \\[m\\] at the spectral coordinate
`λ` \\[1/m\\], with `γ²` \\[1/m²\\], `H` and `y` \\[m\\]. `formula(Val(:Z), Val(:self), i, buffers, input, media, s, k)`
and `formula(Val(:Z), Val(:mutual), i, j, buffers, input, media, s, k)` return an entry of
``\\mathbf Z_g`` \\[Ω/m\\], with `i` the target and `j` the source.
"""
module EarthReturn

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults, validate
using ...Commons: FormulationOptions
using ...DataModel: support
using ..Engine: AbstractEarthReturn, CoaxialInput, LineParametersProblem, SpectralIntegral,
                integrate, nominal, ε₀, μ₀
import ..Engine: impedance!, potential!, to_buffer!, from_buffer, initialize_buffers, Λ
using LinearAlgebra: ldiv!, lu!
using SpecialFunctions: besselix, besselkx, besseljx, besselk

export Formula, FORMULAS

"Identifiers of the earth-return formulas."
const FORMULAS = (
    :martins2026, :ideal, :gary1976, :wedepohl1973, :saad1996, :lucca1994, :ametani2009,
    :carson1926)

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

# A formula that does not read Γ takes none: a result records its formulation.
function validate(Γ::Union{Number, AbstractVector}, formula::Formula{ID}) where {ID}
    all(iszero, Γ) || throw(DomainError(Γ,
        "the earth-return formula :$ID does not read Γ; a nonzero Γ needs a formula that " *
        "models the longitudinal dependence in every earth place"))
    return Γ
end

# A formula treats the pairs it declares, each a target and a source, `Val{1}` in the air and
# `Val{2}` in the earth.
function validate(formula::Formula{ID}, ::Val{t}, ::Val{s}) where {ID, t, s}
    throw(DomainError(formula,
        "the earth-return formula :$ID treats no pair of a receiver in the " *
        "$(t == 1 ? "air" : "earth") and a source in the $(s == 1 ? "air" : "earth"); " *
        "select a formula for that domain"))
end

# Formulas of one identifier are equal when their options are equal.
Base.:(==)(a::Formula{ID}, b::Formula{ID}) where {ID} = a.options == b.options

# A formula's block in the buffers, keyed by its identifier.
from_buffer(::Formula{ID}, buffers) where {ID} = getproperty(buffers, ID)

# The direct and image mean of the family, the paper's Λ_ii and Λ_ij. Below a nominal |κ D′|
# of 0.5 it evaluates the two terms together, with their logarithmic limit at κ = 0. Above, it
# combines the exponents of the exponentially scaled Bessel functions.
function Λ(::Formula, ::Val{:self}, κ, r, h)
    D′ = 2h
    if abs(nominal(κ * D′)) < 0.5
        iszero(κ) && return complex(log(D′ / r))
        return (besselk(0, κ * r) - besselk(0, κ * D′)) - i0m1(κ * r) * besselk(0, κ * D′)
    end
    ξ = κ * r
    return besselkx(0, ξ) * exp(-ξ) -
           besselix(0, ξ) * besselkx(0, κ * D′) * exp(abs(real(ξ)) - κ * D′)
end
function Λ(::Formula, ::Val{:mutual}, κ, r, D, D′)
    if abs(nominal(κ * D′)) < 0.5
        iszero(κ) && return complex(log(D′ / D))
        return besselix(0, κ * r) * exp(abs(real(κ * r))) *
               (besselk(0, κ * D) - besselk(0, κ * D′))
    end
    ξ = κ * r
    return besselix(0, ξ) * (besselkx(0, κ * D) * exp(abs(real(ξ)) - κ * D) -
            besselkx(0, κ * D′) * exp(abs(real(ξ)) - κ * D′))
end

include("martins2026.jl")
include("ideal.jl")
include("gary1976.jl")
include("wedepohl1973.jl")
include("saad1996.jl")
include("lucca1994.jl")
include("ametani2009.jl")
include("carson1926.jl")

end
