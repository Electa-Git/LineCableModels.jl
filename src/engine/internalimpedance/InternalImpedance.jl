"""
    LineCableModels.Engine.InternalImpedance

Compute the inner, outer, and transfer surface impedances of a tubular conductor.

# Formulas

## `Formula{:schelkunoff1934}`

Surface impedances of a homogeneous tubular conductor of inner radius ``r_\\mathrm{in}`` and
outer radius ``r_\\mathrm{ex}``, from the cylindrical solution of the field inside the metal
[Schelkunoff1934](@cite), with the modified Bessel functions ``I_\\nu`` and ``K_\\nu``:

```math
\\begin{aligned}
z_\\mathrm{inner} &= \\frac{s \\mu}{2 \\pi w_\\mathrm{in}}
    \\frac{I_0(w_\\mathrm{in}) K_1(w_\\mathrm{ex}) + K_0(w_\\mathrm{in})
      I_1(w_\\mathrm{ex})}{D}, \\\\
z_\\mathrm{outer} &= \\frac{s \\mu}{2 \\pi w_\\mathrm{ex}}
    \\frac{I_0(w_\\mathrm{ex}) K_1(w_\\mathrm{in}) + K_0(w_\\mathrm{ex})
      I_1(w_\\mathrm{in})}{D}, \\\\
z_\\mathrm{transfer} &= \\frac{1}{2 \\pi r_\\mathrm{in} r_\\mathrm{ex} \\sigma D},
\\end{aligned}
```

where
``D = I_1(w_\\mathrm{ex}) K_1(w_\\mathrm{in}) - I_1(w_\\mathrm{in}) K_1(w_\\mathrm{ex})``,
``w = m r``, ``m = \\sqrt{s \\mu \\sigma}``, ``\\mu = \\mu_0 \\mu_r`` and
``\\sigma = 1 / \\rho``. Here ``s`` is the complex frequency \\[1/s\\], ``r`` is the radius
of the corresponding conductor surface \\[m\\], ``m`` is the radial wavenumber \\[1/m\\],
``\\sigma`` is the conductivity \\[S/m\\], ``\\rho`` is the resistivity \\[Ω·m\\],
``\\mu`` is the permeability \\[H/m\\] and ``\\mu_r`` is the relative permeability
\\[dimensionless\\].
The constant ``\\mu_0`` is [`μ₀`](@ref). The modified Bessel functions
``I_\\nu`` and ``K_\\nu`` have order ``\\nu``.

A solid conductor, with ``r_\\mathrm{in} = 0``, has no inner surface: its inner and transfer
impedances are zero, and its outer impedance is
``s \\mu I_0(w_\\mathrm{ex}) / (2 \\pi w_\\mathrm{ex} I_1(w_\\mathrm{ex}))``. The
implementation evaluates the Bessel functions in their exponentially scaled form.

At a harmonic order ``n``, the field varies around the conductor's axis as
``\\cos n\\varphi`` or ``\\sin n\\varphi``. The record then relates the axial electric field
``E`` to ``J = 2\\pi r H_\\varphi`` at the two surfaces:
``E(r_\\mathrm{in}) = -z_\\mathrm{inner} J_\\mathrm{in} + z_\\mathrm{transfer} J_\\mathrm{ex}``
and ``E(r_\\mathrm{ex}) = -z_\\mathrm{transfer} J_\\mathrm{in} + z_\\mathrm{outer} J_\\mathrm{ex}``.
At ``n = 0``, ``J_\\mathrm{in}`` and ``J_\\mathrm{ex}`` are the currents inside the radii
``r_\\mathrm{in}`` and ``r_\\mathrm{ex}``, and the record is the one above. A solid's
expression at an order is that of [Schelkunoff1934](@cite). The expressions of a finite wall
at an order are the solution of the boundary-value problem of equations (5) and (6) of
[Kane1995](@cite):

```math
\\begin{aligned}
z_\\mathrm{inner} &= \\frac{s \\mu}{2 \\pi w_\\mathrm{in}}
    \\frac{K_n(w_\\mathrm{in}) I_n'(w_\\mathrm{ex})
    - I_n(w_\\mathrm{in}) K_n'(w_\\mathrm{ex})}{D_n}, \\\\
z_\\mathrm{outer} &= \\frac{s \\mu}{2 \\pi w_\\mathrm{ex}}
    \\frac{I_n'(w_\\mathrm{in}) K_n(w_\\mathrm{ex})
    - K_n'(w_\\mathrm{in}) I_n(w_\\mathrm{ex})}{D_n}, \\\\
z_\\mathrm{transfer} &= \\frac{s \\mu}{2 \\pi w_\\mathrm{in} w_\\mathrm{ex} D_n}, \\\\
z_\\mathrm{outer}\\big|_{r_\\mathrm{in}=0} &= \\frac{s \\mu}{2 \\pi}
    \\frac{I_n(w_\\mathrm{ex})}{w_\\mathrm{ex} I_n'(w_\\mathrm{ex})},
\\end{aligned}
```

where
``D_n = I_n'(w_\\mathrm{in}) K_n'(w_\\mathrm{ex}) - K_n'(w_\\mathrm{in}) I_n'(w_\\mathrm{ex})``
and the primes denote derivatives in the argument.

The implementation calls the Bessel functions at the orders 0 and 1 alone, and it forms the
ratios of consecutive orders ``\\alpha_k(w) = w I_{k-1}(w) / I_k(w)`` and
``\\beta_k(w) = w K_{k-1}(w) / K_k(w)``. It runs ``\\alpha_k = 2k + w^2/\\alpha_{k+1}``
backward from ``\\alpha_N = 2N`` and ``\\beta_{k+1} = w^2/(\\beta_k + 2k)`` forward from
``\\beta_1``. The backward run starts at ``N = n + \\lceil |w_\\mathrm{ex}| \\rceil + 80``,
from the nominal value of ``|w_\\mathrm{ex}|``. The orders 0 and 1, in their exponentially scaled
form, give the products
``\\lambda_n = I_n(w_\\mathrm{in}) K_n(w_\\mathrm{ex}) / (I_n(w_\\mathrm{ex}) K_n(w_\\mathrm{in}))``
and ``\\Pi_n = I_n(w_\\mathrm{ex}) K_n(w_\\mathrm{in})``. With
``\\alpha_\\mathrm{in} = \\alpha_n(w_\\mathrm{in})``,
``\\alpha_\\mathrm{ex} = \\alpha_n(w_\\mathrm{ex})``,
``\\beta_\\mathrm{in} = \\beta_n(w_\\mathrm{in})`` and
``\\beta_\\mathrm{ex} = \\beta_n(w_\\mathrm{ex})``, the implementation evaluates

```math
\\begin{aligned}
\\mathrm{den} &= (\\beta_\\mathrm{in} + n)(\\alpha_\\mathrm{ex} - n)
    - \\lambda_n (\\alpha_\\mathrm{in} - n)(\\beta_\\mathrm{ex} + n), \\\\
z_\\mathrm{inner} &= \\frac{s \\mu}{2 \\pi}
    \\frac{(\\alpha_\\mathrm{ex} - n) + \\lambda_n (\\beta_\\mathrm{ex} + n)}{\\mathrm{den}}, \\\\
z_\\mathrm{outer} &= \\frac{s \\mu}{2 \\pi}
    \\frac{(\\beta_\\mathrm{in} + n) + \\lambda_n (\\alpha_\\mathrm{in} - n)}{\\mathrm{den}}, \\\\
z_\\mathrm{transfer} &= \\frac{s \\mu}{2 \\pi \\Pi_n \\mathrm{den}}, \\\\
z_\\mathrm{outer}\\big|_{r_\\mathrm{in}=0} &= \\frac{s \\mu}{2 \\pi (\\alpha_\\mathrm{ex} - n)}.
\\end{aligned}
```

### Callable

`formula(r_in, r_ex, ρ, μ_r, s)` returns the record `(; inner, outer, transfer)`, the
impedances of the inner and outer surfaces and the transfer impedance between them
\\[Ω/m\\].

- `r_in`: inner radius \\[m\\].
- `r_ex`: outer radius, in \\[m\\].
- `ρ`: resistivity \\[Ω·m\\].
- `μ_r`: relative permeability \\[dimensionless\\].
- `s`: complex frequency, `j 2π f` \\[1/s\\].

The radii are finite with `0 ≤ r_in < r_ex`, `ρ` and `μ_r` are positive and finite, and `s` is
finite and nonzero. A computation refuses a terminal whose conductor regions all have infinite
resistivity.

`formula(n, r_in, r_ex, ρ, μ_r, s)` returns that record at the harmonic order `n`
\\[Ω/m\\]. The order `n` \\[dimensionless\\] is a nonnegative integer, and the other
arguments are those of the first form. At `n = 0` it returns
`formula(r_in, r_ex, ρ, μ_r, s)`.

## `Formula{:wedepohl1973}`

Closed forms of the surface impedances of a thin tube, equations (4) to (6) of
[Wedepohl1973](@cite):

```math
\\begin{aligned}
z_\\mathrm{inner} &= \\frac{\\rho m}{2 \\pi r_\\mathrm{in}} \\coth(m \\Delta)
    - \\frac{\\rho}{2 \\pi r_\\mathrm{in} (r_\\mathrm{in} + r_\\mathrm{ex})}, \\\\
z_\\mathrm{outer} &= \\frac{\\rho m}{2 \\pi r_\\mathrm{ex}} \\coth(m \\Delta)
    + \\frac{\\rho}{2 \\pi r_\\mathrm{ex} (r_\\mathrm{in} + r_\\mathrm{ex})}, \\\\
z_\\mathrm{transfer} &= \\frac{\\rho m}{\\pi (r_\\mathrm{in} + r_\\mathrm{ex})
    \\sinh(m \\Delta)},
\\end{aligned}
```

where ``\\Delta = r_\\mathrm{ex} - r_\\mathrm{in}`` is the thickness of the tube \\[m\\]. The
other symbols are as in `Formula{:schelkunoff1934}`. At zero frequency each of the three
impedances equals the direct-current resistance of the annulus. The paper states good accuracy
for ``\\Delta / (r_\\mathrm{in} + r_\\mathrm{ex}) < 1/8``, and a computation warns for a thicker
tube. The closed forms are those of a net current, at order 0 alone. A computation at a higher
`order` warns, and refuses the members of an interior.

A solid conductor, with ``r_\\mathrm{in} = 0``, has no inner surface: its inner and transfer
impedances are zero, and its outer impedance is equation (3) of the paper:

```math
z_\\mathrm{outer} = \\frac{\\rho m}{2 \\pi r_\\mathrm{ex}} \\coth(0.777 m r_\\mathrm{ex})
    + \\frac{0.356 \\rho}{\\pi r_\\mathrm{ex}^2}.
```

### Callable

`formula(r_in, r_ex, ρ, μ_r, s)` returns the record `(; inner, outer, transfer)`, the
impedances of the inner and outer surfaces and the transfer impedance between them
\\[Ω/m\\].

- `r_in`: inner radius \\[m\\].
- `r_ex`: outer radius, in \\[m\\].
- `ρ`: resistivity \\[Ω·m\\].
- `μ_r`: relative permeability \\[dimensionless\\].
- `s`: complex frequency, `j 2π f` \\[1/s\\].
"""
module InternalImpedance

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults, validate
using ...Commons: FormulationOptions, nominal
using ..Engine: AbstractInternalImpedance, LineParametersProblem
using ..Engine: μ₀
using ...Materials: Material
import ..Engine: domain
using ..Engine: impedance!
using IntervalSets: (..)
using SpecialFunctions: besselix, besselkx

export Formula, FORMULAS

"Identifiers of the internal-impedance formulas."
const FORMULAS = (:schelkunoff1934, :wedepohl1973)

"""
$(TYPEDEF)

Represent a formula for the surface impedances of a tubular conductor.
Its options form a record of type `FormulationOptions{Formula{ID}}`, with `ID` the formula
identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, O <: FormulationOptions} <: AbstractInternalImpedance
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
"""
Formula(::Val{ID}; kw...) where {ID} = Formula{ID}(FormulationOptions(Formula{ID}; kw...))

# An identifier without a formula has no defaults.
function defaults(::Type{Formula{ID}}) where {ID}
    throw(ArgumentError(
        "unknown internal-impedance formula :$ID; the formulas are $FORMULAS"))
end

# A computation refuses a terminal whose conductor regions all have infinite resistivity. Its
# tube has no conductance, and a material admits an infinite resistivity.
function validate(problem::LineParametersProblem, formula::Formula)
    (; designs, regions, terminals) = problem.system
    for row in terminals
        rho = [region.source.material.rho
               for region in regions[row.regions]
               if region.source.material isa Material{<:Real, :conductor}]
        any(isfinite, rho) && continue
        design = designs[row.cable]
        throw(DomainError(rho,
            "the internal-impedance formula :$(nameof(formula)) computes the surfaces of a tube of finite resistivity; the terminal $(design.terminals[row.terminal].name) of the cable $(design.id) has no region of finite resistivity"))
    end
    return problem
end

# A formula gives the surfaces at the orders it has a method for. The recipe of the interiors
# above order 0 refuses a formula that gives order 0 alone.
function validate(formula::Formula, ::typeof(impedance!))
    throw(ArgumentError(
        "the formula :$(nameof(formula)) gives the surface impedances at order 0 alone; the coaxial formulation computes the members of an interior above order 0 with a formula that gives every order, or at order = 0"))
end

include("schelkunoff1934.jl")
include("wedepohl1973.jl")

end
