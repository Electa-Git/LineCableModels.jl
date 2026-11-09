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

### Callable

`formula(r_in, r_ex, ρ, μ_r, s)` returns the record `(; inner, outer, transfer)`, the
impedances of the inner and outer surfaces and the transfer impedance between them
\\[Ω/m\\].

- `r_in`: inner radius \\[m\\].
- `r_ex`: outer radius, in \\[m\\].
- `ρ`: resistivity \\[Ω·m\\].
- `μ_r`: relative permeability \\[dimensionless\\].
- `s`: complex frequency, `j 2π f` \\[1/s\\].

It throws `DomainError` unless the radii are finite with `0 ≤ r_in < r_ex`, `ρ` and `μ_r`
are positive and finite, and `s` is finite and nonzero.

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
tube.

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
import ...Commons: defaults
using ...Commons: FormulationOptions
using ..Engine: AbstractInternalImpedance
using ..Engine: μ₀
import ..Engine: domain
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

include("schelkunoff1934.jl")
include("wedepohl1973.jl")

end
