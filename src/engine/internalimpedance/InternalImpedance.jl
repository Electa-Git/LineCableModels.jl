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
The constant ``\\mu_0`` is [`vacuum_permeability`](@ref). The modified Bessel functions
``I_\\nu`` and ``K_\\nu`` have order ``\\nu``.

A solid conductor, with ``r_\\mathrm{in} = 0``, has no inner surface: its inner and transfer
impedances are zero, and its outer impedance is
``s \\mu I_0(w_\\mathrm{ex}) / (2 \\pi w_\\mathrm{ex} I_1(w_\\mathrm{ex}))``. The
implementation evaluates the Bessel functions in their exponentially scaled form.
"""
module InternalImpedance

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults
using ...Commons: FormulationOptions
using ..Engine: AbstractInternalImpedance
using ..Engine: vacuum_permeability
using SpecialFunctions: besselix, besselkx

export Formula, FORMULAS

"Identifiers of the internal-impedance formulas."
const FORMULAS = (:schelkunoff1934,)

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

end
