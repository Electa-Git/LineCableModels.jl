"""
    LineCableModels.Engine.InsulationImpedance

Compute the magnetic series-impedance term across an annular dielectric.

# Formulas

## `Formula{:ametani1980}`

Magnetic term of the series impedance across an annular dielectric layer, from the field
between its radii [Ametani1980, Ametani2015b](@cite):

```math
z = \\frac{s \\mu_0 \\mu_r}{2 \\pi} \\ln \\frac{r_\\mathrm{ex}}{r_\\mathrm{in}},
```

where ``s`` is the complex frequency \\[1/s\\], ``r_\\mathrm{in}`` and
``r_\\mathrm{ex}`` are the inner and outer radii of the layer \\[m\\], and ``\\mu_r`` is its
relative permeability \\[dimensionless\\]. The constant ``\\mu_0`` is
[`vacuum_permeability`](@ref). The term is zero when ``r_\\mathrm{in}`` is zero or equals
``r_\\mathrm{ex}``.

### Callable

`formula(r_in, r_ex, μ_r, s)` returns ``z`` \\[Ω/m\\].

- `r_in`: inner radius \\[m\\].
- `r_ex`: outer radius, in \\[m\\].
- `μ_r`: relative permeability \\[dimensionless\\].
- `s`: complex frequency, `j 2π f` \\[1/s\\].
"""
module InsulationImpedance

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults
using ...Commons: FormulationOptions
using ..Engine: AbstractInsulationImpedance
using ..Engine: vacuum_permeability

export Formula, FORMULAS

"Identifiers of the insulation-impedance formulas."
const FORMULAS = (:ametani1980,)

"""
$(TYPEDEF)

Represent a magnetic series-impedance formula for an annular dielectric.
Its options form a record of type `FormulationOptions{Formula{ID}}`, with `ID` the formula
identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, O <: FormulationOptions} <: AbstractInsulationImpedance
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
        "unknown insulation-impedance formula :$ID; the formulas are $FORMULAS"))
end

include("ametani1980.jl")

end
