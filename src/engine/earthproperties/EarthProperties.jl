"""
    LineCableModels.Engine.EarthProperties

Compute the resistivity, relative permittivity, and relative permeability of earth and air.

# Formulas

## `Formula{:constant}`

Properties of a layer of the earth or of the air independent of the frequency ``f``:

```math
\\rho(f) = \\rho, \\qquad \\varepsilon_r(f) = \\varepsilon_r, \\qquad \\mu_r(f) = \\mu_r,
```

where ``\\rho`` is the resistivity of the layer \\[Ω·m\\], and ``\\varepsilon_r`` and
``\\mu_r`` are its relative permittivity and permeability \\[dimensionless\\].
"""
module EarthProperties

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults
using ...Commons: FormulationOptions
using ..Engine: AbstractEarthProperties
import ..Engine: properties
using ...Earth: EarthLayer

export Formula, FORMULAS

"Identifiers of the earth-properties formulas."
const FORMULAS = (:constant,)

"""
$(TYPEDEF)

Represent a frequency law for the properties of earth and air.
Its options form a record of type `FormulationOptions{Formula{ID}}`, with `ID` the formula
identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, O <: FormulationOptions} <: AbstractEarthProperties
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
        "unknown earth-properties formula :$ID; the formulas are $FORMULAS"))
end

include("constant.jl")

end
