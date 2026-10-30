"""
    LineCableModels.Engine.Resistivity

Compute conductor resistivity as a function of temperature.

# Formulas

## `Formula{:linear}`

Resistivity of a conductor at the temperature ``\\theta``, by the linear correction of its
reference value [IEC60287](@cite):

```math
\\rho(\\theta) = \\rho_0 \\left[1 + \\alpha (\\theta - T_0)\\right],
```

where ``\\rho_0`` is the resistivity `rho` of the material at the reference temperature
``T_0``, the field `T0` \\[°C\\], in \\[Ω·m\\], ``\\alpha`` is the temperature coefficient
`alpha` \\[1/K\\], and ``\\theta`` is the conductor temperature \\[°C\\]. The law holds
within 150 K of ``T_0`` and returns an infinite resistivity unchanged.
"""
module Resistivity

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults
using ...Commons: FormulationOptions
using ..Engine: AbstractResistivity
using ...Materials: Material

export Formula, FORMULAS

"Identifiers of the resistivity formulas."
const FORMULAS = (:linear,)

"""
$(TYPEDEF)

Represent a temperature law for conductor resistivity.
Its options form a record of type `FormulationOptions{Formula{ID}}`, with `ID` the formula
identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, O <: FormulationOptions} <: AbstractResistivity
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
        "unknown resistivity formula :$ID; the formulas are $FORMULAS"))
end

include("linear.jl")

end
