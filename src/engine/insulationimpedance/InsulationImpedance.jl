"""
    LineCableModels.Engine.InsulationImpedance

Compute the magnetic series-impedance term across an annular dielectric.
"""
module InsulationImpedance

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults
using ...Commons: FormulationOptions
using ..Engine: AbstractInsulationImpedance
import ..Engine: impedance
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
