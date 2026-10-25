"""
    LineCableModels.Engine.EarthProperties

The frequency laws of the resistivity, permittivity and permeability of the earth.
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

An earth-properties formula, named by its identifier `ID`. Its options form a record
of type `FormulationOptions{Formula{ID}}`.

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

Build the formula `ID` with the options given as keywords, over the defaults of its type.

# Errors

- Throws `ArgumentError` for an identifier that is not one of `FORMULAS`, or for an
  option that the formula lacks.
"""
Formula(::Val{ID}; kw...) where {ID} = Formula{ID}(FormulationOptions(Formula{ID}; kw...))

# An identifier without a formula has no defaults.
function defaults(::Type{Formula{ID}}) where {ID}
    throw(ArgumentError(
        "unknown earth-properties formula :$ID; the formulas are $FORMULAS"))
end

include("constant.jl")

end
