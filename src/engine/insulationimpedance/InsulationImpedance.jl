"""
    LineCableModels.Engine.InsulationImpedance

The formulas of the magnetic term of the series impedance across an annular dielectric.
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

An insulation-impedance formula, named by its identifier `ID`. Its options form a record
of type `FormulationOptions{Formula{ID}}`.

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

Build the formula `ID` with the options given as keywords, over the defaults of its type.

# Errors

- Throws `ArgumentError` for an identifier that is not one of `FORMULAS`, or for an
  option that the formula lacks.
"""
Formula(::Val{ID}; kw...) where {ID} = Formula{ID}(FormulationOptions(Formula{ID}; kw...))

# An identifier without a formula has no defaults.
function defaults(::Type{Formula{ID}}) where {ID}
    throw(ArgumentError(
        "unknown insulation-impedance formula :$ID; the formulas are $FORMULAS"))
end

include("ametani1980.jl")

end
