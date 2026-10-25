"""
    LineCableModels.Engine.SemiconAdmittance

The admittivity laws of semiconducting layers.
"""
module SemiconAdmittance

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults
using ...Commons: FormulationOptions
using ..Engine: AbstractSemiconAdmittance
import ..Engine: admittivity
using ..Engine: vacuum_permittivity
using ...Materials: Material

export Formula, FORMULAS

"Identifiers of the semicon-admittance formulas."
const FORMULAS = (:lossless,)

"""
$(TYPEDEF)

A semicon-admittance formula, named by its identifier `ID`. Its options form a record
of type `FormulationOptions{Formula{ID}}`.

$(TYPEDFIELDS)
"""
struct Formula{ID, O <: FormulationOptions} <: AbstractSemiconAdmittance
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
        "unknown semicon-admittance formula :$ID; the formulas are $FORMULAS"))
end

include("lossless.jl")

end
