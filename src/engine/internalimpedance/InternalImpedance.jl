"""
    LineCableModels.Engine.InternalImpedance

The formulas of the surface impedances of a tubular conductor: the impedance of its inner
and outer surfaces and the transfer impedance between them.
"""
module InternalImpedance

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults
using ...Commons: FormulationOptions
using ..Engine: AbstractInternalImpedance

export Formula, FORMULAS

"Identifiers of the internal-impedance formulas."
const FORMULAS = ()

"""
$(TYPEDEF)

An internal-impedance formula, named by its identifier `ID`. Its options form a record
of type `FormulationOptions{Formula{ID}}`.

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

Build the formula `ID` with the options given as keywords, over the defaults of its type.

# Errors

- Throws `ArgumentError` for an identifier that is not one of `FORMULAS`, or for an
  option that the formula lacks.
"""
Formula(::Val{ID}; kw...) where {ID} = Formula{ID}(FormulationOptions(Formula{ID}; kw...))

# An identifier without a formula has no defaults.
function defaults(::Type{Formula{ID}}) where {ID}
    throw(ArgumentError(
        "unknown internal-impedance formula :$ID; the formulas are $FORMULAS"))
end

end
