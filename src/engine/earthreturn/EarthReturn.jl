"""
    LineCableModels.Engine.EarthReturn

The formulas of the earth-return impedance and potential coefficients of the cables. A
formula gives the impedance, the potential coefficients or both, by the methods it defines.
"""
module EarthReturn

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults, validate
using ...Commons: FormulationOptions
using ...DataModel: support
using ..Engine: AbstractEarthReturn, CoaxialInput, LineParametersProblem, SpectralIntegral,
                integrate, nominal, ray!, vacuum_permittivity, vacuum_permeability
import ..Engine: impedance!, potential!, prepare!, initialize_buffers
using LinearAlgebra: ldiv!, lu!
using SpecialFunctions: besselix, besselkx, besseljx, besselk

export Formula, FORMULAS

"Identifiers of the earth-return formulas."
const FORMULAS = (:unified,)

"""
$(TYPEDEF)

An earth-return formula, named by its identifier `ID`. Its options form a record
of type `FormulationOptions{Formula{ID}}`.

$(TYPEDFIELDS)
"""
struct Formula{ID, O <: FormulationOptions} <: AbstractEarthReturn
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
        "unknown earth-return formula :$ID; the formulas are $FORMULAS"))
end

# A formula gives the quantities it has a method for, and a slot refuses a formula that does
# not give its quantity.
function validate(formula::Formula{ID}, verb) where {ID}
    throw(DomainError(formula, "the earth-return formula :$ID gives no $(nameof(verb))"))
end

# Formulas of one identifier are equal when their options are equal.
Base.:(==)(a::Formula{ID}, b::Formula{ID}) where {ID} = a.options == b.options

include("unified.jl")

end
