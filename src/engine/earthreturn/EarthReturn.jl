"""
    LineCableModels.Engine.EarthReturn

Compute the earth-return impedances and potential coefficients of cables. A formula gives
the impedance, the potential coefficients, or both, by the methods it defines.
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

Represent an earth-return impedance and potential-coefficient formula.
Its options form a record of type `FormulationOptions{Formula{ID}}`, with `ID` the formula
identifier.

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

Build the formula named by `ID` with its declared defaults and supplied options.

# Arguments

- `Val{ID}`: formula identifier from [`FORMULAS`](@ref).

# Keywords

- `kw...`: options declared by [`defaults`](@ref) for `Formula{ID}`.

# Returns

- A [`Formula`](@ref) with a typed option record.

# Errors

- Throws `ArgumentError` for an unknown identifier or option.
- Throws `ArgumentError` for invalid integration controls.
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
