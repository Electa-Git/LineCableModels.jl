"""
    LineCableModels.Engine.InsulationAdmittance

Compute the admittivity of insulator layers.
"""
module InsulationAdmittance

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults
using ...Commons: FormulationOptions
using ..Engine: AbstractInsulationAdmittance
import ..Engine: admittivity
using ..Engine: vacuum_permittivity
using ...Materials: Material

export Formula, FORMULAS

"Identifiers of the insulation-admittance formulas."
const FORMULAS = (:lossless,)

"""
$(TYPEDEF)

Represent an admittivity law for an insulator layer.
Its options form a record of type `FormulationOptions{Formula{ID}}`, with `ID` the formula
identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, O <: FormulationOptions} <: AbstractInsulationAdmittance
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
        "unknown insulation-admittance formula :$ID; the formulas are $FORMULAS"))
end

include("lossless.jl")

end
