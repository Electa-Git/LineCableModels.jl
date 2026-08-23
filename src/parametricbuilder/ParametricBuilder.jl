module ParametricBuilder

export Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid
export AbstractGrid, AbstractUncertainGrid, UncertainValue
export Gridspace, Configuration, configurations, materialize
export has_uncertainty, configuration_manifest, nominal, standard_uncertainty
export @gridspace, @relax

export Material, Conductor, Insulator, CableBuilder
export at, trifoil, hflat, vflat, Earth, SystemBuilder
export WireEstimate, make_stranded, make_screened

using DocStringExtensions: SIGNATURES, TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
using Random
import ..LineCableModels: add!, nominal, standard_uncertainty, maxfill
import ..Grammar
import ..Grammar: Gridspace, configurations, has_uncertainty
using ..Grammar:
                 Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid,
                 AbstractGrid, AbstractUncertainGrid, UncertainValue,
                 RelativeUncertainty, AbsoluteUncertainty, _gridspace_axis,
                 Configuration, materialize, configuration_manifest,
                 @gridspace, @relax
import ..Materials
import ..Materials: Material
import ..DataModel
import ..EarthProps
import ..Engine

abstract type _AbstractDefinition{Target} end

target_type(::Type{<:_AbstractDefinition{Target}}) where {Target} = Target
target_type(definition::_AbstractDefinition) = target_type(typeof(definition))

Grammar._gridspace_axis(definition::_AbstractDefinition) = definition
Grammar._axis_cases(definition::_AbstractDefinition) = configurations(Gridspace(definition))
function has_uncertainty(definition::_AbstractDefinition)
    any(has_uncertainty, configurations(Gridspace(definition)))
end

include("materialdefinition.jl")
include("cablebuilderdefinition.jl")
include("positiondefinition.jl")
include("systembuilderdefinition.jl")

include("wirepatterns/WirePatterns.jl")
using .WirePatterns: WireEstimate, make_stranded, make_screened

end
