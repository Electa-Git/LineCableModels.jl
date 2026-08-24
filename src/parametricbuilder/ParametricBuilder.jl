module ParametricBuilder

export Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid
export AbstractGrid, AbstractUncertainGrid, UncertainValue
export Gridspace, configurations, materialize
export has_uncertainty, configuration_manifest, nominal, standard_uncertainty
export @gridspace, @relax

export Combinatorial, ParametricProblem, ParametricResult
export CalculationManifest, ConfigurationFailure
export result, manifest

export Material, Conductor, Insulator, CableBuilder
export at, trifoil, hflat, vflat, Earth, SystemBuilder
export WireEstimate, make_stranded, make_screened

using DocStringExtensions: SIGNATURES, TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
using Random
using SHA
import ..LineCableModels: add!, maxfill
import ..Grammar
import ..Grammar: compute, observables, nominal, standard_uncertainty
using ..Grammar:
                 AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult,
                 AbstractParametricResult, AbstractUncertaintyResult
import ..Materials
import ..Materials: Material
import ..DataModel
import ..EarthProps
import ..Engine

include("grid.jl")
include("gridspace.jl")
include("macros.jl")
include("results.jl")

abstract type _AbstractDefinition{Target} end

target_type(::Type{<:_AbstractDefinition{Target}}) where {Target} = Target
target_type(definition::_AbstractDefinition) = target_type(typeof(definition))

_gridspace_axis(definition::_AbstractDefinition) = definition
_axis_cases(definition::_AbstractDefinition) = configurations(Gridspace(definition))
function has_uncertainty(definition::_AbstractDefinition)
    any(has_uncertainty, configurations(Gridspace(definition)))
end

include("materialdefinition.jl")
include("cablebuilderdefinition.jl")
include("positiondefinition.jl")
include("systembuilderdefinition.jl")

include("compute.jl")

include("wirepatterns/WirePatterns.jl")
using .WirePatterns: WireEstimate, make_stranded, make_screened

end
