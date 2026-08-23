"""
    LineCableModels.Grammar

Own the package-local calculation, result, and parameter-space vocabulary used
by LineCableModels. Scientific problem and result types specialize this grammar
in their concept-owning modules.
"""
module Grammar

export AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult
export AbstractParametricResult, AbstractUncertaintyResult
export FormulationOptions, ComputationOptions
export compute, observables, primitives, preprocess
export Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid
export AbstractGrid, AbstractUncertainGrid, UncertainValue
export Gridspace, Configuration, configurations, materialize
export has_uncertainty, configuration_manifest, nominal, standard_uncertainty
export @gridspace, @relax
export Combinatorial, LinearError, MonteCarlo, ParametricProblem
export ParametricResult, LinearErrorResult, MonteCarloResult
export CalculationManifest, ConfigurationFailure, SampleSummary, HistogramDensity, RLCG
export result, statistics, samples, histograms, uncertain_value, manifest

using DocStringExtensions: SIGNATURES, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES

"""
$(TYPEDEF)

Abstract root for complete LineCableModels calculation inputs.
"""
abstract type AbstractProblemDefinition end

"""
$(TYPEDEF)

Abstract root for scientific and higher-order calculation selections.
"""
abstract type AbstractFormulation end

"""
$(TYPEDEF)

Abstract root for completed LineCableModels calculation results.
"""
abstract type AbstractProblemResult end

"""
$(TYPEDEF)

Abstract root for deterministic composite results whose primitive scientific
result type is `T`.
"""
abstract type AbstractParametricResult{T} <: AbstractProblemResult end

"""
$(TYPEDEF)

Abstract root for uncertainty results whose primitive scientific result type
is `T`.
"""
abstract type AbstractUncertaintyResult{T} <: AbstractProblemResult end

"Alias for validated, formulation-owned named-tuple options."
const FormulationOptions = NamedTuple
"Alias for validated, execution-owned named-tuple options."
const ComputationOptions = NamedTuple

"""
$(SIGNATURES)

Calculate a completed result from an explicit problem and formulation.

Concrete methods validate and normalize their supported `options`. Unsupported
problem/formulation pairs fail through ordinary Julia dispatch.
"""
function compute end

"""
$(SIGNATURES)

Return an immutable named tuple containing the scientific quantities published
by a completed result. Array-valued observables use detached storage so that
presentation and reporting cannot mutate the result.
"""
function observables end

"""
$(SIGNATURES)

Convert a completed result through an explicitly selected receiving
formulation. LineCableModels defines no implicit or zero-argument conversion.
"""
function primitives end

"""
$(SIGNATURES)

Construct a subsequent LineCableModels problem through explicitly selected
mathematics. Unsupported calculation orderings have no fallback method.
"""
function preprocess end

"Return the nominal value of a deterministic or uncertain quantity."
nominal(value) = value
nominal(value::Complex) = complex(nominal(real(value)), nominal(imag(value)))
nominal(values::AbstractArray) = nominal.(values)

"Return the standard uncertainty of a quantity; deterministic numbers return zero."
standard_uncertainty(value::Number) = zero(nominal(value))
standard_uncertainty(::Any) = 0.0

include("grammar/grid.jl")
include("grammar/gridspace.jl")
include("grammar/macros.jl")
include("grammar/computation.jl")

end
