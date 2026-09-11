"""
$(TYPEDEF)

Supertype for complete LineCableModels calculation inputs.
"""
abstract type AbstractProblemDefinition end

"""
$(TYPEDEF)

Supertype for scientific and higher-order calculation selections.
"""
abstract type AbstractFormulation end

"""
$(TYPEDEF)

Supertype for completed LineCableModels calculation results.
"""
abstract type AbstractProblemResult end

"""
$(TYPEDEF)

Supertype for direct results owned by LineCableModels computations.
"""
abstract type AbstractCoreResult <: AbstractProblemResult end

"""
$(TYPEDEF)

Supertype for completed finite collections whose element type is `T`.
"""
abstract type AbstractResultSpace{T} <: AbstractProblemResult end

"""
$(TYPEDEF)

Supertype for deterministic result spaces whose element type is `T`.
"""
abstract type AbstractParametricResult{T} <: AbstractResultSpace{T} end

"""
$(TYPEDEF)

Supertype for uncertainty result spaces whose element type is `T`.
"""
abstract type AbstractUncertaintyResult{T} <: AbstractResultSpace{T} end

"Alias identifying a formulation-owned named-tuple option record."
const FormulationOptions = NamedTuple
"Alias identifying an execution-owned named-tuple option record."
const ComputationOptions = NamedTuple
"Alias identifying an immutable computation-owned supplemental-output record."
const ComputationDetails = NamedTuple
