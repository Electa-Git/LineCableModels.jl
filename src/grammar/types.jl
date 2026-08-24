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

Supertype for deterministic composite results whose primitive scientific
result type is `T`.
"""
abstract type AbstractParametricResult{T} <: AbstractProblemResult end

"""
$(TYPEDEF)

Supertype for uncertainty results whose primitive scientific result type
is `T`.
"""
abstract type AbstractUncertaintyResult{T} <: AbstractProblemResult end

"Alias identifying a formulation-owned named-tuple option record."
const FormulationOptions = NamedTuple
"Alias identifying an execution-owned named-tuple option record."
const ComputationOptions = NamedTuple
