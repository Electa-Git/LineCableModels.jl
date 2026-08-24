"""
    LineCableModels.Grammar

Own the shared calculation roots and action generics used by LineCableModels.
Scientific and higher-order types specialize this vocabulary in their
concept-owning modules.
"""
module Grammar

export AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult
export AbstractParametricResult, AbstractUncertaintyResult
export FormulationOptions, ComputationOptions
export formulation_options, computation_options
export compute, observables, primitives, preprocess
export nominal, standard_uncertainty

using DocStringExtensions: SIGNATURES, TYPEDEF

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

"Alias identifying a formulation-owned named-tuple option record."
const FormulationOptions = NamedTuple
"Alias identifying an execution-owned named-tuple option record."
const ComputationOptions = NamedTuple

"""
$(SIGNATURES)

Validate and normalize the options owned by a formulation type.

Backends extend this function with an owner token of the form
`Val(FormulationType)`. LineCableModels defines no broad fallback: a formulation
without an option grammar fails through ordinary Julia dispatch.

# Arguments

- `owner`: Formulation-type dispatch token.
- `options`: Mergeable named tuple supplied by the caller.

# Returns

- A formulation-owned [`FormulationOptions`](@ref) named tuple with a fixed set
  of keys for the selected owner.
"""
function formulation_options end

"""
$(SIGNATURES)

Validate and normalize the options owned by one computation or orchestrator.

Backends and orchestrators extend this function with an owner token of the form
`Val(OwnerType)`. LineCableModels defines no broad fallback: an owner without a
computation-option grammar fails through ordinary Julia dispatch.

# Arguments

- `owner`: Computation-owner dispatch token.
- `options`: Mergeable named tuple supplied by the caller.

# Returns

- A computation-owned [`ComputationOptions`](@ref) named tuple with a fixed set
  of outer keys for the selected owner.
"""
function computation_options end

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

Reserved extension point for converting completed results through explicitly
selected future mathematics. LineCableModels defines no methods, including no
zero-argument or broad fallback method.
"""
function primitives end

"""
$(SIGNATURES)

Reserved extension point for constructing a subsequent problem through
explicitly selected future mathematics. LineCableModels defines no methods;
unsupported calculation orderings fail through ordinary Julia dispatch.
"""
function preprocess end

"Return the nominal value of a deterministic or uncertain quantity."
nominal(value) = value
nominal(value::Complex) = complex(nominal(real(value)), nominal(imag(value)))
nominal(values::AbstractArray) = nominal.(values)

"Return the standard uncertainty of a quantity; deterministic numbers return zero."
standard_uncertainty(value::Number) = zero(nominal(value))
standard_uncertainty(::Any) = 0.0

end
