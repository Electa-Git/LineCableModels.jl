"""
$(SIGNATURES)

Validate and normalise the options owned by a formulation type.

The implementation that owns `FormulationType` defines a method for
`Val(FormulationType)`. No broad fallback exists. An unregistered formulation
raises `MethodError`.

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

Validate and normalise the options owned by one computation.

The implementation that owns `OwnerType` defines a method for `Val(OwnerType)`.
`OwnerType` may identify a primitive solver or a composite calculation such as
a Gauntlet case. No broad fallback exists. An unregistered owner raises
`MethodError`.

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

Concrete methods validate and normalise their supported `options`. Unsupported
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

Convert a completed result through an explicitly selected mathematical
operation.

LineCableModels currently defines no methods. Unsupported operations raise
`MethodError`.
"""
function primitives end

"""
$(SIGNATURES)

Construct a subsequent problem through an explicitly selected mathematical
operation.

LineCableModels currently defines no methods. Unsupported operations raise
`MethodError`.
"""
function preprocess end
