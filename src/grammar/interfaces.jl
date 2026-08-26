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
`OwnerType` may identify a core solver or a composite calculation such as
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

Normalize supplemental output owned by one core or composite computation.

The computation owner defines a method for `Val(OwnerType)` and returns a
fixed-key [`ComputationDetails`](@ref) named tuple. No broad fallback exists;
an unregistered owner raises `MethodError`.
"""
function computation_details end

"""
$(SIGNATURES)

Return the typed supplemental output retained by a completed higher-order
result. Result owners define narrow methods beside their result containers.
"""
function details end

"""
$(SIGNATURES)

Calculate a completed result from an explicit problem and formulation.

Concrete methods validate and normalise their supported `options`. Unsupported
problem/formulation pairs fail through ordinary Julia dispatch.
"""
function compute end

"""
$(SIGNATURES)

Return native numerical values selected from a completed scientific result.

The selector and optional transform are function objects. Result owners define
the supported combinations beside their result representations.
"""
function observe end

"""
    @observe source accessor[i, j, k]

Expand indexed observable syntax into
`observe(source, accessor, i, j, k)`. The first two indices select a matrix
element and the third selects the sample dimension.
"""
macro observe(source, request)
    valid_request = request isa Expr &&
                    request.head === :ref &&
                    length(request.args) == 4
    valid_request || throw(ArgumentError(
        "@observe expects `accessor[i, j, k]`, for example " *
        "`@observe parameters Z[1, 1, :]`; got `$(request)`.",
    ))
    accessor, i, j, k = request.args
    return :(
        observe(
        $(esc(source)),
        $(esc(accessor)),
        $(esc(i)),
        $(esc(j)),
        $(esc(k))
    )
    )
end

"""
$(SIGNATURES)

Publish explicitly requested scientific values for presentation or reporting.

`observables(::Type{T})` declares the selectors supported by `T`.
`observables(source, requests; units)` returns one detached payload for each
named request. Every payload contains only `values`, `quantity`, and `unit`.
"""
function observables end
