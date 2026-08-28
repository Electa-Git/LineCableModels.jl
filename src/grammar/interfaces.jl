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

function _observe_macro_parts(request)
    valid_request = request isa Expr &&
                    request.head === :ref &&
                    length(request.args) >= 2
    valid_request || throw(ArgumentError(
        "@observe expects indexed `accessor[...]` or `(accessor, transform)[...]`; " *
        "got `$(request)`.",
    ))

    selector = first(request.args)
    indices = request.args[2:end]
    if selector isa Expr && selector.head === :tuple
        length(selector.args) == 2 || throw(ArgumentError(
            "transformed @observe syntax requires exactly `(accessor, transform)`; " *
            "got `$(selector)`.",
        ))
        return selector.args[1], selector.args[2], indices
    end
    return selector, nothing, indices
end

"""
    @observe accessor[indices...]
    @observe (accessor, transform)[indices...]

Construct a plain observable-request tuple without reading a result.
"""
macro observe(request)
    selector, transform, indices = _observe_macro_parts(request)
    parts = transform === nothing ?
            (selector, indices...) :
            (selector, transform, indices...)
    return Expr(:tuple, map(esc, parts)...)
end

"""
    @observe source accessor[indices...]
    @observe source (accessor, transform)[indices...]

Expand indexed observable syntax into an immediate [`observe`](@ref) call. The
Indices follow the selected observable. Ordinary line quantities use row,
column, and sample indices; diagonal transforms use mode and sample indices.
"""
macro observe(source, request)
    selector, transform, indices = _observe_macro_parts(request)
    parts = transform === nothing ?
            (source, selector, indices...) :
            (source, selector, transform, indices...)
    escaped = map(esc, parts)
    return :(observe($(escaped...)))
end

"""
$(SIGNATURES)

Publish explicitly requested scientific values for presentation or reporting.

`observables(::Type{T})` declares the selectors supported by `T`.
`observables(source, requests; units, clip)` returns one detached payload for
each positional request. Every payload contains only `values`, `quantity`, and
`unit`. `units` is empty or positionally aligned with `requests`; `clip`
controls display-residue cleanup and defaults to `true`.
"""
function observables end
