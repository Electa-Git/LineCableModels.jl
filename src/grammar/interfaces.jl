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

function _request_identity(request, supported::Tuple)
    request isa Function && return request
    request isa Tuple && !isempty(request) || throw(
        ArgumentError("observable requests must be selector functions or nonempty tuples"),
    )
    pair = length(request) >= 2 ? (request[1], request[2]) : nothing
    return pair !== nothing && pair in supported ? pair : first(request)
end

function _validate_observable_requests(source, requests::NamedTuple, unit_overrides::NamedTuple)
    supported = observables(typeof(source))
    supported isa Tuple || throw(
        ArgumentError("observables($(typeof(source))) must return a tuple of selectors"),
    )
    unsupported_units = Tuple(key for key in keys(unit_overrides) if key ∉ keys(requests))
    isempty(unsupported_units) || throw(
        ArgumentError("unit overrides do not match requests: $(join(unsupported_units, ", "))"),
    )
    identities = map(request -> _request_identity(request, supported), values(requests))
    for identity in identities
        identity in supported || throw(
            ArgumentError("$(typeof(source)) does not publish selector $(repr(identity))"),
        )
    end
    return identities
end

_observe_request(source, request::Function) = observe(source, request)
_observe_request(source, request::Tuple) = observe(source, request...)

function _quantity(identity)
    return identity isa Tuple ? quantity(identity...) : quantity(identity)
end

_detach_and_scale(value::Number, factor) = value * factor
_detach_and_scale(values::AbstractArray, factor) = map(value -> value * factor, values)

function _publish_observable(source, request, identity, target)
    scientific_quantity = _quantity(identity)
    native = native_unit(scientific_quantity, basis(source))
    displayed = target === nothing ? display_unit(scientific_quantity, basis(source)) :
                target
    displayed isa UnitExpr ||
        throw(ArgumentError("observable unit overrides must be UnitExpr values"))
    factor = scale_factor(native, displayed)
    detached = _detach_and_scale(_observe_request(source, request), factor)
    return (; values = detached, quantity = scientific_quantity, unit = displayed)
end

function observables(
        source,
        requests::NamedTuple;
        units::NamedTuple = (;)
)
    identities = _validate_observable_requests(source, requests, units)
    payloads = map(keys(requests), values(requests), identities) do key, request, identity
        target = haskey(units, key) ? getproperty(units, key) : nothing
        _publish_observable(source, request, identity, target)
    end
    return NamedTuple{keys(requests)}(payloads)
end
