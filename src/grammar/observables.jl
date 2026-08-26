function _request_identity(request, supported::Tuple)
    request isa Function && return request
    request isa Tuple && !isempty(request) || throw(
        ArgumentError("observable requests must be selector functions or nonempty tuples"),
    )
    pair = length(request) >= 2 ? (request[1], request[2]) : nothing
    return pair !== nothing && pair in supported ? pair : first(request)
end

function _observable_declaration(source)
    applicable(observables, typeof(source)) || throw(ArgumentError(
        "$(typeof(source)) does not declare observable selectors",
    ))
    supported = observables(typeof(source))
    supported isa Tuple || throw(
        ArgumentError("observables($(typeof(source))) must return a tuple of selectors"),
    )
    return supported
end

"""
$(TYPEDSIGNATURES)

Validate named scientific requests and keyed display-unit overrides against
the observable declaration for `source`.

# Arguments

- `source`: Result or retained product that owns observable selectors.
- `requests`: Named selector functions or selector tuples.
- `unit_overrides`: Display-unit overrides keyed by fields in `requests`.

# Returns

- A tuple containing the normalized selector identity for each request.

# Errors

- Throws `ArgumentError` when the source has no observable declaration, the
  declaration is malformed, a request is unsupported, or an override key has
  no corresponding request.
"""
function validate_observables(
        source,
        requests::NamedTuple,
        unit_overrides::NamedTuple = (;)
)
    supported = _observable_declaration(source)
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

function _quantity(request)
    request isa Function && return quantity(request)
    request isa Tuple && !isempty(request) || throw(
        ArgumentError("scientific requests must be selector functions or nonempty tuples"),
    )
    transformed = length(request) >= 2 && request[2] isa Function
    return transformed ? quantity(request[1], request[2]) : quantity(first(request))
end

function _override_candidates(key, request)
    identity = request isa Function ? request :
               length(request) >= 2 && request[2] isa Function ?
               (request[1], request[2]) : first(request)
    names = identity isa Function ? (nameof(identity),) :
            identity isa Tuple && first(identity) isa Function ?
            (nameof(first(identity)),) : ()
    return key === nothing ? (request, identity, names...) :
           (key, request, identity, names...)
end

function _unit_override(overrides, key, request)
    overrides === nothing && return nothing
    overrides isa Union{Symbol, UnitExpr} && return overrides
    overrides isa Union{NamedTuple, AbstractDict} || throw(ArgumentError(
        "unit overrides must be a prefix, UnitExpr, keyed collection, or nothing",
    ))
    for candidate in _override_candidates(key, request)
        if overrides isa NamedTuple
            candidate isa Symbol && haskey(overrides, candidate) &&
                return overrides[candidate]
        elseif haskey(overrides, candidate)
            return overrides[candidate]
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Resolve display-unit targets aligned with a tuple or named tuple of scientific
requests.

# Arguments

- `requests`: Selector functions or selector tuples.
- `result_basis`: `:pul` or `:total`.

# Keywords

- `length_prefix`: Metric prefix applied to per-length denominators.
- `overrides`: A global prefix or `UnitExpr`, a keyed collection of those values,
  or nothing.

# Returns

- A tuple or named tuple of concrete `UnitExpr` values with the same shape as
  `requests`.
"""
function unit_targets(
        requests::NamedTuple,
        result_basis::Symbol;
        length_prefix::Symbol = :kilo,
        overrides = nothing
)
    targets = map(keys(requests), values(requests)) do key, request
        display_unit(
            _quantity(request),
            result_basis,
            _unit_override(overrides, key, request);
            length_prefix
        )
    end
    return NamedTuple{keys(requests)}(targets)
end

function unit_targets(
        requests::Tuple,
        result_basis::Symbol;
        length_prefix::Symbol = :kilo,
        overrides = nothing
)
    return map(requests) do request
        display_unit(
            _quantity(request),
            result_basis,
            _unit_override(overrides, nothing, request);
            length_prefix
        )
    end
end

"""
$(TYPEDSIGNATURES)

Detach an observed value from its result storage and apply a display-unit scale
factor. Structured result owners extend this method for their published value
types.
"""
detach(value::Number, factor) = value * factor
detach(values::AbstractArray, factor) = map(value -> value * factor, values)

function _publish_observable(source, request, identity, override)
    scientific_quantity = _quantity(identity)
    native = native_unit(scientific_quantity, basis(source))
    displayed = display_unit(scientific_quantity, basis(source), override)
    factor = scale_factor(native, displayed)
    detached = detach(_observe_request(source, request), factor)
    return (; values = detached, quantity = scientific_quantity, unit = displayed)
end

function observables(
        source,
        requests::NamedTuple;
        units::NamedTuple = (;)
)
    identities = validate_observables(source, requests, units)
    payloads = map(keys(requests), values(requests), identities) do key, request, identity
        override = haskey(units, key) ? getproperty(units, key) : nothing
        _publish_observable(source, request, identity, override)
    end
    return NamedTuple{keys(requests)}(payloads)
end
