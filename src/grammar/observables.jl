function _request_identity(request, supported::Tuple)
    request isa Function && return request
    request isa Tuple && !isempty(request) || throw(
        ArgumentError("observable requests must be selector functions or nonempty tuples"),
    )
    pair = length(request) >= 2 ? (request[1], request[2]) : nothing
    return pair !== nothing && pair in supported ? pair : first(request)
end

"""
$(TYPEDSIGNATURES)

Return the selector or selector/transform pair encoded by one observable
request. Positional indices are omitted from the result.
"""
function request_identity(request)
    request isa Function && return request
    request isa Tuple && !isempty(request) || throw(ArgumentError(
        "observable requests must be selector functions or nonempty tuples",
    ))
    return length(request) >= 2 && request[1] isa Function &&
           request[2] isa Function && !(request[2] isa Colon) ?
           (request[1], request[2]) : first(request)
end

"""
$(TYPEDSIGNATURES)

Return the [`LineCableModels.Units.Quantity`](@ref) encoded by a scientific
request.
"""
request_quantity(request) = _quantity(request_identity(request))

"""
$(TYPEDSIGNATURES)

Return the positional indices encoded by a scientific request.
"""
function request_indices(request)
    request isa Function && return ()
    identity = request_identity(request)
    selector_count = identity isa Tuple ? length(identity) : 1
    return request[(selector_count + 1):end]
end

function _observable_declaration(source)
    supported = observables(typeof(source))
    supported isa Tuple || throw(
        ArgumentError("observables($(typeof(source))) must return a tuple of selectors"),
    )
    return supported
end

"""
$(TYPEDSIGNATURES)

Resolve one scientific request into its declared selector identity, physical
quantity, and positional indices.

# Arguments

- `source`: Value whose type declares the supported observable identities.
- `request`: Selector function or plain observable request tuple.

# Returns

- A named tuple containing `identity`, `quantity`, and `indices`.

# Errors

- `ArgumentError`: The declaration is malformed or the request is unsupported.
"""
function observation_request(source, request)
    supported = _observable_declaration(source)
    identity = _request_identity(request, supported)
    identity in supported || throw(ArgumentError(
        "$(typeof(source)) does not publish selector $(repr(identity))",
    ))
    selector_count = identity isa Tuple ? length(identity) : 1
    indices = request isa Function ? () : request[(selector_count + 1):end]
    return (; identity, quantity = _quantity(identity), indices)
end

"""
$(TYPEDSIGNATURES)

Resolve an integer, range, vector, or colon selector against a dimension of
length `count`.

# Returns

- Concrete integer indices in selection order.

# Errors

- `ArgumentError`: The selector form is unsupported.
- `BoundsError`: At least one resolved index is outside `1:count`.
"""
function observation_indices(selector, count::Integer)
    indices = if selector isa Colon
        collect(1:Int(count))
    elseif selector isa Integer
        [Int(selector)]
    elseif selector isa AbstractRange || selector isa AbstractVector
        collect(Int, selector)
    else
        throw(ArgumentError(
            "observable indices must be integers, ranges, vectors, or `:`",
        ))
    end
    all(index -> index in 1:count, indices) || throw(BoundsError(1:count, indices))
    return indices
end

"""
$(TYPEDSIGNATURES)

Return a plain observable request with `indices` and the selector identity from
`resolved`.
"""
function materialize_observation(resolved::NamedTuple, indices::Tuple)
    identity = resolved.identity
    prefix = identity isa Tuple ? identity : (identity,)
    return (prefix..., indices...)
end

"""
$(TYPEDSIGNATURES)

Validate positional scientific requests and aligned display-unit overrides
against the observable declaration for `source`.

# Arguments

- `source`: Result or retained product that owns observable selectors.
- `requests`: Tuple of selector functions or selector tuples.
- `unit_overrides`: Empty tuple or one display unit for each request.

# Returns

- A tuple containing the normalized selector identity for each request.

# Errors

- Throws `ArgumentError` when the declaration is malformed or a request is
  unsupported.
- Throws `DimensionMismatch` when the display-unit tuple is not aligned with
  the requests.
- A source that does not implement the declaration fails through ordinary
  method dispatch.
"""
function validate_observables(
        source,
        requests::Tuple,
        unit_overrides::Tuple = ()
)
    supported = _observable_declaration(source)
    isempty(unit_overrides) || length(unit_overrides) == length(requests) ||
        throw(
            DimensionMismatch("display units must align with observable requests"),
        )
    identities = map(request -> _request_identity(request, supported), requests)
    for identity in identities
        identity in supported || throw(
            ArgumentError("$(typeof(source)) does not publish selector $(repr(identity))"),
        )
    end
    return identities
end

_observe_request(source, request::Function) = observe(source, request)
_observe_request(source, request::Tuple) = observe(source, request...)

_quantity(request::Function) = quantity(request)

function _quantity(request::Tuple{F, Colon, Vararg}) where {F <: Function}
    return quantity(first(request))
end

function _quantity(request::Tuple{F, G, Vararg}) where {F <: Function, G <: Function}
    return quantity(request[1], request[2])
end

function _quantity(request::Tuple)
    isempty(request) && throw(
        ArgumentError("scientific requests must be selector functions or nonempty tuples"),
    )
    return quantity(first(request))
end

function _override_candidates(request)
    identity = request isa Function ? request :
               length(request) >= 2 && request[1] isa Function &&
               request[2] isa Function && !(request[2] isa Colon) ?
               (request[1], request[2]) : first(request)
    names = identity isa Function ? (nameof(identity),) :
            identity isa Tuple && first(identity) isa Function ?
            (nameof(first(identity)),) : ()
    return (request, identity, names...)
end

function _unit_override(overrides, request)
    overrides === nothing && return nothing
    overrides isa Union{Symbol, UnitExpr} && return overrides
    overrides isa Union{NamedTuple, AbstractDict} || throw(ArgumentError(
        "unit overrides must be a prefix, UnitExpr, keyed collection, or nothing",
    ))
    for candidate in _override_candidates(request)
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

Resolve display-unit targets aligned with a tuple of scientific requests.

# Arguments

- `requests`: Selector functions or selector tuples.
- `result_basis`: `:pul` or `:total`.

# Keywords

- `length_prefix`: Metric prefix applied to per-length denominators.
- `overrides`: A global prefix or `UnitExpr`, a keyed collection of those values,
  or nothing.

# Returns

- A tuple of concrete `UnitExpr` values aligned with `requests`.
"""
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
            _unit_override(overrides, request);
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

const _DISPLAY_CLIP_TOLERANCE = eps(Float64)

_clip_detached(value, ::Val{false}) = value
_clip_detached(values::AbstractArray, ::Val{false}) = values

function _clip_detached(value::Real, ::Val{true})
    isfinite(value) || return value
    return abs(value) <= _DISPLAY_CLIP_TOLERANCE ? zero(value) : value
end

_clip_detached(value::Complex, ::Val{true}) = value
_clip_detached(value::Missing, ::Val{true}) = value
function _clip_detached(values::AbstractArray, ::Val{true})
    return map(value -> _clip_detached(value, Val(true)), values)
end
_clip_detached(value, ::Val{true}) = value

"""
$(TYPEDSIGNATURES)

Detach an observed value, convert it by `factor`, and optionally replace
finite scalar display residue with exact zero. Structured value owners may add
narrow methods that preserve their constructor invariants.

# Arguments

- `value`: Observed scalar, array, or supported structured product.
- `factor`: Multiplicative native-to-display unit conversion.
- `clip`: Whether values no larger than machine epsilon after conversion are
  replaced by exact zero.

# Returns

- A detached value in the requested display unit.
"""
function detach(value, factor, clip::Bool)
    return _clip_detached(detach(value, factor), Val(clip))
end

function _publish_observable(source, request, identity, override, clip::Bool)
    scientific_quantity = _quantity(identity)
    native = native_unit(scientific_quantity, basis(source))
    displayed = display_unit(scientific_quantity, basis(source), override)
    factor = scale_factor(native, displayed)
    detached = detach(_observe_request(source, request), factor, clip)
    return (; values = detached, quantity = scientific_quantity, unit = displayed)
end

function observables(
        source,
        requests::Tuple;
        units::Tuple = (),
        clip::Bool = true
)
    identities = validate_observables(source, requests, units)
    overrides = isempty(units) ? ntuple(_ -> nothing, length(requests)) : units
    payloads = map(requests, identities, overrides) do request, identity, override
        _publish_observable(source, request, identity, override, clip)
    end
    return payloads
end
