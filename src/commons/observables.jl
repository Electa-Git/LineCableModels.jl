"""
$(TYPEDSIGNATURES)

Return the function-valued selector prefix encoded by one observable
request. Positional indices are omitted from the result.
"""
function request_identity(request)
    request isa Function && return request
    request isa Tuple && !isempty(request) || throw(ArgumentError(
        "observable requests must be selector functions or nonempty tuples",
    ))
    count = findfirst(value -> !(value isa Function) || value isa Colon, request)
    count = count === nothing ? length(request) : count - 1
    count > 0 || throw(ArgumentError("an observable request must begin with a selector function"))
    return count == 1 ? first(request) : request[1:count]
end

"""
$(TYPEDSIGNATURES)

Normalize a selector prefix to the retained component selector used for lookup.
The input excludes positional indices. [`request_identity`](@ref) extracts that
prefix from a complete request. Scientific owners may extend this method for
accessor aliases. The default leaves selectors unchanged. `Val` selector names
also support keyed display-unit overrides such as `:alpha` and `:beta`.
"""
normalize_observation_selector(selector) = selector

"""
$(TYPEDSIGNATURES)

Return the [`LineCableModels.Units.Quantity`](@ref) encoded by a scientific
request.
"""
function request_quantity(request)
    identity=request_identity(request)
    return identity isa Tuple ? quantity(identity...) : quantity(identity)
end

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

- `source`: value whose type declares the supported observable identities.
- `request`: selector function or plain observable request tuple.

# Returns

- A named tuple containing `identity`, `quantity`, and `indices`.

# Errors

- `ArgumentError`: the declaration is malformed or the request is unsupported.
"""
function observation_request(source, request)
    supported = _observable_declaration(source)
    identity = request_identity(request)
    identity in supported || throw(ArgumentError(
        "$(typeof(source)) does not publish selector $(repr(identity))",
    ))
    selector_count = identity isa Tuple ? length(identity) : 1
    indices = request isa Function ? () : request[(selector_count + 1):end]
    return (; identity, quantity = request_quantity(request), indices)
end

"""
$(TYPEDSIGNATURES)

Resolve an integer, range, vector, or colon selector against a dimension of
length `count`.

# Returns

- Concrete integer indices in selection order.

# Errors

- `ArgumentError`: the selector form is unsupported.
- `BoundsError`: at least one resolved index is outside `1:count`.
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

- `source`: result or retained product that owns observable selectors.
- `requests`: tuple of selector functions or selector tuples.
- `unit_overrides`: empty tuple or one display unit for each request.

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
    isempty(unit_overrides) || length(unit_overrides) == length(requests) ||
        throw(
            DimensionMismatch("display units must align with observable requests"),
        )
    allunique(requests) || throw(ArgumentError("observable requests must be distinct"))
    return map(request -> observation_request(source, request).identity, requests)
end

function _unit_override_keys(request, overrides)
    identity = request_identity(request)
    indices=request_indices(request)
    selector=identity isa Tuple ? first(identity) : identity
    names=selector isa Base.Fix2 ? () : selector isa Function ? (nameof(selector),) : ()
    prefix = identity isa Tuple && length(identity) > 2 ? (identity[1:2],) : ()
    component_selector=normalize_observation_selector(identity)
    indexed_aliases=isempty(indices) ? () : Tuple(key for key in keys(overrides) if
        key isa Tuple && !isempty(key) && first(key) isa Function &&
        isequal(normalize_observation_selector(request_identity(key)),component_selector) &&
        isequal(request_indices(key),indices) && !isequal(key,request))
    function_aliases=Tuple(key for key in keys(overrides) if key isa Function &&
        isequal(normalize_observation_selector(key),component_selector) && !isequal(key,identity))
    name_aliases=Tuple(key for key in keys(overrides) if key isa Symbol &&
        isequal(normalize_observation_selector(Val(key)),component_selector))
    bound=selector isa Base.Fix2 ? (selector,selector.f,nameof(selector.f)) : ()
    physical=identity isa Tuple && length(identity)>1 && identity[2] isa Function &&
        applicable(quantity,identity[2]) ? (identity[2],nameof(identity[2])) : ()
    return (request,indexed_aliases...,identity, prefix..., bound...,function_aliases...,name_aliases...,
        names...,physical...)
end

function _unit_override(overrides, request)
    overrides === nothing && return nothing
    overrides isa Union{Symbol, UnitExpr} && return overrides
    overrides isa Union{NamedTuple, AbstractDict} || throw(ArgumentError(
        "unit overrides must be a prefix, UnitExpr, keyed collection, or nothing",
    ))
    for override_key in _unit_override_keys(request,overrides)
        if overrides isa NamedTuple
            override_key isa Symbol && haskey(overrides, override_key) &&
                return overrides[override_key]
        elseif haskey(overrides, override_key)
            return overrides[override_key]
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Resolve display-unit targets aligned with a tuple of scientific requests.

# Arguments

- `requests`: selector functions or selector tuples.
- `result_basis`: `:pul` or `:total`.

# Keywords

- `length_prefix`: metric prefix applied to per-length denominators.
- `overrides`: a global prefix or `UnitExpr`, a keyed collection of those values or `nothing`.

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
            request_quantity(request),
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
detach(value::AbstractFloat, factor::Real) = value * oftype(value, factor)
detach(value::Complex{T}, factor::Real) where {T<:AbstractFloat} = value * T(factor)
detach(value::Missing, factor) = missing
detach(value::Nothing, factor) = nothing
detach(values::NamedTuple, factor) = map(value -> detach(value, factor), values)
detach(values::AbstractArray, factor) = map(value -> detach(value, factor), values)

"""
$(TYPEDSIGNATURES)

Resolve the declared physical reporting resolution for one scientific request.
Result owners extend this operation. The fallback leaves precision unspecified.

# Arguments

- `source`: result defining the requested values and native physical basis.
- `request`: an observable selector or indexed request.

# Keywords

- `atol`: optional absolute reporting cutoff in the requested native units.
- `frequencies`: optional frequency context \\[Hz\\] for standalone tensors.

# Returns

- A record containing `kind`, native `atol` and `unit`, and detached
  `unresolved` and `available` masks aligned with the requested values.
  An unassessed request returns `nothing` for its cutoff and masks. Reporting cutoffs
  are not certified numerical forward-error bounds.
"""
function observation_resolution(source, request; atol=nothing, frequencies=nothing)
    return (kind=:unassessed, atol=nothing, unit=nothing,
        unresolved=nothing, available=nothing)
end

"""
$(TYPEDSIGNATURES)

Construct detached atomic observations. For an ordinary collection, apply the
atomic constructor to each element.
"""
observables(source,requests::Tuple=();kwargs...) = ObservedResult(source,requests;kwargs...)
