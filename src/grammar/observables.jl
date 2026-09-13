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
    identity = request_identity(request)
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
    isempty(unit_overrides) || length(unit_overrides) == length(requests) ||
        throw(
            DimensionMismatch("display units must align with observable requests"),
        )
    allunique(requests) || throw(ArgumentError("observable requests must be distinct"))
    return map(request -> observation_request(source, request).identity, requests)
end

_observe_request(source, request::Function) = observe(source, request)
_observe_request(source, request::Tuple) = observe(source, request...)

_quantity(request::Function) = quantity(request)

function _quantity(request::Tuple)
    identity = request_identity(request)
    return identity isa Tuple ? quantity(identity...) : quantity(identity)
end

function _override_candidates(request)
    identity = request_identity(request)
    names = identity isa Function ? (nameof(identity),) :
            identity isa Tuple && first(identity) isa Function ?
            (nameof(first(identity)),) : ()
    prefix = identity isa Tuple && length(identity) > 2 ? (identity[1:2],) : ()
    return (request, identity, prefix..., names...)
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
detach(value::AbstractFloat, factor::Real) = value * oftype(value, factor)
detach(value::Complex{T}, factor::Real) where {T<:AbstractFloat} = value * T(factor)
detach(value::Missing, factor) = missing
detach(values::AbstractArray, factor) = map(value -> detach(value, factor), values)

"""
$(TYPEDSIGNATURES)

Detach an observed value and convert it by `factor`. Without an observation
owner and physical quantity, no numerical-resolution threshold is inferred.
Structured value owners preserve their constructor and uncertainty invariants.

# Arguments

- `value`: Observed scalar, array, or supported structured product.
- `factor`: Multiplicative native-to-display unit conversion.
- `clip`: Retained call compatibility. Physical clipping is performed by
  [`observables`](@ref), before conversion, using [`observation_resolution`](@ref).

# Returns

- A detached value in the requested display unit.
"""
function detach(value, factor, clip::Bool)
    return detach(value, factor)
end

"""
$(TYPEDSIGNATURES)

Resolve the declared physical reporting resolution for one scientific request.
Result owners extend this operation; the fallback makes no precision claim.

# Arguments

- `source`: Result owning the requested values and native physical basis.
- `request`: An observable selector or indexed request.

# Keywords

- `atol`: Optional absolute reporting cutoff in the requested native units.
- `frequencies`: Optional frequency context \\[Hz\\] for standalone tensors.

# Returns

- A record containing `kind`, semantic `revision`, native `atol` and `unit`, and
  a detached `unresolved` mask aligned with the request. An unassessed request
  returns `nothing` for its cutoff, unit and mask. Declared reporting cutoffs
  are not certified numerical forward-error bounds.
"""
function observation_resolution(source, request; atol=nothing, frequencies=nothing)
    return (kind=:unassessed, revision=0, atol=nothing, unit=nothing, unresolved=nothing)
end

_resolved_observation(value, ::Nothing, phase) = value
_resolved_observation(value, unresolved::Bool, ::Val{false}) = unresolved ? zero(value) : value
_resolved_observation(value, unresolved::Bool, ::Val{true}) = unresolved ? missing : value
function _resolved_observation(values::AbstractArray, unresolved::AbstractArray, phase)
    return map((value, masked) -> _resolved_observation(value, masked, phase), values, unresolved)
end

function _publish_observable(source, request, identity, override, clip::Bool, atol, frequencies)
    scientific_quantity = _quantity(identity)
    native = native_unit(scientific_quantity, basis(source))
    displayed = display_unit(scientific_quantity, basis(source), override)
    factor = scale_factor(native, displayed)
    resolution = observation_resolution(source, request; atol, frequencies)
    values = _observe_request(source, request)
    phase = Val(identity isa Tuple && last(identity) === angle)
    resolved = clip ? _resolved_observation(values, resolution.unresolved, phase) : values
    detached = detach(resolved, factor)
    masked = resolution.unresolved
    unresolved_count = masked === nothing ? 0 : masked isa Bool ? Int(masked) : count(masked)
    return (
        observation=(; values=detached, quantity=scientific_quantity, unit=displayed),
        resolution=(; resolution.kind, resolution.revision, resolution.atol, resolution.unit,
            clip, unresolved_count),
    )
end

"""
$(TYPEDEF)

Hold detached scientific observations and their column-oriented table view.

The publication is the sole Tables.jl boundary for result observations. It
does not retain the source result and cannot reopen result storage.
"""
struct ObservationPublication{P <: Tuple, C <: NamedTuple, M <: NamedTuple}
    "Detached observations in request order."
    observations::P
    "Validated equal-length table columns."
    columns::C
    "Basis, row order, and quantity and unit metadata for each column."
    metadata::M

    function ObservationPublication(
            observations::P,
            columns::C,
            metadata::M
    ) where {P <: Tuple, C <: NamedTuple, M <: NamedTuple}
        lengths = map(length, values(columns))
        isempty(lengths) || all(==(first(lengths)), lengths) || throw(
            DimensionMismatch("observation publication columns must have equal lengths"),
        )
        keys(metadata) == (:basis, :row_order, :observation_columns) || throw(
            ArgumentError(
                "observation publication metadata must contain basis, row_order, and observation_columns",
            ),
        )
        return new{P, C, M}(observations, columns, metadata)
    end
end

Base.length(publication::ObservationPublication) = length(publication.observations)
Base.firstindex(publication::ObservationPublication) = firstindex(publication.observations)
Base.lastindex(publication::ObservationPublication) = lastindex(publication.observations)
Base.getindex(publication::ObservationPublication, index::Integer) =
    publication.observations[index]
Base.iterate(publication::ObservationPublication, state...) =
    iterate(publication.observations, state...)
Base.tail(publication::ObservationPublication) = Base.tail(publication.observations)

function Base.summary(io::IO, publication::ObservationPublication)
    rows = isempty(publication.columns) ? 0 : length(first(values(publication.columns)))
    print(io, "Observation publication with $rows rows")
end
function Base.show(io::IO, publication::ObservationPublication)
    rows = isempty(publication.columns) ? 0 : length(first(values(publication.columns)))
    print(io, "ObservationPublication(", rows, " rows × ", length(publication.columns), " columns)")
end
function Base.show(io::IO, ::MIME"text/plain", publication::ObservationPublication)
    show(io, publication)
end

#! explicit-imports: off
# Tables' interface functions are deliberately qualified protocol extensions;
# the package does not claim local ownership of those dependency bindings.
Tables.istable(::Type{<:ObservationPublication}) = true
Tables.columnaccess(::Type{<:ObservationPublication}) = true
Tables.columns(publication::ObservationPublication) = publication.columns
Tables.schema(publication::ObservationPublication) = Tables.schema(publication.columns)
Tables.columnnames(publication::ObservationPublication) = keys(publication.columns)
Tables.getcolumn(publication::ObservationPublication, index::Int) =
    getfield(publication.columns, index)
Tables.getcolumn(publication::ObservationPublication, name::Symbol) =
    getproperty(publication.columns, name)
#! explicit-imports: on

_publication_column(value::Number) = [value]
_publication_column(value::AbstractArray) = collect(vec(value))
_publication_column(value) = [value]

function _publication_names(observations::Tuple)
    names = map(payload -> Symbol(Units.symbol(payload.quantity)), observations)
    all(name -> !isempty(string(name)), names) || throw(ArgumentError(
        "every published quantity must define a nonempty table symbol",
    ))
    length(unique(names)) == length(names) || throw(ArgumentError(
        "published quantities must have distinct table symbols",
    ))
    return names
end

function _publication_contract(names::Tuple, observations::Tuple)
    records = map(observations) do payload
        (; quantity = payload.quantity, unit = payload.unit)
    end
    return NamedTuple{names}(records)
end

"""
$(TYPEDSIGNATURES)

Construct the detached column layout for one observation owner.

Result owners add methods when scientific coordinates such as row, column, or
frequency accompany the requested quantity columns.

# Arguments

- `source`: Result that owns the observations.
- `requests`: Positional scientific requests.
- `observations`: Detached values, quantities, and display units in request
  order.
- `options`: Display options used to detach the values.

# Returns

- A named tuple containing equal-length `columns`, `row_order`, and the
  quantity and unit metadata in `observation_columns`.

# Errors

- Throws `DimensionMismatch` when the generic observation columns do not have
  equal lengths.
- Throws `ArgumentError` when two requests would create the same scientific
  column name.
"""
function publication_table(source, requests::Tuple, observations::Tuple, options::NamedTuple)
    names = _publication_names(observations)
    columns = map(payload -> _publication_column(payload.values), observations)
    isempty(columns) || all(length(column) == length(first(columns)) for column in columns) ||
        throw(DimensionMismatch("published observation columns must have equal lengths"))
    return (
        columns = NamedTuple{names}(columns),
        row_order = names,
        observation_columns = _publication_contract(names, observations),
    )
end

function observables(
        source,
        requests::Tuple;
        units::Tuple = (),
        length_unit::Symbol = :kilo,
        frequency_unit::Symbol = :base,
        quantity_units = nothing,
        clip::Bool = true,
        atol = nothing,
        frequencies = nothing
)
    identities = validate_observables(source, requests, units)
    atol isa Real && length(unique(request_quantity.(requests))) > 1 && throw(ArgumentError(
        "a scalar atol requires one observable quantity; use keyed native-unit tolerances",
    ))
    isempty(units) || quantity_units === nothing || throw(ArgumentError(
        "use either aligned units or quantity_units, not both",
    ))
    overrides = if isempty(units)
        map(requests) do request
            scientific_quantity = _quantity(request)
            scientific_quantity isa Units.Quantity{:frequency} ?
                Units.units(frequency_unit, :hertz) :
                display_unit(
                    scientific_quantity,
                    basis(source),
                    _unit_override(quantity_units, request);
                    length_prefix = length_unit
                )
        end
    else
        units
    end
    publications = map(requests, identities, overrides) do request, identity, override
        _publish_observable(source, request, identity, override, clip, atol, frequencies)
    end
    payloads = map(publication -> publication.observation, publications)
    table = publication_table(
        source,
        requests,
        payloads,
        (; length_unit, frequency_unit, quantity_units, clip, atol)
    )
    contracts = map(keys(table.observation_columns), values(table.observation_columns)) do name, contract
        selected = findall(eachindex(requests)) do index
            haskey(contract, :requests) ? requests[index] in contract.requests :
                Symbol(Units.symbol(payloads[index].quantity)) == name
        end
        isempty(selected) && return contract
        resolutions = Tuple(publications[index].resolution for index in selected)
        merge(contract, (requests=Tuple(requests[index] for index in selected),
            observation_indices=Tuple(selected),
            resolution=length(resolutions) == 1 ? only(resolutions) : resolutions))
    end
    metadata = (
        basis = basis(source),
        row_order = table.row_order,
        observation_columns = NamedTuple{keys(table.observation_columns)}(contracts),
    )
    return ObservationPublication(payloads, table.columns, metadata)
end

basis(publication::ObservationPublication) = publication.metadata.basis

function observation_request(publication::ObservationPublication,request)
    identity=request_identity(request)
    any(contract -> any(stored -> request_identity(stored)==identity,
        get(contract,:requests,())),values(publication.metadata.observation_columns)) ||
        throw(ArgumentError("the requested product was not retained in this publication"))
    return (;identity,quantity=request_quantity(request),indices=request_indices(request))
end

"""
$(TYPEDSIGNATURES)

Read a retained publication product in native units. The publication retains
no source object and performs no reconstruction of absent statistical products.
An ambiguous retained selection must be selected explicitly before reuse.
"""
function observe(publication::ObservationPublication, selectors...)
    request=length(selectors)==1 ? only(selectors) : selectors
    resolved=observation_request(publication,request)
    matching=[(contract,index) for contract in Base.values(publication.metadata.observation_columns)
        for (stored,index) in zip(get(contract,:requests,()),get(contract,:observation_indices,()))
        if request_identity(stored)==resolved.identity]
    length(matching)==1 || throw(ArgumentError("publication product selection is ambiguous"))
    _,index=only(matching)
    payload=publication[index]
    factor=scale_factor(payload.unit,native_unit(payload.quantity,basis(publication)))
    values=detach(payload.values,factor)
    indices=resolved.indices
    if resolved.identity isa Tuple && length(resolved.identity)==3
        isempty(indices) || first(indices)==1 || throw(ArgumentError("this publication retains one selected point"))
        isempty(indices) || (indices=Base.tail(indices))
    end
    return isempty(indices) ? values : values[indices...]
end
