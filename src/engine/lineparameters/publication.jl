function _publication_contract(names::Tuple, observations::Tuple)
    records = map(observations) do payload
        (; quantity = payload.quantity, unit = payload.unit)
    end
    return NamedTuple{names}(records)
end

function _publication_quantity_names(observations::Tuple)
    names = map(payload -> Symbol(Units.symbol(payload.quantity)), observations)
    all(name -> !isempty(string(name)), names) || throw(ArgumentError(
        "every published line quantity must define a nonempty symbol",
    ))
    length(unique(names)) == length(names) || throw(ArgumentError(
        "published line quantities must be distinct",
    ))
    return names
end

function _line_publication_coordinates(source::LineParameters, requests::Tuple)
    isempty(requests) && throw(ArgumentError(
        "line observation tables require at least one result quantity",
    ))
    resolved = map(request -> observation_request(source, request), requests)
    rank = length(first(resolved).indices)
    all(request -> length(request.indices) == rank, resolved) || throw(DimensionMismatch(
        "line observations in one table must use the same coordinate rank",
    ))
    rank in (2, 3) || throw(ArgumentError(
        "line observation tables require matrix and frequency indices",
    ))
    rank_value = Val(rank)
    dimensions = size(source.Z)
    coordinates = map(resolved) do request
        _line_observation_coordinates(request, dimensions, rank_value)
    end
    all(==(first(coordinates)), coordinates) || throw(DimensionMismatch(
        "line observations in one table must select the same coordinates",
    ))
    return rank_value, first(coordinates), resolved
end

function _line_observation_coordinates(request, dimensions, ::Val{3})
    return map(observation_indices, request.indices, dimensions)
end

function _line_observation_coordinates(request, dimensions, ::Val{2})
    return (
        observation_indices(request.indices[1], dimensions[1]),
        observation_indices(request.indices[2], dimensions[3]),
    )
end

function _detached_at(values, selectors::Tuple, local_indices::Tuple)
    retained = _retained_indices(selectors, local_indices)
    return isempty(retained) ? values : getindex(values, retained...)
end

_retained_indices(::Tuple{}, ::Tuple{}) = ()

function _retained_indices(selectors::Tuple, local_indices::Tuple)
    tail = _retained_indices(Base.tail(selectors), Base.tail(local_indices))
    return _retain_index(first(selectors), first(local_indices), tail)
end

_retain_index(::Integer, local_index, tail) = tail
_retain_index(selector, local_index, tail) = (local_index, tail...)

_line_quantity_pairs(::Tuple{}, ::Tuple{}) = ()

function _line_quantity_pairs(requests::Tuple, observations::Tuple)
    tail = _line_quantity_pairs(Base.tail(requests), Base.tail(observations))
    return _line_quantity_pair(first(requests), first(observations), tail)
end

function _line_quantity_pair(
        request,
        payload::NamedTuple{
            (:values, :quantity, :unit),
            Tuple{V, Units.Quantity{:frequency}, U}
        },
        tail
) where {V, U}
    return tail
end

_line_quantity_pair(request, payload, tail) = ((request, payload), tail...)

function _frequency_only_publication(observations::Tuple)
    names = _publication_quantity_names(observations)
    columns = map(payload -> collect(vec(payload.values)), observations)
    return (
        columns = NamedTuple{names}(columns),
        row_order = names,
        observation_columns = _publication_contract(names, observations),
    )
end

function publication_table(
        source::LineParameters,
        requests::Tuple,
        observations::Tuple,
        options::NamedTuple
)
    pairs = _line_quantity_pairs(requests, observations)
    isempty(pairs) && return _frequency_only_publication(observations)
    quantity_requests = map(first, pairs)
    quantity_observations = map(last, pairs)
    names = _publication_quantity_names(quantity_observations)
    rank, coordinates, resolved = _line_publication_coordinates(source, quantity_requests)
    return _line_quantity_publication(
        source,
        quantity_observations,
        Val(names),
        coordinates,
        resolved,
        options,
        rank
    )
end

function _line_frequency_unit(source, options)
    frequency_unit = Units.units(options.frequency_unit, :hertz)
    frequency_factor = Units.scale_factor(
        Units.native_unit(Units.Quantity{:frequency}(), basis(source)),
        frequency_unit
    )
    return frequency_unit, frequency_factor
end

function _line_publication_contract(::Val{Names}, observations, frequency_unit) where {Names}
    frequency_contract = (
        quantity = Units.Quantity{:frequency}(),
        unit = frequency_unit,
    )
    return merge(
        (frequency = frequency_contract,),
        _publication_contract(Names, observations)
    )
end

function _line_quantity_publication(
        source,
        observations,
        ::Val{Names},
        coordinates,
        resolved,
        options,
        ::Val{3}
) where {Names}
    rows, columns, samples = coordinates
    frequency_unit, frequency_factor = _line_frequency_unit(source, options)
    frequency = [source.f[k] * frequency_factor
                 for k in samples for _ in rows for _ in columns]
    row_column = [row for _ in samples for row in rows for _ in columns]
    column_column = [column for _ in samples for _ in rows for column in columns]
    quantity_values = map(observations, resolved) do payload, request
        [_detached_at(
             payload.values,
             request.indices,
             (local_row, local_column, local_frequency)
         )
         for local_frequency in eachindex(samples)
         for local_row in eachindex(rows)
         for local_column in eachindex(columns)]
    end
    table_columns = merge(
        (; frequency, row = row_column, column = column_column),
        NamedTuple{Names}(quantity_values)
    )
    return (
        columns = table_columns,
        row_order = (:frequency, :row, :column),
        observation_columns = _line_publication_contract(
            Val(Names), observations, frequency_unit
        ),
    )
end


function _line_quantity_publication(
        source,
        observations,
        ::Val{Names},
        coordinates,
        resolved,
        options,
        ::Val{2}
) where {Names}
    rows, samples = coordinates
    frequency_unit, frequency_factor = _line_frequency_unit(source, options)
    frequency = [source.f[k] * frequency_factor for k in samples for _ in rows]
    row_column = [row for _ in samples for row in rows]
    quantity_values = map(observations, resolved) do payload, request
        [_detached_at(
             payload.values,
             request.indices,
             (local_row, local_frequency)
         )
         for local_frequency in eachindex(samples)
         for local_row in eachindex(rows)]
    end
    table_columns = merge(
        (; frequency, row = row_column, column = copy(row_column)),
        NamedTuple{Names}(quantity_values)
    )
    return (
        columns = table_columns,
        row_order = (:frequency, :row, :column),
        observation_columns = _line_publication_contract(
            Val(Names), observations, frequency_unit
        ),
    )
end

function publication_table(
        source::LineParametersBenchmark,
        requests::Tuple,
        observations::Tuple,
        options::NamedTuple
)
    dimensions = size(first(observations).values)
    all(payload -> size(payload.values) == dimensions, observations) || throw(
        DimensionMismatch("benchmark observations must share one matrix shape"),
    )
    rows = [row for row in 1:dimensions[1] for _ in 1:dimensions[2]]
    columns = [column for _ in 1:dimensions[1] for column in 1:dimensions[2]]
    names = _publication_quantity_names(observations)
    values = map(payload -> collect(vec(transpose(payload.values))), observations)
    return (
        columns = merge((; row = rows, column = columns), NamedTuple{names}(values)),
        row_order = (:row, :column),
        observation_columns = _publication_contract(names, observations),
    )
end
