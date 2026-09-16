const _SAMPLE_STATISTICS = (:mean, :std, :min, :q05, :median, :q95, :max)

function _statistics_request(request)
    identity = request_identity(request)
    return identity isa Tuple && first(identity) === statistics
end

function _uq_product_publication(requests::Tuple, observations::Tuple)
    names = map(requests, observations) do request, payload
        identity = request_identity(request)
        product = identity isa Tuple ? first(identity) : identity
        Symbol(nameof(product), :_, Units.symbol(payload.quantity))
    end
    allunique(names) || throw(ArgumentError("UQ product columns must be distinct"))
    columns = map(payload -> [payload.values], observations)
    records = map(requests, observations) do request, payload
        (; payload.quantity, payload.unit, requests = (request,))
    end
    return (columns = NamedTuple{names}(columns), row_order = names,
        observation_columns = NamedTuple{names}(records))
end

function _statistics_point(request)
    indices = request_indices(request)
    isempty(indices) &&
        throw(ArgumentError("UQ observations require an explicit point index"))
    point = first(indices)
    point isa Integer && !(point isa Bool) ||
        throw(ArgumentError("UQ point indices must be integers"))
    return Int(point), Base.tail(indices)
end

"""
$(TYPEDSIGNATURES)

Publish statistics on common point and physical coordinates. Statistics are
rows, not extra result-space dimensions. MC full summaries retain their seven
statistics and trial counts; LEP publishes only supported projections.
Different statistics of one quantity share a physical-unit column.
Missing cells indicate quantity/statistic combinations not requested.
"""
function publication_table(source::Union{MonteCarloResult, LinearErrorResult},
        requests::Tuple, observations::Tuple, options::NamedTuple)
    all(_statistics_request, requests) ||
        return _uq_product_publication(requests, observations)
    selections = map(requests) do request
        point, indices = _statistics_point(request)
        core = source[point]
        dimensions = core isa Engine.LineParameters ? size(observe(core, Engine.Z)) :
                     (length(core),)
        coordinates = isempty(indices) ? map(count -> collect(1:count), dimensions) :
                      length(indices) == length(dimensions) ?
                      map(observation_indices, indices, dimensions) :
                      throw(ArgumentError("statistic indices must select every physical dimension"))
        (point, coordinates)
    end
    all(==(first(selections)), selections) || throw(DimensionMismatch(
        "one UQ publication requires common point and physical coordinates"))
    point, coordinates = first(selections)
    core = source[point]
    selected_statistics = map(requests) do request
        identity = request_identity(request)
        length(identity) == 2 && return _SAMPLE_STATISTICS
        transform = last(identity)
        name = transform === minimum ? :min :
               transform === maximum ? :max :
               transform isa Base.Fix2{typeof(Statistics.quantile)} ?
               (transform.x == 0 ? :min :
                transform.x == 0.05 ? :q05 :
                transform.x == 0.5 ? :median : transform.x == 0.95 ? :q95 : :max) :
               nameof(transform)
        return (name,)
    end
    statistics_order = unique(collect(Iterators.flatten(selected_statistics)))
    names = Tuple(unique(Symbol(Units.symbol(payload.quantity))
    for payload in observations))
    physical_keys = core isa Engine.LineParameters ?
                    [(row, column, frequency) for frequency in eachindex(coordinates[3])
                     for row in eachindex(coordinates[1])
                     for column in eachindex(coordinates[2])] :
                    [(index,) for index in eachindex(only(coordinates))]
    entries = [(index, statistic) for index in physical_keys
               for statistic in statistics_order]
    columns = map(names) do name
        positions = findall(payload -> Symbol(Units.symbol(payload.quantity)) == name, observations)
        unit = observations[first(positions)].unit
        all(index -> observations[index].unit == unit, positions) || throw(ArgumentError(
            "statistics of one quantity require the same display unit"))
        offered = collect(Iterators.flatten(selected_statistics[index]
        for index in positions))
        allunique(offered) ||
            throw(ArgumentError("a quantity/statistic was requested more than once"))
        arrays = map(positions) do index
            payload = observations[index].values
            reshape(payload isa AbstractArray ? payload : [payload], length.(coordinates)...)
        end
        [begin
             position = findfirst(index -> statistic in selected_statistics[index], positions)
             if position === nothing
                 missing
             else
                 value = arrays[position][index...]
                 value isa SampleSummary ? getproperty(value, statistic) : value
             end
         end
         for (index, statistic) in entries]
    end
    records = map(names) do name
        positions = findall(payload -> Symbol(Units.symbol(payload.quantity)) == name, observations)
        payload = observations[first(positions)]
        (; payload.quantity, payload.unit,
            requests = Tuple(requests[index] for index in positions))
    end
    row_count = length(entries)
    table = if core isa Engine.LineParameters
        frequency_unit = Units.units(options.frequency_unit, :hertz)
        factor = Units.scale_factor(
            Units.native_unit(Units.Quantity{:frequency}(), basis(source)), frequency_unit)
        (point = fill(point, row_count),
            frequency = [frequencies(core)[coordinates[3][index[3]]] * factor
                         for (index, _) in entries],
            row = [coordinates[1][index[1]] for (index, _) in entries],
            column = [coordinates[2][index[2]] for (index, _) in entries],
            statistic = [statistic for (_, statistic) in entries])
    else
        (point = fill(point, row_count),
            core = [core.cores[coordinates[1][index[1]]] for (index, _) in entries],
            statistic = [statistic for (_, statistic) in entries])
    end
    row_order = keys(table)
    contract = NamedTuple{names}(records)
    if core isa Engine.LineParameters
        contract = merge(
            (frequency = (quantity = Units.Quantity{:frequency}(),
                unit = Units.units(options.frequency_unit, :hertz)),),
            contract)
    end
    table = merge(table, NamedTuple{names}(columns))
    if source isa MonteCarloResult
        table = merge(table,
            (trials = fill(trial_count(source, point), row_count),
                point_seed = fill(point_seed(source, point), row_count)))
    end
    return (columns = table, row_order, observation_columns = contract)
end
