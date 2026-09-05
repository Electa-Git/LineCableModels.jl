const _SAMPLE_STATISTICS = (:mean, :std, :min, :q05, :median, :q95, :max)

function _uq_publication_names(observations::Tuple)
    names = map(payload -> Symbol(Units.symbol(payload.quantity)), observations)
    length(unique(names)) == length(names) || throw(ArgumentError(
        "Monte Carlo table quantities must be distinct",
    ))
    return names
end

function _uq_publication_contract(names::Tuple, observations::Tuple)
    records = map(observations) do payload
        (; quantity = payload.quantity, unit = payload.unit)
    end
    return NamedTuple{names}(records)
end

function _statistics_request(request)
    identity = request_identity(request)
    return identity isa Tuple && first(identity) === statistics
end

function _polynomial_chaos_publication_name(request, payload)
    identity = request_identity(request)
    if identity isa Tuple && first(identity) === statistics
        selector_count = length(identity)
        remainder = request isa Tuple ? request[(selector_count + 1):end] : ()
        prefix = !isempty(remainder) && first(remainder) isa _PolynomialChaosStatisticSelector ?
                 nameof(first(remainder)) : :statistics
        return Symbol(prefix, :_, Units.symbol(payload.quantity))
    end
    return Symbol(Units.symbol(payload.quantity))
end

_polynomial_chaos_publication_column(value::Number) = [value]
_polynomial_chaos_publication_column(value::AbstractArray) = collect(vec(value))
_polynomial_chaos_publication_column(value) = [value]

function publication_table(
        source::PolynomialChaosResult,
        requests::Tuple,
        observations::Tuple,
        options::NamedTuple
)
    names = map(_polynomial_chaos_publication_name, requests, observations)
    length(unique(names)) == length(names) || throw(ArgumentError(
        "polynomial-chaos publication columns must be distinct",
    ))
    columns = map(
        payload -> _polynomial_chaos_publication_column(payload.values),
        observations
    )
    isempty(columns) || all(length(column) == length(first(columns)) for column in columns) ||
        throw(DimensionMismatch(
            "polynomial-chaos publication columns must have equal lengths",
        ))
    return (
        columns = NamedTuple{names}(columns),
        row_order = names,
        observation_columns = _uq_publication_contract(names, observations),
    )
end

function _uq_product_publication(requests::Tuple, observations::Tuple)
    names = map(requests, observations) do request, payload
        identity = request_identity(request)
        product = identity isa Tuple ? first(identity) : identity
        Symbol(nameof(product), :_, Units.symbol(payload.quantity))
    end
    length(unique(names)) == length(names) || throw(ArgumentError(
        "Monte Carlo product columns must be distinct",
    ))
    columns = map(payload -> [payload.values], observations)
    contract = _uq_publication_contract(names, observations)
    return (
        columns = NamedTuple{names}(columns),
        row_order = names,
        observation_columns = contract,
    )
end

function _statistics_point(request)
    indices = request_indices(request)
    isempty(indices) && throw(ArgumentError(
        "Monte Carlo table observations require an explicit point index",
    ))
    point = first(indices)
    point isa Integer || throw(ArgumentError(
        "Monte Carlo table point indices must be integers",
    ))
    return Int(point), Base.tail(indices)
end

function _summary_value(summary::SampleSummary, statistic::Symbol)
    return getproperty(summary, statistic)
end

function publication_table(
        source::MonteCarloResult{<:Engine.CableConstants},
        requests::Tuple,
        observations::Tuple,
        options::NamedTuple
)
    all(_statistics_request, requests) ||
        return _uq_product_publication(requests, observations)
    points = map(request -> first(_statistics_point(request)), requests)
    all(==(first(points)), points) || throw(DimensionMismatch(
        "one Monte Carlo publication must select one common point",
    ))
    point = first(points)
    selections = map(requests) do request
        indices = last(_statistics_point(request))
        if isempty(indices)
            collect(eachindex(source.values[point].cores))
        elseif length(indices) == 1
            observation_indices(only(indices), length(source.values[point].cores))
        else
            throw(ArgumentError(
                "cable-constant summaries accept at most one assembly index",
            ))
        end
    end
    all(==(first(selections)), selections) || throw(DimensionMismatch(
        "cable-constant summary requests must select common assemblies",
    ))
    assemblies = first(selections)
    names = _uq_publication_names(observations)
    values = map(observations) do payload
        summaries = payload.values isa AbstractVector ? payload.values : [payload.values]
        length(summaries) == length(assemblies) || throw(DimensionMismatch(
            "cable-constant summary values must align with selected assemblies",
        ))
        [
            _summary_value(summaries[position], statistic)
            for position in eachindex(assemblies) for statistic in _SAMPLE_STATISTICS
        ]
    end
    row_count = length(assemblies) * length(_SAMPLE_STATISTICS)
    columns = merge(
        (
            point = fill(point, row_count),
            core = repeat(source.values[point].cores[assemblies];
                inner = length(_SAMPLE_STATISTICS)),
            statistic = repeat(collect(_SAMPLE_STATISTICS); outer = length(assemblies)),
        ),
        NamedTuple{names}(values),
        (
            trials = fill(trial_count(source, point), row_count),
            point_seed = fill(point_seed(source, point), row_count),
        )
    )
    return (
        columns,
        row_order = (:point, :core, :statistic),
        observation_columns = _uq_publication_contract(names, observations),
    )
end

function _line_summary_coordinates(source, requests::Tuple)
    selections = map(requests) do request
        point, indices = _statistics_point(request)
        coordinates = if isempty(indices)
            (collect(axes(source.values[point].Z, 1)),
             collect(axes(source.values[point].Z, 2)),
             collect(axes(source.values[point].Z, 3)))
        elseif length(indices) == 3
            dimensions = size(source.values[point].Z)
            map(observation_indices, indices, dimensions)
        else
            throw(ArgumentError(
                "line-summary requests require point, row, column, and frequency indices",
            ))
        end
        return point, coordinates
    end
    all(==(first(selections)), selections) || throw(DimensionMismatch(
        "one Monte Carlo publication must select common point and matrix coordinates",
    ))
    return first(selections)
end

function publication_table(
        source::MonteCarloResult{<:Engine.LineParameters},
        requests::Tuple,
        observations::Tuple,
        options::NamedTuple
)
    all(_statistics_request, requests) ||
        return _uq_product_publication(requests, observations)
    point, coordinates = _line_summary_coordinates(source, requests)
    rows, columns, samples = coordinates
    names = _uq_publication_names(observations)
    frequency_unit = Units.units(options.frequency_unit, :hertz)
    frequency_factor = Units.scale_factor(
        Units.native_unit(Units.Quantity{:frequency}(), basis(source)),
        frequency_unit
    )
    keys = [(frequency_index, row, column, statistic)
            for frequency_index in eachindex(samples)
            for row in eachindex(rows)
            for column in eachindex(columns)
            for statistic in _SAMPLE_STATISTICS]
    values = map(observations) do payload
        [_summary_value(
             payload.values[entry[2], entry[3], entry[1]],
             entry[4]
         ) for entry in keys]
    end
    table_columns = merge(
        (
            point = fill(point, length(keys)),
            frequency = [source.values[point].f[samples[entry[1]]] * frequency_factor
                         for entry in keys],
            row = [rows[entry[2]] for entry in keys],
            column = [columns[entry[3]] for entry in keys],
            statistic = [entry[4] for entry in keys],
        ),
        NamedTuple{names}(values),
        (
            trials = fill(trial_count(source, point), length(keys)),
            point_seed = fill(point_seed(source, point), length(keys)),
        )
    )
    frequency_contract = (
        quantity = Units.Quantity{:frequency}(),
        unit = frequency_unit,
    )
    contract = merge(
        (frequency = frequency_contract,),
        _uq_publication_contract(names, observations)
    )
    return (
        columns = table_columns,
        row_order = (:point, :frequency, :row, :column, :statistic),
        observation_columns = contract,
    )
end
