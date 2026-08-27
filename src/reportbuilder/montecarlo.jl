"""
$(TYPEDEF)

Define display units for one wide Monte Carlo summary table.

$(TYPEDFIELDS)
"""
struct MonteCarloTableDefinition{U} <: AbstractReportDefinition
    "Length prefix used for per-length quantities."
    length_unit::Symbol
    "Optional display-unit overrides resolved from the published quantities."
    quantity_units::U
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
end
function MonteCarloTableDefinition(length_unit::Symbol, quantity_units)
    MonteCarloTableDefinition(length_unit, quantity_units, true)
end

illustrate(::MonteCarloTableDefinition, source, published, table) = nothing
encode(::MonteCarloTableDefinition, source, published, table, ::Nothing) = nothing
write(::MonteCarloTableDefinition, source, published, table, ::Nothing, ::Nothing) = nothing

entitle(::MonteCarloTableDefinition, source::UQ.MonteCarloResult) = source

function _monte_carlo_context(source, point::Int)
    return (
        point,
        trials = UQ.trial_count(source, point),
        seed = UQ.point_seed(source, point)
    )
end

function select(
        definition::MonteCarloTableDefinition,
        source::UQ.MonteCarloResult{<:DataModel.CableConstants}
)
    return map(1:length(source)) do point
        requests = (
            (UQ.statistics, R, point),
            (UQ.statistics, L, point),
            (UQ.statistics, C, point)
        )
        targets = unit_targets(
            requests,
            basis(source);
            length_prefix = definition.length_unit,
            overrides = definition.quantity_units
        )
        observations = observables(
            source,
            requests;
            units = targets,
            clip = definition.clip
        )
        return (; frequency = nothing, observations,
            context = _monte_carlo_context(source, point))
    end
end

function select(
        definition::MonteCarloTableDefinition,
        source::UQ.MonteCarloResult{<:Engine.LineParameters}
)
    return map(1:length(source)) do point
        requests = (
            (UQ.statistics, R, point),
            (UQ.statistics, L, point),
            (UQ.statistics, G, point),
            (UQ.statistics, C, point)
        )
        targets = unit_targets(
            requests,
            basis(source);
            length_prefix = definition.length_unit,
            overrides = definition.quantity_units
        )
        all_requests = ((frequencies, point, Colon()), requests...)
        all_targets = (Units.units(:base, :hertz), targets...)
        published = observables(
            source,
            all_requests;
            units = all_targets,
            clip = definition.clip
        )
        return (
            frequency = first(published),
            observations = Base.tail(published),
            context = _monte_carlo_context(source, point)
        )
    end
end

const _SUMMARY_STATISTICS = (:mean, :std, :min, :q05, :median, :q95, :max)

function _monte_carlo_metadata!(table::DataFrame, source, selected)
    metadata!(table, "basis", basis(source), style = :note)
    metadata!(
        table,
        "monte_carlo",
        (
            confidence = UQ.confidence(source),
            cdf_tolerance = UQ.cdf_tolerance(source),
            distribution = UQ.sampling_distribution(source),
            root_seed = UQ.root_seed(source),
            point_seeds = [item.context.seed for item in selected],
            trial_counts = [item.context.trials for item in selected]
        );
        style = :note
    )
    return table
end

function _cable_summary_table(source, selected)
    first_observations = first(selected).observations
    names = _observation_names(first_observations)
    all(item -> _observation_names(item.observations) == names, selected) || throw(
        DimensionMismatch("Monte Carlo points publish different cable quantities"),
    )
    point = [item.context.point for item in selected for _ in _SUMMARY_STATISTICS]
    statistic = [statistic for _ in selected for statistic in _SUMMARY_STATISTICS]
    values = ntuple(length(names)) do quantity_index
        [getproperty(item.observations[quantity_index].values, statistic)
         for item in selected for statistic in _SUMMARY_STATISTICS]
    end
    trials = [item.context.trials for item in selected for _ in _SUMMARY_STATISTICS]
    point_seed = [item.context.seed for item in selected for _ in _SUMMARY_STATISTICS]
    table = DataFrame(merge(
        (; point, statistic),
        NamedTuple{names}(values),
        (; trials, point_seed)
    ))
    metadata!(table, "row_order", (:point, :statistic), style = :note)
    _observation_columns!(table, names, first_observations)
    return _monte_carlo_metadata!(table, source, selected)
end

function _validate_line_summary(selected)
    names = _observation_names(first(selected).observations)
    dimensions = size(first(first(selected).observations).values)
    for item in selected
        _observation_names(item.observations) == names || throw(DimensionMismatch(
            "Monte Carlo points publish different line quantities",
        ))
        all(payload -> size(payload.values) == dimensions, item.observations) || throw(
            DimensionMismatch("Monte Carlo line summaries do not share one matrix shape"),
        )
        length(item.frequency.values) == dimensions[3] || throw(DimensionMismatch(
            "frequency count does not match Monte Carlo line summaries",
        ))
    end
    return names, dimensions
end

function _line_summary_table(source, selected)
    names, dimensions = _validate_line_summary(selected)
    row_keys = [(item, frequency_index, row, column, statistic)
                for item in selected
                for frequency_index in 1:dimensions[3]
                for row in 1:dimensions[1]
                for column in 1:dimensions[2]
                for statistic in _SUMMARY_STATISTICS]
    point = [entry[1].context.point for entry in row_keys]
    frequency = [entry[1].frequency.values[entry[2]] for entry in row_keys]
    row = [entry[3] for entry in row_keys]
    column = [entry[4] for entry in row_keys]
    statistic = [entry[5] for entry in row_keys]
    values = ntuple(length(names)) do quantity_index
        [getproperty(
             entry[1].observations[quantity_index].values[entry[3], entry[4], entry[2]],
             entry[5]
         ) for entry in row_keys]
    end
    trials = [entry[1].context.trials for entry in row_keys]
    point_seed = [entry[1].context.seed for entry in row_keys]
    table = DataFrame(merge(
        (; point, frequency, row, column, statistic),
        NamedTuple{names}(values),
        (; trials, point_seed)
    ))
    metadata!(
        table,
        "row_order",
        (:point, :frequency, :row, :column, :statistic),
        style = :note
    )
    frequency_observation = first(selected).frequency
    observation_names = (:frequency, names...)
    observations = (frequency_observation, first(selected).observations...)
    _observation_columns!(table, observation_names, observations)
    return _monte_carlo_metadata!(table, source, selected)
end

function tabulate(::MonteCarloTableDefinition, source, selected)
    isempty(selected) && throw(ArgumentError(
        "Monte Carlo tables require at least one Gridspace point",
    ))
    first(selected).frequency === nothing && return _cable_summary_table(source, selected)
    return _line_summary_table(source, selected)
end

"""
$(TYPEDSIGNATURES)

Return one wide table of Monte Carlo summaries for every Gridspace point.
"""
function DataFrame(
        source::UQ.MonteCarloResult;
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        clip::Bool = true
)
    return report(
        MonteCarloTableDefinition(length_unit, quantity_units, clip),
        source
    ).table::DataFrame
end
