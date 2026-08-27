"""
$(TYPEDEF)

Define display units for one long-form Monte Carlo summary table.

$(TYPEDFIELDS)
"""
struct MonteCarloTableDefinition{U} <: AbstractReportDefinition
    "Length prefix used for per-length quantities."
    length_unit::Symbol
    "Optional display-unit overrides aligned with the published requests."
    quantity_units::U
end

illustrate(::MonteCarloTableDefinition, source, published, table) = nothing
encode(::MonteCarloTableDefinition, source, published, table, ::Nothing) = nothing
write(::MonteCarloTableDefinition, source, published, table, ::Nothing, ::Nothing) = nothing

entitle(::MonteCarloTableDefinition, source::UQ.MonteCarloResult) = source

function select(
        definition::MonteCarloTableDefinition,
        source::UQ.MonteCarloResult{<:DataModel.CableConstants}
)
    return map(1:length(source)) do point
        requests = (
            R = (UQ.statistics, R, point),
            L = (UQ.statistics, L, point),
            C = (UQ.statistics, C, point)
        )
        targets = unit_targets(
            requests,
            basis(source);
            length_prefix = definition.length_unit,
            overrides = definition.quantity_units
        )
        published = observables(source, requests; units = targets)
        return (
            frequency = nothing,
            published,
            context = (
                point,
                trials = UQ.trial_count(source, point),
                seed = UQ.point_seed(source, point),
                confidence = UQ.confidence(source),
                cdf_tol = UQ.cdf_tolerance(source)
            )
        )
    end
end

function select(
        definition::MonteCarloTableDefinition,
        source::UQ.MonteCarloResult{<:Engine.LineParameters}
)
    return map(1:length(source)) do point
        requests = (
            R = (UQ.statistics, R, point),
            L = (UQ.statistics, L, point),
            C = (UQ.statistics, C, point),
            G = (UQ.statistics, G, point)
        )
        targets = unit_targets(
            requests,
            basis(source);
            length_prefix = definition.length_unit,
            overrides = definition.quantity_units
        )
        all_requests = merge(
            (frequency = (frequencies, point, Colon()),),
            requests
        )
        all_targets = merge(
            (frequency = Units.units(:base, :hertz),),
            targets
        )
        all_published = observables(source, all_requests; units = all_targets)
        published = NamedTuple{keys(requests)}(Tuple(
            getproperty(all_published, key) for key in keys(requests)
        ))
        return (
            frequency = all_published.frequency,
            published,
            context = (
                point,
                trials = UQ.trial_count(source, point),
                seed = UQ.point_seed(source, point),
                confidence = UQ.confidence(source),
                cdf_tol = UQ.cdf_tolerance(source)
            )
        )
    end
end

function _monte_carlo_metadata!(table::DataFrame, source, selected)
    metadata!(
        table,
        "monte_carlo",
        (
            confidence = UQ.confidence(source),
            cdf_tol = UQ.cdf_tolerance(source),
            distribution = UQ.sampling_distribution(source),
            root_seed = UQ.root_seed(source),
            point_seeds = [item.context.seed for item in selected],
            trial_counts = [item.context.trials for item in selected]
        );
        style = :note
    )
    first_published = first(selected).published
    metadata!(
        table,
        "units",
        Dict(key => Units.label(payload.unit)
        for (key, payload) in pairs(first_published)),
        style = :note
    )
    metadata!(
        table,
        "headings",
        Dict(key => Units.label(payload.quantity, payload.unit)
        for (key, payload) in pairs(first_published)),
        style = :note
    )
    return table
end

function _cable_summary_rows(selected)
    return [(
                point = item.context.point,
                quantity = key,
                mean = payload.values.mean,
                std = payload.values.std,
                min = payload.values.min,
                q05 = payload.values.q05,
                median = payload.values.median,
                q95 = payload.values.q95,
                max = payload.values.max,
                n = payload.values.n,
                unit = Units.label(payload.unit),
                trials = item.context.trials,
                confidence = item.context.confidence,
                cdf_tol = item.context.cdf_tol
            )
            for item in selected
            for (key, payload) in pairs(item.published)]
end

function _line_summary_rows(selected)
    for item in selected
        frequency_values = item.frequency.values
        first_payload = first(values(item.published))
        row_count, column_count, frequency_count = size(first_payload.values)
        frequency_count == length(frequency_values) || throw(DimensionMismatch(
            "frequency count does not match Monte Carlo line summaries",
        ))
    end
    return [(
                point = item.context.point,
                row,
                column,
                frequency = item.frequency.values[frequency_index],
                quantity = key,
                mean = payload.values[row, column, frequency_index].mean,
                std = payload.values[row, column, frequency_index].std,
                min = payload.values[row, column, frequency_index].min,
                q05 = payload.values[row, column, frequency_index].q05,
                median = payload.values[row, column, frequency_index].median,
                q95 = payload.values[row, column, frequency_index].q95,
                max = payload.values[row, column, frequency_index].max,
                n = payload.values[row, column, frequency_index].n,
                unit = Units.label(payload.unit),
                trials = item.context.trials,
                confidence = item.context.confidence,
                cdf_tol = item.context.cdf_tol
            )
            for item in selected
            for row in axes(first(values(item.published)).values, 1)
            for column in axes(first(values(item.published)).values, 2)
            for frequency_index in axes(first(values(item.published)).values, 3)
            for (key, payload) in pairs(item.published)]
end

function _monte_carlo_table(
        ::Nothing,
        source,
        selected
)
    table = DataFrame(_cable_summary_rows(selected))
    metadata!(table, "row_order", (:point, :quantity), style = :note)
    return _monte_carlo_metadata!(table, source, selected)
end

function _monte_carlo_table(
        frequency::NamedTuple,
        source,
        selected
)
    table = DataFrame(_line_summary_rows(selected))
    metadata!(table, "frequency_unit", Units.label(frequency.unit), style = :note)
    metadata!(
        table,
        "row_order",
        (:point, :row, :column, :frequency, :quantity),
        style = :note
    )
    return _monte_carlo_metadata!(table, source, selected)
end

function tabulate(::MonteCarloTableDefinition, source, selected)
    isempty(selected) && throw(ArgumentError(
        "Monte Carlo tables require at least one Gridspace point",
    ))
    return _monte_carlo_table(first(selected).frequency, source, selected)
end

"""
$(TYPEDSIGNATURES)

Return one long-form table of published Monte Carlo summaries for every
Gridspace point.
"""
function DataFrame(
        source::UQ.MonteCarloResult;
        length_unit::Symbol = :kilo,
        quantity_units = nothing
)
    return report(
        MonteCarloTableDefinition(length_unit, quantity_units),
        source
    ).table::DataFrame
end
