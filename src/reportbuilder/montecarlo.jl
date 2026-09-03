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

function _monte_carlo_requests(
        ::UQ.MonteCarloResult{<:Engine.CableConstants},
        point::Int
)
    return (
        (UQ.statistics, R, point),
        (UQ.statistics, L, point),
        (UQ.statistics, C, point),
        (UQ.statistics, G, point)
    )
end

function _monte_carlo_requests(
        ::UQ.MonteCarloResult{<:Engine.LineParameters},
        point::Int
)
    return (
        (UQ.statistics, R, point),
        (UQ.statistics, L, point),
        (UQ.statistics, G, point),
        (UQ.statistics, C, point)
    )
end

function select(
        definition::MonteCarloTableDefinition,
        source::UQ.MonteCarloResult
)
    return [observables(
                source,
                _monte_carlo_requests(source, point);
                length_unit = definition.length_unit,
                quantity_units = definition.quantity_units,
                clip = definition.clip
            )
            for point in 1:length(source)]
end

function _monte_carlo_metadata!(table::DataFrame, source, publications)
    metadata!(table, "basis", basis(source), style = :note)
    metadata!(
        table,
        "monte_carlo",
        (
            confidence = UQ.confidence(source),
            cdf_tolerance = UQ.cdf_tolerance(source),
            distribution = UQ.sampling_distribution(source),
            root_seed = UQ.root_seed(source),
            point_seeds = [UQ.point_seed(source, point) for point in 1:length(source)],
            trial_counts = [UQ.trial_count(source, point) for point in 1:length(source)]
        );
        style = :note
    )
    first_publication = first(publications)
    metadata!(table, "row_order", first_publication.metadata.row_order, style = :note)
    metadata!(
        table,
        "observation_columns",
        first_publication.metadata.observation_columns,
        style = :note
    )
    return table
end

function tabulate(::MonteCarloTableDefinition, source, publications)
    isempty(publications) && throw(ArgumentError(
        "Monte Carlo tables require at least one Gridspace point",
    ))
    contract = first(publications).metadata.observation_columns
    row_order = first(publications).metadata.row_order
    all(publication -> publication.metadata.observation_columns == contract,
        publications) || throw(DimensionMismatch(
        "Monte Carlo points publish different scientific columns",
    ))
    all(publication -> publication.metadata.row_order == row_order,
        publications) || throw(DimensionMismatch(
        "Monte Carlo points publish different row coordinates",
    ))
    table = reduce(vcat, (DataFrame(publication) for publication in publications))
    return _monte_carlo_metadata!(table, source, publications)
end
