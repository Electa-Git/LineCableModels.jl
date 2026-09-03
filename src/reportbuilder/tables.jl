"""
$(TYPEDEF)

Publish an [`Engine.CableConstants`](@ref) result as one R/L/C/G table with
one row per concentric assembly.

$(TYPEDFIELDS)
"""
struct CableConstantsTableDefinition <: AbstractReportDefinition
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
end
CableConstantsTableDefinition() = CableConstantsTableDefinition(true)

"""
$(TYPEDEF)

Define the observable requests and display units for one wide line-parameter
table.

$(TYPEDFIELDS)
"""
struct LineParametersTableDefinition{Q <: Tuple, U} <: AbstractReportDefinition
    "Explicit observable requests in output-column order."
    requests::Q
    "SI prefix used to display frequency."
    frequency_unit::Symbol
    "Length prefix used for per-length quantities."
    length_unit::Symbol
    "Optional display-unit overrides resolved from the requests."
    quantity_units::U
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
end

"""
$(TYPEDEF)

Publish a [`Engine.LineParametersBenchmark`](@ref) as one wide comparison
table.

$(TYPEDFIELDS)
"""
struct BenchmarkTableDefinition <: AbstractReportDefinition
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
end
BenchmarkTableDefinition() = BenchmarkTableDefinition(true)

function select(definition::CableConstantsTableDefinition, source::Engine.CableConstants)
    return observables(source, (R, L, C, G); clip = definition.clip)
end

function tabulate(
        ::CableConstantsTableDefinition,
        source,
        published::ObservationPublication
)
    return DataFrame(published)
end

function _line_definition(
        requests::Tuple,
        frequency_unit::Symbol,
        length_unit::Symbol,
        quantity_units,
        clip::Bool
)
    isempty(requests) && throw(ArgumentError(
        "line tables require at least one explicit observable request",
    ))
    all(request -> request isa Tuple, requests) || throw(ArgumentError(
        "line tables require requests constructed with @observe",
    ))
    return LineParametersTableDefinition(
        requests,
        frequency_unit,
        length_unit,
        quantity_units,
        clip
    )
end

function select(definition::LineParametersTableDefinition, source::Engine.LineParameters)
    return observables(
        source,
        definition.requests;
        frequency_unit = definition.frequency_unit,
        length_unit = definition.length_unit,
        quantity_units = definition.quantity_units,
        clip = definition.clip
    )
end

function tabulate(
        ::LineParametersTableDefinition,
        source,
        published::ObservationPublication
)
    return DataFrame(published)
end

function select(
        definition::BenchmarkTableDefinition,
        comparison::Engine.LineParametersBenchmark
)
    requests = (
        (Z, Engine.absolute_error),
        (Z, Engine.relative_error),
        (Y, Engine.absolute_error),
        (Y, Engine.relative_error)
    )
    return observables(comparison, requests; clip = definition.clip)
end

function tabulate(
        ::BenchmarkTableDefinition,
        source,
        published::ObservationPublication
)
    return DataFrame(published)
end
