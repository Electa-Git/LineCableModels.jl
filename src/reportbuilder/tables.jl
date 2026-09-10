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

Select and publish per-term comparisons of scalar or formulation-space results.
Retained products supply detailed tables and a compact summary.

$(TYPEDFIELDS)
"""
struct BenchmarkTableDefinition{S <: NamedTuple, I, O <: NamedTuple} <: AbstractReportDefinition
    "Normalized quantities, frequency bands and numerical comparison controls."
    settings::S
    "Explicitly supplied comparison controls, used when selecting retained products."
    requested::Tuple{Vararg{Symbol}}
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
    "Optional explicit illustration request; nothing produces no figure."
    illustration::I
    "Options for the explicitly requested illustration."
    plot_options::O
end

"""
$(TYPEDSIGNATURES)

Request per-term comparisons grouped by formulation and frequency band.
Quantities use the observation grammar. The default bands are the entire range,
near DC, harmonic, narrowband and wideband. `fundamental` is in Hz. No figure is
created unless `illustration` is explicitly supplied. `pairing` maps each candidate
point to an explicit reference point when both operands are result spaces.
"""
function BenchmarkTableDefinition(; clip::Bool=false, illustration=nothing, plot_options=(;), kwargs...)
    moments=Tuple(get(kwargs,:statistics,(:value,))) == (:mean,:std)
    defaults=(quantities=moments ? (R,L,C,G) : (Z,Y,R,L,G,C), statistics=(:value,),
        bands=moments ? (:all,) : (:all,:dc,:harmonic,:narrow,:wide),
        normalizations=(:reference_rms,), atol=nothing, fundamental=50.0, harmonics=50,
        unsupported=(;), pairing=nothing)
    isempty(setdiff(keys(kwargs),keys(defaults))) || throw(ArgumentError("unknown benchmark comparison controls"))
    supplied=merge(defaults,(;kwargs...))
    settings=merge(supplied,(
        quantities=Tuple(q isa Symbol ? q : Symbol(nameof(q)) for q in supplied.quantities),
        statistics=Tuple(supplied.statistics), bands=Tuple(supplied.bands),
        normalizations=Tuple(supplied.normalizations)))
    return validate(BenchmarkTableDefinition(settings,Tuple(keys(kwargs)),clip,illustration,plot_options))
end

function BenchmarkTableDefinition(clip::Bool; kwargs...)
    return BenchmarkTableDefinition(; clip, kwargs...)
end

"""Validate comparison requests before numerical execution."""
function validate(definition::BenchmarkTableDefinition)
    settings=definition.settings
    !isempty(settings.quantities) && allunique(settings.quantities) &&
        all(q -> q in (:Z, :Y, :R, :L, :G, :C), settings.quantities) ||
        throw(ArgumentError("benchmark quantities must select distinct Z, Y, R, L, G, or C"))
    if settings.pairing !== nothing
        settings.pairing isa Union{Tuple,AbstractVector} && !isempty(settings.pairing) &&
            all(pair -> pair isa Tuple{Integer,Integer} && all(index -> !(index isa Bool) && index>0,pair),settings.pairing) &&
            sort(last.(collect(settings.pairing))) == collect(1:length(settings.pairing)) ||
            throw(ArgumentError("pairing must list positive reference/candidate indices with each candidate exactly once"))
    end
    settings.statistics in ((:value,), (:mean, :std)) ||
        throw(ArgumentError("benchmark statistics must be (:value,) or (:mean, :std)"))
    !isempty(settings.bands) && allunique(settings.bands) ||
        throw(ArgumentError("benchmark needs distinct frequency bands"))
    !isempty(settings.normalizations) && allunique(settings.normalizations) ||
        throw(ArgumentError("benchmark needs distinct RMS normalizations"))
    for band in settings.bands, normalization in settings.normalizations
        validate(Engine.compare; band, normalization, atol=settings.atol,
            fundamental=settings.fundamental, harmonics=settings.harmonics,
            unsupported=settings.unsupported)
    end
    if settings.statistics == (:mean, :std)
        settings.quantities == (:R, :L, :C, :G) && settings.bands == (:all,) &&
            settings.normalizations == (:reference_rms,) && settings.atol === nothing &&
            isempty(settings.unsupported) || throw(ArgumentError(
                "moment comparisons support full-band R/L/C/G means and standard deviations with reference-RMS normalization"))
    end
    return definition
end

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
