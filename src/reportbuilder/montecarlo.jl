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
    physical_contract=map(record -> (;record.quantity,record.unit),contract)
    all(publication -> map(record -> (;record.quantity,record.unit),publication.metadata.observation_columns) == physical_contract,
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

"""Publish retained UQ statistics and sampling evidence without inspecting result fields."""
function tabulate(definition::BenchmarkTableDefinition, operands::NamedTuple{(:reference,:candidate)};
        labels=(reference="Reference",candidate="Candidate"))
    statistics=DataFrame()
    statistic_frames=DataFrame[]
    sampling=DataFrame()
    mean_sampling_precision=DataFrame()
    for (role,operand) in pairs(operands)
        result=operand.result
        if result isa ObservationPublication
            # Historical mean/std rows are already a detached owned table;
            # no trial count, empirical quantile or dependence is invented.
            selected_quantities=unique(request_quantity.(definition.settings.requests))
            coordinates=result.metadata.row_order
            columns=Tuple(name for (name,contract) in pairs(result.metadata.observation_columns)
                if contract.quantity in selected_quantities)
            selected_names=unique((coordinates...,columns...))
            frame=DataFrame(result)[!,collect(selected_names)]
            frame[!,:role]=fill(role,size(frame,1))
            frame[!,:estimator]=fill(:retained_statistic,size(frame,1))
            push!(statistic_frames,frame)
            continue
        end
        result isa AbstractUncertaintyResult || continue
        quantities=unique(identity[2] for identity in map(request_identity,definition.settings.requests)
            if identity isa Tuple && first(identity)===UQ.statistics)
        declared=observables(typeof(result))
        products=Tuple(Iterators.flatten((UQ.statistics,quantity) in declared ?
            ((UQ.statistics,quantity),) :
            Tuple(identity for identity in declared if identity isa Tuple && length(identity)==3 &&
                first(identity)===UQ.statistics && identity[2]===quantity) for quantity in quantities))
        for point in eachindex(result)
            if !isempty(products)
                requests=Tuple((identity...,point) for identity in products)
                frame=DataFrame(observables(result,requests;length_unit=:base,clip=false))
                frame[!,:role]=fill(role,size(frame,1))
                frame[!,:estimator]=fill(result isa UQ.MonteCarloResult ? :empirical : :first_order,size(frame,1))
                push!(statistic_frames,frame)
            end
            result isa UQ.MonteCarloResult || continue
            precision=UQ.confidence(result,point)
            method=getproperty(labels,role)
            push!(sampling,(;role,method,point,trials=precision.trials,
                spread_estimated=precision.spread_estimated,
                marginal_count=precision.marginal_count,confidence=precision.confidence,
                target_cdf=precision.target_cdf,cdf_bound=precision.cdf_bound,
                target_supported=precision.target_supported,scope=precision.scope,
                distribution=precision.distribution,conditioning=precision.conditioning,
                samples_retained=precision.samples_retained,histograms_retained=precision.histograms_retained,
                std_sampling_precision=missing);cols=:union)
            ports=operand.metadata.port_order
            frequency=Engine.frequencies(result[point])
            for quantity in quantities
                # Standard error of the sampled mean is sigma_hat/sqrt(n),
                # not a confidence bound for sigma_hat itself.
                spread=observe(result,UQ.statistics,quantity,Statistics.std,point)
                unit=Units.native_unit(Units.quantity(quantity),basis(result[point]))
                for k in eachindex(frequency),i in eachindex(ports),j in eachindex(ports)
                    value=precision.trials>1 ? spread[i,j,k]/sqrt(precision.trials) : missing
                    push!(mean_sampling_precision,(;role,method,point,
                        quantity=Symbol(Units.symbol(Units.quantity(quantity))),row=i,column=j,
                        response=ports[i],excitation=ports[j],frequency_Hz=frequency[k],
                        mean_standard_error=value,unit=Units.label(unit));cols=:union)
                end
            end
        end
    end
    # Appending publications can discard their metadata, and historical tables
    # may use different prefixes. Normalize detached columns before combining
    # them; retain units, not point-specific request/resolution metadata.
    columns=(;)
    for frame in statistic_frames
        frame_basis=metadata(frame,"basis")
        for (name,contract) in pairs(observation_columns(frame))
            name in propertynames(frame) || continue
            unit=Units.native_unit(contract.quantity,frame_basis)
            column=(;contract.quantity,unit)
            haskey(columns,name) && columns[name]!=column && throw(ArgumentError(
                "retained statistics have inconsistent physical columns"))
            frame[!,name]=detach(frame[!,name],Units.scale_factor(contract.unit,unit))
            columns=merge(columns,NamedTuple{(name,)}((column,)))
        end
        append!(statistics,frame;cols=:union)
    end
    if !isempty(statistic_frames)
        metadata!(statistics,"observation_columns",columns;style=:note)
        metadata!(statistics,"basis",metadata(first(statistic_frames),"basis");style=:note)
        metadata!(statistics,"row_order",metadata(first(statistic_frames),"row_order");style=:note)
    end
    return (;statistics,sampling,mean_sampling_precision)
end
