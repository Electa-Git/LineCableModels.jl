struct MonteCarloTable{U} <: AbstractReportDefinition
    length_unit::Symbol
    quantity_units::U
end

illustrate(::MonteCarloTable, source, published, table) = nothing
encode(::MonteCarloTable, source, published, table, ::Nothing) = nothing
write(::MonteCarloTable, source, published, table, ::Nothing, ::Nothing) = nothing

entitle(::MonteCarloTable, source::UQ.MonteCarloResult) = source

function _monte_carlo_requests(::DataModel.CableConstants)
    return (R = R, L = L, C = C)
end

function _monte_carlo_requests(::Engine.LineParameters)
    return (R = R, L = L, C = C, G = G)
end

function _monte_carlo_requests(representation)
    throw(ArgumentError(
        "Monte Carlo table presentation is not defined for $(typeof(representation))",
    ))
end

function select(definition::MonteCarloTable, source::UQ.MonteCarloResult)
    length(source) == 1 || throw(ArgumentError(
        "DataFrame requires one Monte Carlo point; select one result explicitly",
    ))
    representation = only(source)
    product = only(UQ.statistics(source))
    requests = _monte_carlo_requests(representation)
    target_values = map(keys(requests), values(requests)) do key, request
        _request_target(
            request,
            key,
            basis(product),
            definition.length_unit,
            definition.quantity_units
        )
    end
    targets = NamedTuple{keys(requests)}(target_values)
    published = observables(product, requests; units = targets)
    context = (
        trials = UQ.trial_count(source, 1),
        confidence = UQ.confidence(source),
        cdf_tol = UQ.cdf_tolerance(source),
        distribution = UQ.sampling_distribution(source),
        root_seed = UQ.root_seed(source),
        seed = UQ.point_seed(source, 1)
    )
    return (; representation, published, context)
end

function _summary_table(published::NamedTuple, context)
    payloads = values(published)
    summaries = getproperty.(payloads, :values)
    table = DataFrame(
        quantity = [Units.symbol(payload.quantity) for payload in payloads],
        mean = getproperty.(summaries, :mean),
        std = getproperty.(summaries, :std),
        min = getproperty.(summaries, :min),
        q05 = getproperty.(summaries, :q05),
        median = getproperty.(summaries, :median),
        q95 = getproperty.(summaries, :q95),
        max = getproperty.(summaries, :max),
        n = getproperty.(summaries, :n),
        unit = [Units.label(payload.unit) for payload in payloads],
        trials = fill(context.trials, length(payloads)),
        confidence = fill(context.confidence, length(payloads)),
        cdf_tol = fill(context.cdf_tol, length(payloads))
    )
    metadata!(
        table,
        "monte_carlo",
        (
            trials = context.trials,
            confidence = context.confidence,
            cdf_tol = context.cdf_tol,
            distribution = string(context.distribution),
            seed = context.seed,
            root_seed = context.root_seed
        );
        style = :note
    )
    headings = Dict(
        key => Units.label(payload.quantity, payload.unit)
    for (key, payload) in pairs(published)
    )
    metadata!(table, "headings", headings, style = :note)
    return table
end

function tabulate(::MonteCarloTable, source, selected)
    published = selected.published
    if selected.representation isa DataModel.CableConstants
        return _summary_table(published, selected.context)
    end
    shape = size(first(values(published)).values)
    tables = Array{DataFrame, length(shape)}(undef, shape)
    for index in CartesianIndices(tables)
        payloads = map(published) do payload
            (; values = payload.values[index],
                quantity = payload.quantity, unit = payload.unit)
        end
        tables[index] = _summary_table(payloads, selected.context)
    end
    return tables
end

"""
$(TYPEDSIGNATURES)

Return published Monte Carlo summaries for one selected Gridspace point.
"""
function DataFrame(
        source::UQ.MonteCarloResult;
        length_unit::Symbol = :kilo,
        quantity_units = nothing
)
    return report(MonteCarloTable(length_unit, quantity_units), source).table
end
