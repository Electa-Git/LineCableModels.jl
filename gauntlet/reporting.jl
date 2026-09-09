import DataFrames: DataFrame
import LineCableModels.ReportBuilder: BenchmarkTableDefinition, select, tabulate
import LineCableModels.PlotBuilder: plot

"""
    select(definition::BenchmarkTableDefinition, benchmark)

Select the saved analyses and explicit operands returned by `read_benchmark` for
REPL reporting. Every analysis must name the loaded operands. No comparison is
recalculated and no formulation is inferred from a calculation label.
"""
function select(::BenchmarkTableDefinition,
        benchmark::NamedTuple{(:id, :reference, :candidate, :analyses)})
    isempty(benchmark.analyses) && throw(ArgumentError("benchmark has no saved analyses"))
    for (role, operand) in pairs((reference=benchmark.reference, candidate=benchmark.candidate))
        validate(read_calculation, operand.result, operand.metadata)
    end
    for record in benchmark.analyses, role in (:reference, :candidate)
        recorded = getproperty(record["calculations"], role)
        loaded = getproperty(benchmark, role).metadata
        recorded.sha256 == loaded.sha256 || throw(ArgumentError(
            "analysis $role differs from the loaded operand"))
        recorded.port_order == loaded.port_order &&
        recorded.frequencies == loaded.frequencies &&
        recorded.basis == loaded.basis && recorded.domain == loaded.domain ||
            throw(ArgumentError("analysis $role coordinates differ from the loaded operand"))
    end
    return benchmark
end

"""
    tabulate(definition::BenchmarkTableDefinition, benchmark, selected)

Expose saved benchmark data as three DataFrames: `calculations` contains operand
identities and retained configuration; `comparisons` contains complete RMS
matrices; `terms` contains one row per matrix entry. Absolute errors use the native
quantity units and recorded basis; relative errors are percentages. Per-entry
statuses and reasons remain available, including when only relative RMS is
unavailable. No matrix triangle or recorded quantity is suppressed.

Use `BenchmarkTableDefinition(false)` to preserve even the smallest stored
values. With clipping enabled, only detached numerical display values are clipped
according to the existing ReportBuilder settings. Stored data remain unchanged.
"""
function tabulate(definition::BenchmarkTableDefinition,
        benchmark::NamedTuple{(:id, :reference, :candidate, :analyses)}, selected)
    calculations = DataFrame([
        (role, backend = operand.metadata.backend,
            selection = operand.metadata.selection, formulation = operand.metadata.formulation,
            sha256 = operand.metadata.sha256, input_sha256 = operand.metadata.input_sha256,
            implementation = operand.metadata.implementation,
            path = operand.metadata.path, port_order = copy(operand.metadata.port_order),
            frequencies = copy(operand.metadata.frequencies), basis = operand.metadata.basis,
            domain = operand.metadata.domain, timing = operand.metadata.timing,
            axes = operand.metadata.axes)
        for (role, operand) in pairs((reference = selected.reference, candidate = selected.candidate))])
    comparisons = NamedTuple[]
    terms = NamedTuple[]
    for (analysis, record) in enumerate(selected.analyses), error in record["reference_comparison"]
        ports = record["port_order"]
        dimensions = (length(ports), length(ports))
        size(error.absolute) == size(error.relative) == dimensions ||
            throw(DimensionMismatch("saved RMS matrices differ from the terminal order"))
        unit = Units.native_unit(
            getproperty(LineCableModels, error.quantity), record["basis"])
        absolute = Grammar.detach(error.absolute, 1, definition.clip)
        relative = Grammar.detach(error.relative, 100, definition.clip)
        details = error.details
        statuses = get(details, :status, fill(:not_recorded, dimensions))
        local_reasons = get(details, :normalization_reason, fill(nothing, dimensions))
        size(statuses) == size(local_reasons) == dimensions ||
            throw(DimensionMismatch("saved RMS classifications differ from the terminal order"))
        reasons = map(local_reasons) do reason
            reason === nothing ? get(details, :reason, nothing) : reason
        end
        identity = (analysis, benchmark = record["benchmark_id"], collection = record["collection"],
            quantity = error.quantity, statistic = error.statistic,
            reference_point = error.reference_index, candidate_point = error.candidate_index,
            band = details.band, normalization = details.normalization,
            absolute_unit = Units.label(unit),
            samples = details.sample_count, requested_bounds_Hz = get(details, :requested_bounds, missing),
            actual_bounds_Hz = details.actual_bounds)
        push!(comparisons, merge(identity,
            (absolute_rms = absolute, relative_rms_percent = relative,
                status = copy(statuses), reason = reasons, port_order = copy(ports),
                sample_indices = copy(get(details, :indices, Int[])),
                tolerance = get(details, :atol, missing))))
        for row in eachindex(ports), column in eachindex(ports)
            push!(terms, merge(identity,
                (row, column, response = ports[row], excitation = ports[column],
                    absolute_rms = absolute[row, column], relative_rms_percent = relative[row, column],
                    status = statuses[row, column], reason = reasons[row, column])))
        end
    end
    return (calculations, comparisons = DataFrame(comparisons), terms = DataFrame(terms))
end

"""
    plot(benchmark, requests; pair=nothing, kwargs...)

Overlay one explicitly compared pair of saved line-parameter results using the
existing Makie matrix-cell recipe. `benchmark` is the loaded result bundle from
`read_benchmark`; `requests` uses the ordinary PlotBuilder observation grammar.
With multiple saved point pairs, supply `pair=(reference_index, candidate_index)`.
The pair must occur in a retained analysis, including for zipped result spaces.

Plot keywords pass to the existing multi-result recipe. Defaults preserve small
values (`clip=false`), use a logarithmic frequency axis and label the declared
reference/candidate IDs and point indices. Backends remain explicit in the configuration table. Loading GLMakie, CairoMakie
or WGLMakie enables plotting; tables require none of them. This method performs
no solve or RMS calculation.
"""
function plot(benchmark::NamedTuple{(:id, :reference, :candidate, :analyses)}, requests;
        pair = nothing, xscale = :log10, clip::Bool = false, series_labels = nothing,
        legend_position = :bottom,
        panel_titles = facet -> "$(Units.symbol(facet.quantity))[$(facet.row),$(facet.column)]",
        kwargs...)
    selected = select(BenchmarkTableDefinition(false), benchmark)
    declared = unique([(error.reference_index, error.candidate_index)
                       for record in selected.analyses for error in record["reference_comparison"]])
    if pair === nothing
        length(declared) == 1 || throw(ArgumentError(
            "benchmark has $(length(declared)) point pairs; choose pair=(reference_index, candidate_index)"))
        pair = only(declared)
    end
    pair isa Tuple{Integer, Integer} && all(index -> !(index isa Bool), pair) && pair in declared ||
        throw(ArgumentError("pair must identify an explicitly saved reference/candidate comparison"))
    pair = Int.(pair)
    sources = map((selected.reference, selected.candidate), pair) do operand, index
        select(operand.result, index)
    end
    labels = if series_labels === nothing
        map((:reference, :candidate), (selected.reference, selected.candidate), pair) do role, operand, index
            id = getproperty(first(selected.analyses)["calculations"], role).id
            "$role: $id · point $index"
        end
    else
        series_labels
    end
    return plot(sources..., requests; series_labels = labels, xscale, clip,
        legend_position, panel_titles, kwargs...)
end

function select(result::AbstractCoreResult, index::Integer)
    index == 1 || throw(ArgumentError("a scalar operand has only point 1"))
    return result
end

function select(result::AbstractParametricResult, index::Integer)
    return result[index]
end

function select(result::MomentResult, index::Integer)
    throw(ArgumentError("matrix-curve overlays require line-parameter operands; moment errors remain available through report"))
end
