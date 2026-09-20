import LineCableModels.PlotBuilder: plot

"""
Read the most recent retained observed benchmark, or an explicitly selected
snapshot, into the common report workflow. No comparison or source extraction
occurs during report construction.
"""
function report(definition::BenchmarkTableDefinition,
        benchmark::NamedTuple{(:id,:reference,:candidate,:analyses,:measurements)})
    isempty(benchmark.analyses) && throw(ArgumentError("benchmark has no saved observations"))
    selected=argmax(record -> record["recorded_at_utc"],benchmark.analyses)
    retained=ImportExport.deserialize_value(selected["observed_data"])
    return report(definition,retained.observed;reference=retained.reference)
end

"""Plot the same observed inputs retained by the selected benchmark report."""
function plot(benchmark::NamedTuple{(:id,:reference,:candidate,:analyses,:measurements)},
        selection=nothing;ydata=nothing,kwargs...)
    artifact=report(BenchmarkTableDefinition(),benchmark)
    return plot(artifact,selection;ydata,kwargs...)
end
