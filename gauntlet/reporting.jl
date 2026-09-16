import LineCableModels.PlotBuilder: plot

"""
    select(definition::BenchmarkTableDefinition, benchmark)

Verify saved operands and decode retained comparisons into the shared report
representation. Reading performs no numerical comparison or solver execution.
"""
function select(definition::BenchmarkTableDefinition,
        benchmark::NamedTuple{(:id, :reference, :candidate, :analyses, :measurements)})
    isempty(benchmark.analyses) && throw(ArgumentError("benchmark has no saved analyses"))
    for operand in (benchmark.reference,benchmark.candidate)
        validate(read_calculation, operand.result, operand.metadata)
    end
    comparisons=NamedTuple[]
    for record in benchmark.analyses
        for role in (:reference,:candidate)
            recorded=getproperty(record["calculations"],role)
            loaded=getproperty(benchmark,role).metadata
            recorded.sha256 == loaded.sha256 || throw(ArgumentError("analysis $role differs from the loaded operand"))
            recorded.port_order == loaded.port_order && recorded.frequencies == loaded.frequencies &&
                recorded.basis == loaded.basis && recorded.domain == loaded.domain ||
                throw(ArgumentError("analysis $role coordinates differ from the loaded operand"))
        end
        snapshot=get(record,"analysis_id",semantic_sha256(record["comparison_settings"]))
        for row in record["reference_comparison"]
            T=Base.nonmissingtype(eltype(row.absolute))
            error=RMSError{T}(row.absolute,row.relative;details=Grammar.ComputationDetails(row.details))
            push!(comparisons,(snapshot,request=row.request,quantity=row.quantity,statistic=row.statistic,
                reference_index=row.reference_index,candidate_index=row.candidate_index,error))
        end
    end
    first_record=first(benchmark.analyses)
    settings=first_record["comparison_settings"]
    if length(benchmark.analyses)>1
        selected_keys=(:requests,:bands,:normalizations)
        settings=merge(settings,NamedTuple{selected_keys}(Tuple(Tuple(unique(Iterators.flatten(
            getproperty(record["comparison_settings"],key) for record in benchmark.analyses))) for key in selected_keys)))
        for key in definition.requested
            key in (:requests,:bands,:normalizations) && continue
            all(record -> isequal(getproperty(record["comparison_settings"],key),getproperty(definition.settings,key)),benchmark.analyses) ||
                throw(ArgumentError("saved analyses use different $key; select a snapshot file explicitly"))
        end
    end
    measurements=benchmark.measurements
    return select(definition,(reference=benchmark.reference,candidate=benchmark.candidate,
        context=(id=benchmark.id,case_id=Symbol(first_record["case_id"]),collection=Symbol(first_record["collection"]),
            analyses=Tuple((id=get(record,"analysis_id",semantic_sha256(record["comparison_settings"])),
                settings=record["comparison_settings"]) for record in benchmark.analyses)),
        settings,comparisons,measurements))
end

"""Plot explicitly requested saved comparisons through the shared report recipe."""
function plot(benchmark::NamedTuple{(:id, :reference, :candidate, :analyses, :measurements)},
        selection=nothing; ydata=nothing, kwargs...)
    artifact=report(BenchmarkTableDefinition(),benchmark)
    selection === nothing || ydata === nothing || throw(ArgumentError(
        "use either positional ydata or the ydata keyword, not both"))
    selected_ydata=ydata === nothing ? something(selection,(Z,Y)) : ydata
    return plot(artifact,selected_ydata;kwargs...)
end
