"""
$(TYPEDSIGNATURES)

Tabulate completed timing records as scalar columns. No timings are collected
and no files are opened. Allocated bytes are Julia allocations, not peak memory.
Native backend timing scopes are preserved rather than relabelled as wall time.

The optional `labels` are supplied by the same formulation descriptions used
by the comparison report. Raw workload/session records stay in its publication.
"""
function tabulate(::BenchmarkTableDefinition, measurements::Union{Nothing, NamedTuple};
        labels=(reference="Reference",candidate="Candidate"))
    execution=DataFrame()
    source_timings=DataFrame()
    performance=DataFrame()
    performance_samples=DataFrame()
    performance_environment=DataFrame()
    performance_policy=DataFrame()
    performance_comparison=DataFrame()
    if measurements !== nothing
        for (role, record) in pairs(get(measurements,:execution,(;)))
            method=getproperty(labels,role)
            timing=get(record,:timing,(;))
            session=get(record,:session,nothing)
            push!(execution,(;role,method,point=missing,scope=get(timing,:scope,missing),
                seconds=get(timing,:seconds,missing),reused=get(record,:reused,missing),
                session_id=session===nothing ? missing : get(session,:id,missing),
                execution_wall_seconds=get(record,:execution_wall_seconds,missing),
                reused_points=get(timing,:reused_points,missing));cols=:union)
            native=get(timing,:source_timings,(;))
            sources=haskey(native,:points) ? enumerate(native.points) : ((1,native),)
            for (point,measured) in sources
                isempty(measured) && continue
                # Source timing records are already scoped by their measurement
                # owner. Non-scalar annotations are text, never hidden arrays.
                scalar=(; (key => (value===nothing ? missing :
                    value isa Union{Number,Bool,Symbol,AbstractString,Missing} ? value : string(value))
                    for (key,value) in pairs(measured))...)
                push!(source_timings,merge((;role,method,point),scalar);cols=:union)
            end
            for (point,measured) in enumerate(get(timing,:points,()))
                push!(execution,(;role,method,point,scope=get(measured,:scope,missing),
                    seconds=get(measured,:seconds,missing),reused=get(measured,:reused,missing));cols=:union)
            end
        end
        recorded=get(measurements,:performance,nothing)
        if recorded !== nothing
            session=get(measurements,:session,nothing)
            for role in (:reference,:candidate)
                method=getproperty(labels,role)
                measured=getproperty(recorded,role)
                observations=get(measured,:observations,())
                policy=get(measured,:policy,(;))
                push!(performance,(;role,method,scope=measured.scope,
                    median_seconds=measured.median_seconds,allocated_bytes=measured.bytes,
                    allocated_MiB=measured.bytes/2.0^20,
                    allocation_statistic=get(policy,:allocation_statistic,
                        isempty(observations) ? :not_recorded : :maximum),
                    allocation_scope=get(policy,:allocation_scope,missing),
                    samples=measured.samples,requested_samples=recorded.settings.samples,
                    reused=any(row -> get(row,:reused,false),observations),
                    session_id=session===nothing ? missing : get(session,:id,missing),
                    checksum_verified=get(measurements,:checksum_verified,missing),
                    workload_verified=get(measurements,:workload_verified,missing));cols=:union)
                for (sample,observation) in enumerate(observations)
                    push!(performance_samples,(;role,method,sample,scope=measured.scope,
                        seconds=get(observation,:seconds,missing),
                        allocated_bytes=get(observation,:bytes,missing),
                        reused=get(observation,:reused,missing));cols=:union)
                end
                for (table,values) in ((performance_environment,measured.environment),
                        (performance_policy,policy))
                    for (key,value) in pairs(values)
                        # Full calculation/workload bindings remain retained,
                        # not printed as one opaque default table cell.
                        key in (:calculation,:workload,:settings) && continue
                        value isa NamedTuple && continue
                        push!(table,(;role,method,setting=string(key),
                            value=value===nothing ? missing :
                                value isa Union{Number,Bool,Symbol,AbstractString,Missing} ? value : string(value));cols=:union)
                    end
                end
            end
            push!(performance_comparison,(
                reference_over_candidate=recorded.speedup,comparable=recorded.comparable,
                passes=something(recorded.passes,missing),
                requested_samples=recorded.settings.samples,time_budget_seconds=recorded.settings.seconds,
                minimum_speedup=recorded.settings.minimum_speedup))
        end
    end
    return (;execution,source_timings,performance,performance_samples,
        performance_environment,performance_policy,performance_comparison)
end
