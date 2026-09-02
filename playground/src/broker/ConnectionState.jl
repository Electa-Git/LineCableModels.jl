const BROKER_STATES = (:offline, :connecting, :online, :degraded, :closed)

struct WorkerPresence
    heartbeat::WorkerHeartbeat
    received_at::Float64
end

function worker_capabilities(workers::Dict{String,WorkerPresence})
    result = Dict{String,Vector{String}}()
    for (worker_id, presence) in workers
        result[worker_id] = copy(presence.heartbeat.capabilities)
    end
    return result
end

function supports_operation(
        workers::Dict{String,WorkerPresence},
        operation::AbstractString;
        engine_constraint::Union{Nothing,AbstractString}=nothing
    )
    return any(
        operation in presence.heartbeat.capabilities &&
        (isnothing(engine_constraint) ||
         presence.heartbeat.engine_version == engine_constraint)
        for presence in values(workers)
    )
end
