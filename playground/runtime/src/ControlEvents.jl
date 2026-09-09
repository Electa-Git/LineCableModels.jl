"""
    ControlEvents(; capacity=512)

Keep bounded, structured diagnostics. Only registered event codes and validated
identities are accepted; arbitrary exception text, credentials, and terminal input
cannot enter this buffer. The sequence exposes eviction/reconnect gaps.
"""
mutable struct ControlEvents
    "Maximum retained records."
    capacity::Int
    "Fresh buffer incarnation, so a reconnect cannot mistake restarted sequences for continuity."
    epoch::String
    "Next monotonically increasing event sequence."
    sequence::Int
    "Structured records, oldest first."
    records::Vector{NamedTuple}
    "Serialize appends and snapshots."
    lock::ReentrantLock
    function ControlEvents(; capacity::Integer=512)
        !(capacity isa Bool) && 1 <= capacity <= 4096 || throw(ArgumentError("invalid event capacity"))
        new(capacity, string(uuid4()), 0, NamedTuple[], ReentrantLock())
    end
end

const CONTROL_EVENT_CODES = (
    :control_connected, :control_unavailable, :control_rejected, :control_stopped,
    :worker_reported, :worker_enrolled, :worker_registration_changed,
    :assignment_reserved, :assignment_releasing, :assignment_acknowledged,
    :job_queued, :job_submitted, :job_succeeded, :job_failed, :job_uncertain,
    :job_canceled, :job_revoked, :job_cancel_requested,
)

"""
    record_event!(events, code; fence=nothing, job=nothing, worker_id=nothing, owner=nothing)

Append a fixed diagnostic code and optional ownership/fence context. Return its
sequence. Callers must supply context already authorized by the root orchestrator.
"""
function record_event!(events::ControlEvents, code::Symbol;
        fence::Union{Nothing,Protocol.AssignmentFence}=nothing,
        job::Union{Nothing,Protocol.AssignedJob}=nothing,
        worker_id::Union{Nothing,AbstractString}=nothing,
        owner::Union{Nothing,AbstractString}=nothing)
    code in CONTROL_EVENT_CODES || throw(ArgumentError("unregistered control event"))
    startswith(String(code),"job_") == (job!==nothing) || throw(ArgumentError("Job events require exact job context"))
    if job!==nothing
        Protocol.validate(job)
        fence===nothing || fence==job.fence || throw(ArgumentError("Job event fence differs"))
        fence=job.fence
    end
    worker_id === nothing || Protocol.runtime_token(worker_id)
    owner === nothing || Principal(owner)
    fence === nothing || Protocol.validate(fence)
    return lock(events.lock) do
        events.sequence < typemax(Int) || throw(ArgumentError("event sequence exhausted"))
        events.sequence += 1
        item = (sequence=events.sequence, at=string(now(UTC)), code=String(code),
            owner=fence === nothing ? owner : fence.owner,
            worker_id=fence === nothing ? worker_id : fence.worker_id,
            run_id=fence === nothing ? nothing : fence.run_id,
            worker_boot=fence === nothing ? nothing : fence.worker_boot,
            lease_id=fence === nothing ? nothing : fence.lease_id,
            generation=fence === nothing ? nothing : fence.generation,
            job_id=job===nothing ? nothing : job.request.job_id,
            executor_id=job===nothing ? nothing : job.execution.executor_id,
            executor_generation=job===nothing ? nothing : job.execution.executor_generation,
            stage=job===nothing ? nothing : String(code)[5:end])
        length(events.records) == events.capacity && popfirst!(events.records)
        push!(events.records, item)
        return events.sequence
    end
end

"""
    control_events(events, principal; after=0, epoch=nothing)

Return only public worker events and the caller's run events (all for an
administrator). Report buffer eviction with a gap flag, never silently imply a
complete historical log. Owner identity is not part of the public event payload.
"""
function control_events(events::ControlEvents, principal::Principal; after::Integer=0,
        epoch::Union{Nothing,AbstractString}=nothing)
    !(after isa Bool) && after >= 0 || throw(ArgumentError("invalid event cursor"))
    epoch === nothing || Protocol.runtime_uuid(epoch)
    return lock(events.lock) do
        replaced = epoch !== nothing && epoch != events.epoch
        cursor = replaced ? 0 : after
        oldest = isempty(events.records) ? events.sequence + 1 : first(events.records).sequence
        # Do not let the first nonempty/event-specific batch specialize the
        # response container and JSON writer while a lease ACK is in flight.
        records = NamedTuple[(; (k=>v for (k,v) in pairs(item) if k != :owner)...)
            for item in events.records if item.sequence > cursor &&
                (item.owner === nothing || item.owner == principal.id || principal.administrator)]
        return (epoch=events.epoch, cursor=events.sequence,
            gap=replaced || cursor < oldest - 1 || cursor > events.sequence,
            evicted=max(0, oldest - 1), records=records)
    end
end
