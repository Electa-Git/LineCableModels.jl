"""
    AgentService(config, resources)

Own worker control and lease expiry independently of resource cleanup. Approved
resource hooks must be complete at construction. No user operation is evaluated
in this scheduler; execution belongs to the separate owned resource supervisor.
"""
mutable struct AgentService{L<:AgentLeaseLedger,R<:AbstractAgentResources}
    "Operator-owned agent settings."
    config::AgentConfig
    "Worker-incarnation lease authority."
    ledger::L
    "Restricted physical-resource owner."
    resources::R
    "Separate preparation/status/cancel owner, absent for non-scientific resource fixtures."
    science::Union{Nothing,AgentScientificService}
    "Separate bounded durable job consumer, absent for non-scientific fixtures."
    jobs::Union{Nothing,AgentJobService}
    "Independent private terminal relay, absent when the resource owner has none."
    terminals::Union{Nothing,AgentTerminalService}
    "Worker-role control link."
    link::ControlLink{WorkerIdentity}
    "Owned control scheduler."
    task::Union{Nothing,Task}
    "Independent finite connection attempt."
    connector::Union{Nothing,Task}
    "Cleanup tasks keyed by exact lease UUID, bounded by occupied capacity."
    cleanups::Dict{String,Task}
    "Next cleanup retry times for unresolved resources."
    cleanup_retry::Dict{String,Float64}
    "Monotonic announcement sequence for this worker boot."
    sequence::Int
    "Next initial broker connection attempt."
    retry_at::Float64
    "Current diagnostic state, not executor readiness."
    state::Symbol
    "Number of rejected control records or unresolved cleanup attempts."
    rejected::Int
    "Stop accepting control and revoke owned authority."
    closed::Bool
    "Whether physical teardown completed, distinct from rejecting new work."
    cleanup_complete::Bool
    "Serialize lifecycle and connection ownership."
    lock::ReentrantLock
    "Join repeated concurrent teardown calls through complete cleanup."
    shutdown_lock::ReentrantLock
end

function AgentService(config::AgentConfig, resources::AbstractAgentResources; clock=()->time_ns()/1e9)
    profiles = verified_agent_profiles(config, resources)
    ledger = AgentLeaseLedger(config.worker_id, profiles; capacity=config.capacity, clock)
    bind_agent!(resources, ledger) === nothing || throw(ArgumentError("agent resource authority binding failed"))
    scientific=scientific_resources(resources)
    science = scientific === nothing ? nothing : AgentScientificService(config.endpoint,scientific,ledger)
    jobs = scientific === nothing ? nothing : AgentJobService(config.endpoint,scientific,ledger;artifacts=config.artifacts)
    terminal=terminal_resources(resources)
    terminals=terminal === nothing ? nothing : AgentTerminalService(config.endpoint,terminal,ledger)
    return AgentService(config, ledger, resources, science, jobs, terminals, ControlLink(WorkerIdentity), nothing, nothing,
        Dict{String,Task}(), Dict{String,Float64}(), 0, 0.0, :offline, 0, false, false, ReentrantLock(), ReentrantLock())
end

function reject_agent_record!(agent::AgentService)
    agent.rejected == typemax(Int) || (agent.rejected += 1)
end

function connect_agent!(agent::AgentService)
    lock(agent.lock) do
        agent.closed && return
        agent.link.control === nothing || return
        agent.connector !== nothing && !istaskdone(agent.connector) && return
        agent.ledger.clock() >= agent.retry_at || return
        agent.retry_at = agent.ledger.clock() + 2
        agent.connector = @async begin
            connection = nothing
            try
                connection = BrokerControl(agent.config.endpoint, WorkerIdentity(agent.config.worker_id))
                lock(agent.lock) do
                    if !agent.closed
                        lock(agent.link.lock) do
                            agent.link.control = connection
                        end
                        connection = nothing
                    end
                end
            catch
                agent.state = :offline
            finally
                connection === nothing || close(connection)
            end
        end
    end
end

function announce_agent!(agent::AgentService, probe::Protocol.WorkerProbe)
    agent.sequence = Protocol.runtime_sequence(agent.sequence + 1)
    profiles = [Protocol.ProfileAdvertisement(p.id, string(p.version), p.fingerprint)
        for p in sort!(collect(values(agent.ledger.profiles.definitions)); by=p -> p.id)]
    report = Protocol.WorkerAnnouncement("2.0", agent.config.worker_id, agent.ledger.boot_id,
        probe.coordinator_id, probe.challenge, agent.sequence, agent.config.capacity, profiles)
    send_control!(agent.link, report)
end

function receive_agent_record!(agent::AgentService, record::Protocol.WorkerProbe)
    if receive_probe!(agent.ledger, record)
        announce_agent!(agent, record)
    else
        agent.state = :reconciling
    end
end
function receive_agent_record!(agent::AgentService, record::Protocol.LeaseControl)
    acknowledgement = handle_lease_control!(agent.ledger, record)
    acknowledgement === nothing || send_control!(agent.link, acknowledgement)
end

function poll_agent_cleanup!(agent::AgentService)
    # Expiry first revokes new starts, independently of whether cleanup can run.
    pending = expire_agent_leases!(agent.ledger)
    for (id, task) in collect(agent.cleanups)
        istaskdone(task) || continue
        delete!(agent.cleanups, id)
        success = !istaskfailed(task) && fetch(task) === true
        fence = findfirst(f -> f.lease_id == id, pending)
        fence === nothing && continue
        if success
            acknowledgement = complete_agent_cleanup!(agent.ledger, pending[fence])
            delete!(agent.cleanup_retry, id)
            if acknowledgement !== nothing
                try
                    send_control!(agent.link, acknowledgement)
                catch error
                    error isa BrokerUnavailable || rethrow()
                    # Closed lease retains its exact ACK for a retried release.
                end
            end
        else
            reject_agent_record!(agent)
            agent.cleanup_retry[id] = agent.ledger.clock() + 1
        end
    end
    for fence in pending
        haskey(agent.cleanups, fence.lease_id) && continue
        # A just-confirmed cleanup no longer occupies the ledger.
        lease = get(agent.ledger.leases, agent_key(fence), nothing)
        lease !== nothing && lease.state == :closing || continue
        agent.ledger.clock() >= get(agent.cleanup_retry, fence.lease_id, 0.0) || continue
        length(agent.cleanups) < agent.config.capacity || continue
        agent.cleanups[fence.lease_id] = @async try
            release_owned!(agent.resources, fence) === true
        catch
            false # no exception object or arbitrary text enters control diagnostics
        end
    end
    probe = agent.ledger.pending_probe
    if probe !== nothing && all(l -> l.state == :closed, values(agent.ledger.leases))
        receive_probe!(agent.ledger, probe) && announce_agent!(agent, probe)
    end
end

"""
    tick_agent!(agent)

Run one bounded control/expiry pass. Physical cleanup proceeds in separately
owned tasks; no acknowledgement or replacement announcement can precede it.
"""
function tick_agent!(agent::AgentService)
    poll_agent_cleanup!(agent)
    control = lock(() -> agent.link.control, agent.link.lock)
    control === nothing && return
    connected = !control.closed && NATS.status(control.connection) == NATS.CONNECTED
    agent.state = connected ? (agent.ledger.pending_probe === nothing ? :online : :reconciling) : :offline
    for envelope in poll_control!(control)
        try
            receive_agent_record!(agent, envelope.record)
        catch error
            error isa AccessDenied || error isa ArgumentError || error isa BrokerUnavailable || rethrow()
            reject_agent_record!(agent)
        end
    end
end

"""
    start_agent!(agent) -> AgentService

Recover owned resources before any announcement, then start independent control
scheduling. Repeated starts do not repeat recovery. A closed incarnation cannot
restart; construct a new agent so old leases and warm state remain fenced.
"""
function start_agent!(agent::AgentService)
    lock(agent.lock) do
        agent.closed && throw(ArgumentError("agent is closed"))
        agent.task !== nothing && return agent
        recover_owned!(agent.resources) === nothing ||
            throw(ArgumentError("agent resource recovery did not complete"))
        compile_runtime_paths(agent)
        agent.science===nothing || start_science!(agent.science)
        agent.jobs===nothing || start_jobs!(agent.jobs)
        agent.terminals===nothing || start_terminals!(agent.terminals)
        agent.task = @async begin
            while !agent.closed
                try
                    connect_agent!(agent)
                    tick_agent!(agent)
                catch error
                    error isa InterruptException && rethrow()
                    agent.state = :offline
                    reject_agent_record!(agent)
                end
                sleep(0.05)
            end
        end
    end
    return agent
end

Base.close(agent::AgentService) = lock(() -> close_agent_service!(agent), agent.shutdown_lock)

function close_agent_service!(agent::AgentService)
    first_close = lock(agent.lock) do
        agent.closed && return false
        agent.closed = true
        return true
    end
    if !first_close
        agent.cleanup_complete && return nothing
        # A previous failure keeps ownership. Retry the same owned teardown
        # when the physical driver can now finish; never reopen admission.
    end
    try
        # A failed scheduler is already joined; its exception must not prevent
        # resource retirement. The CLI still reports non-interrupt failures.
        agent.task === nothing || wait(agent.task;throw=false)
        agent.connector === nothing || wait(agent.connector;throw=false)
        lock(agent.ledger.lock) do
            for lease in values(agent.ledger.leases)
                lease.state == :active && (lease.state = :closing)
            end
        end
        agent.science===nothing || stop_science!(agent.science)
        agent.jobs===nothing || stop_jobs!(agent.jobs)
        agent.terminals===nothing || stop_terminals!(agent.terminals)
        # The resource owner has the same bounded stop primitive as release;
        # close also covers partial starts and resources predating this boot.
        close(agent.resources)
        agent.terminals===nothing || close(agent.terminals)
        agent.jobs===nothing || close(agent.jobs)
        agent.science===nothing || close(agent.science)
        foreach(wait, values(agent.cleanups))
        empty!(agent.cleanups)
        for fence in expire_agent_leases!(agent.ledger)
            if release_owned!(agent.resources, fence) === true
                acknowledgement = complete_agent_cleanup!(agent.ledger, fence)
                if acknowledgement !== nothing
                    try
                        send_control!(agent.link, acknowledgement)
                    catch error
                        error isa BrokerUnavailable || rethrow()
                    end
                end
            end
        end
        all(lease -> lease.state == :closed, values(agent.ledger.leases)) ||
            throw(ArgumentError("Agent resource cleanup remains unresolved"))
        agent.cleanup_complete = true
    finally
        control = lock(agent.link.lock) do
            previous = agent.link.control
            agent.link.control = nothing
            previous
        end
        try
            control === nothing || close(control)
        finally
            agent.state = :stopped
        end
    end
    return nothing
end
