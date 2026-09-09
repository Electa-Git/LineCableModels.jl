"""Connect a stable lease coordinator to a replaceable, owned control connection."""
mutable struct ControlLink{I<:AbstractBrokerIdentity}
    "Current role-specific control connection, or nothing while unavailable."
    control::Union{Nothing,BrokerControl{I}}
    "Protect connection publication and replacement."
    lock::ReentrantLock
end
ControlLink() = ControlLink(CoordinatorIdentity)
ControlLink(::Type{I}) where {I<:AbstractBrokerIdentity} = ControlLink{I}(nothing, ReentrantLock())

function send_control!(link::ControlLink, record::Protocol.RuntimeRecord)
    lock(link.lock) do
        link.control === nothing && throw(BrokerUnavailable())
        send_control!(link.control, record)
    end
end

"""
    ControlService(config, store, applications)

Construct inert worker control ownership around the existing SQLite store.
Call start_control! explicitly to schedule bounded connection/probe/lease work.
No scientific import, preparation, evaluation or container launch occurs here.
"""
mutable struct ControlService{I<:WorkerInventory,C<:LeaseCoordinator}
    "Validated server-owned worker/profile configuration."
    config::ControlConfig
    "Volatile challenged presence, separate from approval."
    inventory::I
    "Durable reservations and acknowledged lease authority."
    coordinator::C
    "Stable transient sender shared by the coordinator."
    link::ControlLink{CoordinatorIdentity}
    "Bounded redacted diagnostic history."
    events::ControlEvents
    "Independent preparation/status/cancel channel sharing the same lease authority."
    science::ScientificCoordinator
    "Independent durable job/cancellation/result owner."
    jobs::JobCoordinator
    "Independent private terminal request/reply owner."
    terminals::TerminalCoordinator
    "Owned scheduling task."
    task::Union{Nothing,Task}
    "Separate finite initial-connect task; it never holds the state lock while connecting."
    connector::Union{Nothing,Task}
    "Current connection dimension, not scientific readiness."
    state::Symbol
    "Next allowed initial connection attempt, in monotonic seconds."
    retry_at::Float64
    "Next inventory probe cycle, in monotonic seconds."
    probe_at::Float64
    "Idempotent shutdown flag."
    closed::Bool
    "Serialize start, stop and connection ownership."
    lock::ReentrantLock
    "Join repeated concurrent teardown calls through complete cleanup."
    shutdown_lock::ReentrantLock
end

function ControlService(config::ControlConfig, store::RuntimeStore, applications::ApplicationRegistry;
        clock=()->time_ns()/1e9)
    # Existing trust is immutable here. Configuration cannot silently replace a
    # persisted credential binding or increase previously approved capacity.
    for worker in list_registrations(store, Principal("runtime-control"; administrator=true))
        provisioned = get(config.workers, worker.worker_id, nothing)
        provisioned === nothing && continue # retained history; no subscriptions or new presence
        (worker.credential_ref, worker.profiles, worker.capacity) ==
            (provisioned.credential_ref, provisioned.profiles, provisioned.capacity) ||
            throw(ArgumentError("persisted worker trust differs from operator configuration"))
    end
    inventory = WorkerInventory(store, config.profiles; clock)
    link = ControlLink()
    coordinator = LeaseCoordinator(AssignmentManager(inventory, applications; limits=config.limits), link)
    science = ScientificCoordinator(config.endpoint,coordinator,keys(config.workers))
    events = ControlEvents()
    jobs = JobCoordinator(config,coordinator,science,events)
    terminals=TerminalCoordinator(config.endpoint,coordinator,keys(config.workers))
    return ControlService(config, inventory, coordinator, link, events, science, jobs, terminals,
        nothing, nothing, :unavailable, 0.0, 0.0, false, ReentrantLock(), ReentrantLock())
end

function control_state!(service::ControlService, state::Symbol)
    lock(service.lock) do
        service.state == state && return
        service.state = state
        record_event!(service.events, state == :online ? :control_connected :
            state == :stopped ? :control_stopped : :control_unavailable)
    end
end

function connect_control!(service::ControlService)
    lock(service.lock) do
        service.closed && return
        service.link.control === nothing || return
        service.connector !== nothing && !istaskdone(service.connector) && return
        service.inventory.clock() >= service.retry_at || return
        service.retry_at = service.inventory.clock() + 2
        service.connector = @async begin
            connection = nothing
            try
                connection = BrokerControl(service.config.endpoint, CoordinatorIdentity();
                    worker_ids=keys(service.config.workers))
                lock(service.lock) do
                    if !service.closed
                        lock(service.link.lock) do
                            service.link.control = connection
                        end
                        connection = nothing # ownership transferred exactly once
                        service.probe_at = 0.0
                    end
                end
            catch
                control_state!(service, :unavailable)
            finally
                connection === nothing || close(connection)
            end
        end
    end
end

function accept_control_record!(service::ControlService, envelope::ControlEnvelope{Protocol.WorkerAnnouncement})
    prior = get(service.inventory.presence, envelope.worker_id, nothing)
    reconcile_worker_report!(service.coordinator, envelope.worker_id, envelope.record)
    if prior === nothing || prior.report.boot_id != envelope.record.boot_id
        record_event!(service.events, :worker_reported; worker_id=envelope.worker_id)
    end
end
function accept_control_record!(service::ControlService, envelope::ControlEnvelope{Protocol.LeaseAcknowledgement})
    if accept_lease_ack!(service.coordinator, envelope.worker_id, envelope.record)
        record_event!(service.events, :assignment_acknowledged; fence=envelope.record.fence)
    end
end

function renew_live_assignments!(service::ControlService)
    coordinator = service.coordinator
    lock(coordinator.lock) do
        timestamp = lease_clock(coordinator)
        for (id, flight) in collect(coordinator.flights)
            flight.pending && continue
            flight.command.action == "release" && continue
            # Renew at half-life, only for a still-owned live run and unchanged
            # challenged worker. Draining rejects new work, not existing leases.
            0 < flight.authority_until - timestamp <= coordinator.duration_ms / 2000 || continue
            fence = flight.command.fence
            principal = Principal(fence.owner)
            lease = get_assignment(service.inventory.store, principal, id)
            run = get_run(service.inventory.store, principal, UUID(fence.run_id))
            run.state in (:reserved, :starting, :running) && lease.state == :active || continue
            current_assignment_worker(coordinator, fence) || continue
            renew_assignment!(coordinator, principal, id)
        end
    end
end

"""
    tick_control!(service; renew=true)

Perform one bounded control poll, probe cycle and lease-expiry pass. Every record
is validated independently. Broker unavailability does not suspend local expiry;
scientific and terminal output are never handled by this scheduler.
"""
function tick_control!(service::ControlService; renew::Bool=true)
    # Local revocation cannot depend on a successful connection poll or report.
    tick_leases!(service.coordinator)
    control = lock(() -> service.link.control, service.link.lock)
    if control !== nothing
        connected = !control.closed && NATS.status(control.connection) == NATS.CONNECTED
        control_state!(service, connected ? :online : :unavailable)
        for envelope in poll_control!(control)
            try
                accept_control_record!(service, envelope)
            catch error
                error isa AccessDenied || error isa ArgumentError || rethrow()
                record_event!(service.events, :control_rejected; worker_id=envelope.worker_id)
            end
        end
        if connected && service.inventory.clock() >= service.probe_at
            service.probe_at = service.inventory.clock() + service.inventory.heartbeat_seconds
            for registration in list_registrations(service.inventory.store, Principal("runtime-control"))
                registration.state == :pending && continue
                haskey(service.config.workers, registration.worker_id) || continue
                try
                    send_control!(service.link, probe_worker!(service.inventory, registration.worker_id))
                catch error
                    error isa BrokerUnavailable || rethrow()
                    control_state!(service, :unavailable)
                    break
                end
            end
        end
        connected && renew && renew_live_assignments!(service)
    end
    return nothing
end

"""
    start_control!(service) -> ControlService

Start the owned coordinator scheduler without waiting for NATS. The public
gateway can start immediately. Repeated starts are inert; a closed service cannot
be restarted with its retired coordinator incarnation.
"""
function start_control!(service::ControlService)
    lock(service.lock) do
        service.closed && throw(ArgumentError("control service is closed"))
        service.task !== nothing && return service
        compile_runtime_paths(service)
        start_science!(service.science)
        start_jobs!(service.jobs)
        start_terminals!(service.terminals)
        service.task = @async begin
            while !service.closed
                try
                    connect_control!(service)
                    tick_control!(service)
                catch
                    # No raw broker/HTTP/command object is exposed in diagnostics.
                    control_state!(service, :unavailable)
                end
                sleep(0.05)
            end
        end
    end
    return service
end

Base.close(service::ControlService) = lock(() -> close_control_service!(service), service.shutdown_lock)

function close_control_service!(service::ControlService)
    first_close = lock(service.lock) do
        service.closed && return false
        service.closed = true
        return true
    end
    first_close || return nothing
    service.task === nothing || wait(service.task)
    service.connector === nothing || wait(service.connector)
    close(service.jobs)
    close(service.science)
    close(service.terminals)
    # Revocation is immediate; persisted capacity remains occupied if the
    # remote supervisor cannot acknowledge its physical cleanup during shutdown.
    try
        for lease in list_assignments(service.inventory.store, Principal("runtime-control"; administrator=true))
            lease.fence.coordinator_id == service.inventory.coordinator_id || continue
            lease.state in (:reserving, :active, :releasing, :reconciling) || continue
            release_assignment!(service.coordinator, Principal(lease.fence.owner), UUID(lease.fence.lease_id))
        end
        control = lock(() -> service.link.control, service.link.lock)
        deadline = time_ns() / 1e9 + service.coordinator.ack_seconds
        while control !== nothing && NATS.status(control.connection) == NATS.CONNECTED &&
                any(f -> f.pending, values(service.coordinator.flights)) && time_ns() / 1e9 < deadline
            for envelope in poll_control!(control)
                envelope.record isa Protocol.LeaseAcknowledgement || continue
                accept_control_record!(service, envelope)
            end
            sleep(0.02)
        end
    finally
        control = lock(service.link.lock) do
            previous = service.link.control
            service.link.control = nothing
            previous
        end
        try
            control === nothing || close(control)
        finally
            control_state!(service, :stopped)
        end
    end
    return nothing
end
