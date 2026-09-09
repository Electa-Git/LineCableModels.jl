"""Track one bounded pending control and its acknowledged local authority."""
mutable struct LeaseFlight
    "Last exact command; retries keep its request ID and revision."
    command::Protocol.LeaseControl
    "Local monotonic issuance time."
    issued_at::Float64
    "Deadline derived from issuance, never a browser or worker clock."
    candidate_until::Float64
    "End of previously acknowledged authority; zero means unusable."
    authority_until::Float64
    "Whether the current command still needs an acknowledgement."
    pending::Bool
    "Next release retry time."
    retry_at::Float64
end

"""
    LeaseCoordinator(assignments, transport; duration_ms=10000, ack_seconds=2)

Join durable reservations to the control transport. A lease is usable only after
an exact acknowledgement and while its local authority deadline and challenged
worker incarnation remain valid. Lost acknowledgements retain capacity until an
explicit release acknowledgement or reconciled replacement report arrives.
"""
mutable struct LeaseCoordinator{M<:AssignmentManager,T}
    "Transactional allocation authority."
    assignments::M
    "Owned control transport implementing send_control!."
    transport::T
    "Bounded grant duration in milliseconds."
    duration_ms::Int
    "Maximum acknowledgement wait before release, in seconds."
    ack_seconds::Float64
    "Volatile current-incarnation authority, never restored from SQLite."
    flights::Dict{UUID,LeaseFlight}
    "First local observation of a reservation not yet issued to a worker."
    unissued::Dict{UUID,Float64}
    "Serialize controls and acknowledgements."
    lock::ReentrantLock
end

function LeaseCoordinator(assignments::AssignmentManager, transport; duration_ms=10000, ack_seconds=2)
    duration_ms isa Integer && !(duration_ms isa Bool) && 100 <= duration_ms <= 60_000 ||
        throw(ArgumentError("lease duration must be in 100:60000 milliseconds"))
    ack_seconds isa Real && !(ack_seconds isa Bool) && isfinite(ack_seconds) &&
        0 < ack_seconds < duration_ms / 1000 ||
        throw(ArgumentError("acknowledgement bound must be shorter than the lease"))
    return LeaseCoordinator(assignments, transport, duration_ms, Float64(ack_seconds),
        Dict{UUID,LeaseFlight}(), Dict{UUID,Float64}(), ReentrantLock())
end

lease_clock(coordinator::LeaseCoordinator) = coordinator.assignments.inventory.clock()

function update_lease_state!(db, lease::LeaseRecord, state::Symbol; revision=lease.revision)
    sql_rows(db, "UPDATE leases SET state=?,revision=?,updated_at=? WHERE lease_id=? AND revision=?",
        (String(state), revision, string(now(UTC)), lease.fence.lease_id, lease.revision))
    only(sql_rows(db, "SELECT changes() AS count")).count == 1 ||
        throw(AccessDenied(409, "Assignment control revision changed"))
end

function current_assignment_worker(coordinator::LeaseCoordinator, fence::Protocol.AssignmentFence)
    inventory = coordinator.assignments.inventory
    return lock(inventory.lock) do
        presence = get(inventory.presence, fence.worker_id, nothing)
        fence.coordinator_id == inventory.coordinator_id && presence !== nothing &&
            presence.report.boot_id == fence.worker_boot && presence_state(inventory, presence) == :online &&
            any(profile -> profile.profile_id == fence.profile_id && profile.version == fence.profile_version &&
                profile.fingerprint == fence.fingerprint, presence.report.profiles)
    end
end

function send_lease_command!(coordinator::LeaseCoordinator, flight::LeaseFlight)
    try
        send_control!(coordinator.transport, flight.command)
        return true
    catch error
        error isa BrokerUnavailable || rethrow()
        # Transport uncertainty does not free the reservation or extend time.
        return false
    end
end

function issue_lease_control!(coordinator::LeaseCoordinator, principal::Principal, id::UUID,
        action::String; request_id::UUID=uuid4())
    action in ("grant", "renew", "release") || throw(ArgumentError("invalid lease action"))
    return lock(coordinator.lock) do
        inventory = coordinator.assignments.inventory
        lock(inventory.lock) do
            result = transaction(inventory.store) do
                db = inventory.store.db
                lease = owned_lease(db, principal, id)
                lease.fence.coordinator_id == inventory.coordinator_id ||
                    throw(AccessDenied(409, "Assignment requires worker reconciliation"))
                previous = get(coordinator.flights, id, nothing)
                if previous !== nothing && previous.command.request_id == string(request_id)
                    previous.command.action == action ||
                        throw(AccessDenied(409, "Control request ID was reused"))
                    return previous
                end
                if action == "release"
                    lease.state in (:released, :expired, :failed) && return nothing
                    previous !== nothing && previous.command.action == "release" && return previous
                else
                    run = owned_run(db, principal, UUID(lease.fence.run_id))
                    run.state in (:reserved, :starting, :running) ||
                        throw(AccessDenied(409, "Run is not accepting work"))
                    current_assignment_worker(coordinator, lease.fence) ||
                        throw(AccessDenied(409, "Assigned worker is not currently available"))
                    action == "grant" && registered_worker(db, lease.fence.worker_id).state != :approved &&
                        throw(AccessDenied(409, "Worker is not accepting new assignments"))
                    if action == "grant"
                        lease.state == :reserving && lease.revision == 0 ||
                            throw(AccessDenied(409, "Assignment was already granted or revoked"))
                    else
                        lease.state == :active && previous !== nothing && !previous.pending &&
                            previous.authority_until > lease_clock(coordinator) ||
                            throw(AccessDenied(409, "Assignment has no renewable authority"))
                    end
                end
                timestamp = Float64(lease_clock(coordinator))
                revision = Protocol.runtime_sequence(lease.revision + 1)
                duration = action == "release" ? 0 : coordinator.duration_ms
                command = Protocol.LeaseControl("2.0", string(request_id), action, lease.fence, revision, duration)
                Protocol.validate(command)
                state = action == "release" ? :releasing : lease.state
                update_lease_state!(db, lease, state; revision)
                authority = action == "renew" ? previous.authority_until : 0.0
                return LeaseFlight(command, timestamp, timestamp + duration / 1000,
                    authority, true, timestamp + coordinator.ack_seconds)
            end
            result === nothing && return nothing
            coordinator.flights[id] = result
            # Only initial issuance or a still-current explicit retry is sent.
            result.pending && send_lease_command!(coordinator, result)
            return result.command
        end
    end
end

"""Send a reserved assignment's first bounded grant; it is not usable yet."""
grant_assignment!(coordinator::LeaseCoordinator, principal::Principal, id::UUID; request_id=uuid4()) =
    issue_lease_control!(coordinator, principal, id, "grant"; request_id)

"""Renew acknowledged authority without changing its run/role/generation fence."""
renew_assignment!(coordinator::LeaseCoordinator, principal::Principal, id::UUID; request_id=uuid4()) =
    issue_lease_control!(coordinator, principal, id, "renew"; request_id)

"""Revoke local use immediately and await resource-confirmed release."""
release_assignment!(coordinator::LeaseCoordinator, principal::Principal, id::UUID; request_id=uuid4()) =
    issue_lease_control!(coordinator, principal, id, "release"; request_id)

"""
    accept_lease_ack!(coordinator, authenticated_worker_id, acknowledgement) -> Bool

Accept only the outstanding command's exact subject identity, fence, request and
revision. Delayed grant/renewal acknowledgements cannot revive expired authority.
Release acknowledgements may finish cleanup after a transport outage.
"""
function accept_lease_ack!(coordinator::LeaseCoordinator, worker_id::AbstractString,
        acknowledgement::Protocol.LeaseAcknowledgement)
    Protocol.validate(acknowledgement)
    worker = Protocol.runtime_token(worker_id)
    return lock(coordinator.lock) do
        id = UUID(acknowledgement.fence.lease_id)
        flight = get(coordinator.flights, id, nothing)
        flight !== nothing && flight.pending || return false
        command = flight.command
        command.fence.worker_id == worker && acknowledgement.fence == command.fence &&
            acknowledgement.request_id == command.request_id && acknowledgement.revision == command.revision ||
            return false
        inventory = coordinator.assignments.inventory
        return lock(inventory.lock) do
            transaction(inventory.store) do
                db = inventory.store.db
                lease = owned_lease(db, Principal(command.fence.owner), id)
                lease.revision == command.revision && lease.fence == command.fence || return false
                timestamp = lease_clock(coordinator)
                if command.action != "release"
                    timestamp - flight.issued_at <= coordinator.ack_seconds &&
                        timestamp < flight.candidate_until &&
                        current_assignment_worker(coordinator, command.fence) || return false
                end
                if acknowledgement.accepted
                    state = command.action == "release" ? :released : :active
                    update_lease_state!(db, lease, state)
                    flight.authority_until = command.action == "release" ? 0.0 : flight.candidate_until
                else
                    # Even a rejected/uncertain grant is reconciled via release;
                    # no negative message alone claims physical cleanup.
                    update_lease_state!(db, lease, :reconciling)
                    flight.authority_until = 0.0
                end
                flight.pending = false
                return true
            end
        end
    end
end

"""
    assignment_usable(coordinator, principal, lease_id) -> Bool

Require owned durable state, acknowledged unexpired local authority, an exact
current worker boot and an approved registration immediately before a new start.
This check does not imply that the assigned executor is prepared.
"""
function assignment_usable(coordinator::LeaseCoordinator, principal::Principal, id::UUID)
    return lock(coordinator.lock) do
        inventory = coordinator.assignments.inventory
        lock(inventory.lock) do
            lease = get_assignment(inventory.store, principal, id)
            flight = get(coordinator.flights, id, nothing)
            lease.state == :active && flight !== nothing &&
                flight.authority_until > lease_clock(coordinator) &&
                current_assignment_worker(coordinator, lease.fence) || return false
            run = get_run(inventory.store, principal, UUID(lease.fence.run_id))
            return run.state in (:reserved, :starting, :running) &&
                inventory_registration(inventory, lease.fence.worker_id).state == :approved
        end
    end
end

"""
    tick_leases!(coordinator)

Expire local authority, request cleanup on missing acknowledgements, and retry
only the exact release command. This tick does not automatically prepare work
or renew indefinitely after a run closes. Abandoned unissued reservations have
a bounded local grace; closed runs request release immediately.
"""
function tick_leases!(coordinator::LeaseCoordinator)
    lock(coordinator.lock) do
        timestamp = lease_clock(coordinator)
        inventory = coordinator.assignments.inventory
        occupied = lock(inventory.store.lock) do
            sql_rows(inventory.store.db, """
                SELECT leases.*,runs.state AS run_state FROM leases JOIN runs USING(run_id)
                WHERE leases.coordinator_id=? AND leases.state IN $OCCUPIED_LEASE_SQL
            """, (inventory.coordinator_id,))
        end
        occupied_ids = Set{UUID}(UUID(row.lease_id) for row in occupied)
        filter!(pair -> first(pair) in occupied_ids, coordinator.unissued)
        for row in occupied
            id = UUID(row.lease_id)
            flight = get(coordinator.flights, id, nothing)
            run_closed = !(row.run_state in ("reserved", "starting", "running"))
            if flight === nothing
                observed = get!(coordinator.unissued, id, Float64(timestamp))
                if run_closed || row.revision > 0 || timestamp - observed >= coordinator.ack_seconds || timestamp < observed
                    release_assignment!(coordinator, Principal(row.owner), id)
                    delete!(coordinator.unissued, id)
                end
            elseif run_closed && flight.command.action != "release"
                release_assignment!(coordinator, Principal(row.owner), id)
            end
        end
        for (id, flight) in collect(coordinator.flights)
            lease = get_assignment(coordinator.assignments.inventory.store,
                Principal(flight.command.fence.owner), id)
            if lease.state in (:released, :expired, :failed)
                delete!(coordinator.flights, id)
                continue
            end
            if flight.command.action == "release"
                if timestamp >= flight.retry_at
                    flight.retry_at = timestamp + coordinator.ack_seconds
                    flight.pending = true
                    send_lease_command!(coordinator, flight)
                end
                continue
            end
            expired = (flight.pending && timestamp - flight.issued_at >= coordinator.ack_seconds) ||
                (!flight.pending && timestamp >= flight.authority_until) ||
                lease.state == :reconciling
            if expired
                @debug "Lease authority revoked after its finite control bound" worker=flight.command.fence.worker_id lease_id=flight.command.fence.lease_id action=flight.command.action revision=flight.command.revision pending=flight.pending age_seconds=timestamp-flight.issued_at authority_remaining_seconds=flight.authority_until-timestamp
                release_assignment!(coordinator, Principal(lease.fence.owner), id)
            end
        end
    end
    return nothing
end

"""
    reconcile_worker_report!(coordinator, authenticated_worker_id, report)

Accept a challenged report, then retire allocations from replaced worker/coordinator
incarnations. The agent protocol requires completed owned-resource reconciliation
before announcing a new incarnation. Same-incarnation reports never erase leases.
"""
function reconcile_worker_report!(coordinator::LeaseCoordinator, worker_id::AbstractString,
        report::Protocol.WorkerAnnouncement)
    lock(coordinator.lock) do
        inventory = coordinator.assignments.inventory
        lock(inventory.lock) do
            accept_report!(inventory, worker_id, report)
            transaction(inventory.store) do
                sql_rows(inventory.store.db, """
                    UPDATE leases SET state='expired',updated_at=?
                    WHERE worker_id=? AND state IN $OCCUPIED_LEASE_SQL
                        AND (worker_boot<>? OR coordinator_id<>?)
                """, (string(now(UTC)), report.worker_id, report.boot_id, report.coordinator_id))
            end
        end
    end
    return nothing
end
