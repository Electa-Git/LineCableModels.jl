"""
    AgentLease

Track authority for one run/role on this worker incarnation. Closed records
retain a generation fence; they do not retain a process or prepared model.
"""
mutable struct AgentLease
    "Exact accepted assignment identity."
    fence::Protocol.AssignmentFence
    "Latest accepted control command."
    command::Protocol.LeaseControl
    "Acknowledgement for an idempotent retry, or nothing during cleanup."
    acknowledgement::Union{Nothing,Protocol.LeaseAcknowledgement}
    "Local monotonic authority deadline in seconds."
    expires_at::Float64
    "Active, closing or closed; closing still occupies capacity."
    state::Symbol
end

"""
    AgentLeaseLedger(worker_id, profiles; capacity=1, clock=()->time_ns()/1e9,
        presence_seconds=10, max_history=16384)

Enforce bounded worker-side lease authority without executing code. Expiry and
release make a lease unusable immediately, but keep its slot occupied until the
resource supervisor confirms cleanup. Generation tombstones are retained for
this boot; exhausting the finite history rejects new run/role keys rather than
forgetting fences. Restart creates a new boot, never remembered warmth.
"""
mutable struct AgentLeaseLedger{F}
    "Provisioned worker identity."
    worker_id::String
    "Fresh process incarnation."
    boot_id::String
    "Approved profiles, not loaded environments."
    profiles::ProfileRegistry
    "Maximum occupied leases."
    capacity::Int
    "Local monotonic clock."
    clock::F
    "Maximum time without a current coordinator probe, in seconds."
    presence_seconds::Float64
    "Maximum retained distinct run/role generation fences."
    max_history::Int
    "Committed coordinator incarnation, if reconciled."
    coordinator_id::Union{Nothing,String}
    "Most recent probe awaiting old-resource reconciliation."
    pending_probe::Union{Nothing,Protocol.WorkerProbe}
    "Monotonic time of the latest accepted probe."
    probed_at::Float64
    "Retired coordinator incarnations which cannot be revived."
    retired_coordinators::Set{String}
    "Current or closed assignment for each run/role."
    leases::Dict{Tuple{String,String},AgentLease}
    "Serialize controls, expiry and cleanup acknowledgements."
    lock::ReentrantLock
end

function AgentLeaseLedger(worker_id::AbstractString, profiles::ProfileRegistry;
        capacity=1, clock=()->time_ns()/1e9, presence_seconds=10, max_history=16384)
    capacity isa Integer && !(capacity isa Bool) && 1 <= capacity <= 256 ||
        throw(ArgumentError("agent capacity must be in 1:256"))
    max_history isa Integer && !(max_history isa Bool) && capacity <= max_history <= 65536 ||
        throw(ArgumentError("agent history bound must be between capacity and 65536"))
    presence_seconds isa Real && !(presence_seconds isa Bool) && isfinite(presence_seconds) &&
        0 < presence_seconds <= 60 || throw(ArgumentError("invalid agent coordinator-presence bound"))
    return AgentLeaseLedger(Protocol.runtime_token(worker_id), string(uuid4()), profiles,
        capacity, clock, Float64(presence_seconds), max_history, nothing, nothing,
        -Inf, Set{String}(), Dict{Tuple{String,String},AgentLease}(), ReentrantLock())
end

agent_key(fence::Protocol.AssignmentFence) = (fence.run_id, fence.role)
agent_ack(command::Protocol.LeaseControl, accepted::Bool, reason::String) =
    Protocol.LeaseAcknowledgement("2.0", command.request_id, command.fence, command.revision, accepted, reason)

function coordinator_current(ledger::AgentLeaseLedger)
    age = ledger.clock() - ledger.probed_at
    return ledger.coordinator_id !== nothing && ledger.pending_probe === nothing &&
        0 <= age <= ledger.presence_seconds
end

"""
    expire_agent_leases!(ledger) -> Vector{AssignmentFence}

Revoke starts on expired authority or lost coordinator presence. Return occupied
resources awaiting cleanup; do not claim their processes have stopped.
"""
function expire_agent_leases!(ledger::AgentLeaseLedger)
    return lock(ledger.lock) do
        timestamp = ledger.clock()
        current = coordinator_current(ledger)
        for lease in values(ledger.leases)
            if lease.state == :active && (!current || timestamp >= lease.expires_at)
                lease.state = :closing
            end
        end
        return [lease.fence for lease in values(ledger.leases) if lease.state == :closing]
    end
end

"""
    receive_probe!(ledger, probe) -> Bool

Accept a coordinator freshness probe. A new incarnation revokes old leases
before committing its identity; return false while owned cleanup is outstanding.
The agent may announce presence only after this returns true. A new coordinator
does not inherit prior leases or warm processes.
"""
function receive_probe!(ledger::AgentLeaseLedger, probe::Protocol.WorkerProbe)
    Protocol.validate(probe)
    probe.worker_id == ledger.worker_id || throw(AccessDenied(403, "Worker probe identity mismatch"))
    return lock(ledger.lock) do
        probe.coordinator_id in ledger.retired_coordinators &&
            throw(AccessDenied(409, "Coordinator incarnation was retired"))
        if ledger.coordinator_id == probe.coordinator_id && ledger.pending_probe === nothing
            # Expire against the previous probe before refreshing its clock.
            expire_agent_leases!(ledger)
            ledger.probed_at = ledger.clock()
            return true
        end
        if ledger.coordinator_id !== nothing && ledger.coordinator_id != probe.coordinator_id
            ledger.coordinator_id in ledger.retired_coordinators || length(ledger.retired_coordinators) < 1024 ||
                throw(AccessDenied(409, "Coordinator history requires an agent restart"))
            # Retire immediately, not after cleanup: a delayed old probe must
            # not cancel a newer incarnation's reconciliation.
            push!(ledger.retired_coordinators, ledger.coordinator_id)
            for lease in values(ledger.leases)
                lease.state == :active && (lease.state = :closing)
            end
        end
        ledger.pending_probe = probe
        any(lease -> lease.state != :closed, values(ledger.leases)) && return false
        if ledger.coordinator_id !== nothing && ledger.coordinator_id != probe.coordinator_id
            push!(ledger.retired_coordinators, ledger.coordinator_id)
        end
        ledger.coordinator_id = probe.coordinator_id
        ledger.pending_probe = nothing
        ledger.probed_at = ledger.clock()
        return true
    end
end

function compatible_agent_profile(ledger::AgentLeaseLedger, fence::Protocol.AssignmentFence)
    profile = get(ledger.profiles.definitions, fence.profile_id, nothing)
    return profile !== nothing && string(profile.version) == fence.profile_version &&
        profile.fingerprint == fence.fingerprint
end

"""
    handle_lease_control!(ledger, command) -> Union{Nothing,LeaseAcknowledgement}

Apply an authenticated coordinator command. Grants and renewals require a fresh
coordinator, exact worker boot/profile and non-stale generation/revision. Duplicate
current commands return their original acknowledgement without extending expiry.

Release returns nothing while resources require teardown. The resource-owning
supervisor must subsequently call complete_agent_cleanup! for the exact fence.
Older revisions are rejected; they never restore authority after release.
"""
function handle_lease_control!(ledger::AgentLeaseLedger, command::Protocol.LeaseControl)
    Protocol.validate(command)
    return lock(ledger.lock) do
        expire_agent_leases!(ledger)
        fence = command.fence
        reject(reason) = agent_ack(command, false, reason)
        fence.worker_id == ledger.worker_id && fence.worker_boot == ledger.boot_id ||
            return reject("wrong-worker-boot")
        fence.coordinator_id == ledger.coordinator_id && ledger.pending_probe === nothing ||
            return reject("wrong-coordinator")
        key = agent_key(fence)
        previous = get(ledger.leases, key, nothing)
        if previous !== nothing
            if previous.fence == fence
                previous.command.request_id == command.request_id && previous.command != command &&
                    return reject("request-id-reused")
                previous.command == command && return previous.acknowledgement
                command.revision > previous.command.revision || return reject("stale-revision")
            else
                fence.generation > previous.fence.generation || return reject("stale-generation")
                previous.state == :closed || return reject("role-busy")
            end
        elseif length(ledger.leases) >= ledger.max_history
            return reject("history-capacity")
        end
        same_lease = previous !== nothing && previous.fence == fence
        if command.action == "release"
            if same_lease
                previous.command = command
                if previous.state == :closed
                    previous.acknowledgement = agent_ack(command, true, "released")
                    return previous.acknowledgement
                end
                previous.state = :closing
                previous.acknowledgement = nothing
                return nothing
            end
            # An out-of-order release fences a grant that has not arrived. No
            # resource exists for this exact assignment, so cleanup is complete.
            ack = agent_ack(command, true, "released")
            ledger.leases[key] = AgentLease(fence, command, ack, 0.0, :closed)
            return ack
        end
        coordinator_current(ledger) || return reject("coordinator-stale")
        compatible_agent_profile(ledger, fence) || return reject("profile-mismatch")
        if command.action == "grant"
            same_lease && return reject("already-granted")
            command.revision == 1 || return reject("invalid-grant-revision")
            count(lease -> lease.state != :closed, values(ledger.leases)) < ledger.capacity ||
                return reject("capacity-unavailable")
            ack = agent_ack(command, true, "granted")
            ledger.leases[key] = AgentLease(fence, command, ack,
                Float64(ledger.clock() + command.duration_ms / 1000), :active)
            return ack
        end
        same_lease && previous.state == :active || return reject("lease-unavailable")
        command.revision == previous.command.revision + 1 || return reject("stale-revision")
        previous.command = command
        previous.expires_at = Float64(ledger.clock() + command.duration_ms / 1000)
        previous.acknowledgement = agent_ack(command, true, "renewed")
        return previous.acknowledgement
    end
end

"""
    complete_agent_cleanup!(ledger, fence) -> Union{Nothing,LeaseAcknowledgement}

Record supervisor-confirmed teardown for an exact closing assignment. This
function does not stop a process itself and must never precede actual cleanup.
Expiry has no reply; explicit release returns its matching acknowledgement.
"""
function complete_agent_cleanup!(ledger::AgentLeaseLedger, fence::Protocol.AssignmentFence)
    return lock(ledger.lock) do
        lease = get(ledger.leases, agent_key(fence), nothing)
        lease !== nothing && lease.fence == fence || throw(AccessDenied(409, "Cleanup fence is stale"))
        lease.state == :closed && return lease.command.action == "release" ? lease.acknowledgement : nothing
        lease.state == :closing || throw(AccessDenied(409, "Assignment is not awaiting cleanup"))
        lease.state = :closed
        if lease.command.action == "release"
            lease.acknowledgement = agent_ack(lease.command, true, "released")
            return lease.acknowledgement
        end
        return nothing
    end
end

"""
    agent_lease_usable(ledger, fence) -> Bool

Check live authority immediately before preparation or job execution. This does
not establish readiness, validate job inputs or authorize a terminal byte stream.
"""
function agent_lease_usable(ledger::AgentLeaseLedger, fence::Protocol.AssignmentFence)
    return lock(ledger.lock) do
        expire_agent_leases!(ledger)
        lease = get(ledger.leases, agent_key(fence), nothing)
        return lease !== nothing && lease.fence == fence && lease.state == :active &&
            coordinator_current(ledger)
    end
end
