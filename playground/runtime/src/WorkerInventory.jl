"""
    WorkerPresence

Retain a freshly challenged report and its local monotonic arrival time.
Persistence is deliberately excluded: coordinator restart begins with unknown
presence and does not restore warmth from a saved record.
"""
struct WorkerPresence
    "Validated report from an authenticated identity subject."
    report::Protocol.WorkerAnnouncement
    "Monotonic report arrival time in seconds."
    received_at::Float64
end

"""
    WorkerInventory(store, profiles; clock=()->time_ns()/1e9,
        heartbeat_seconds=2, presence_seconds=10, max_workers=128)

Track bounded, fresh worker presence separately from durable approval. The
transport supplies the broker-authenticated subject identity to accept_report!;
payload identity alone is never sufficient. The injectable monotonic clock is
for deterministic lifetime tests, not browser time.
"""
struct WorkerInventory{F}
    "Persistent operator approvals."
    store::RuntimeStore
    "Approved environment descriptions."
    profiles::ProfileRegistry
    "Fresh coordinator incarnation UUID."
    coordinator_id::String
    "Local monotonic clock in seconds."
    clock::F
    "Expected probe interval in seconds."
    heartbeat_seconds::Float64
    "Maximum usable report age in seconds."
    presence_seconds::Float64
    "Maximum tracked worker identities."
    max_workers::Int
    "Outstanding single-use challenge and its issue time, by worker."
    challenges::Dict{String,Tuple{String,Float64}}
    "Last accepted live report, by worker."
    presence::Dict{String,WorkerPresence}
    "Retired boot identities which cannot become current again."
    retired::Dict{String,Set{String}}
    "Serialize identity, freshness and incarnation changes."
    lock::ReentrantLock
end

function WorkerInventory(store::RuntimeStore, profiles::ProfileRegistry;
        clock=()->time_ns()/1e9, heartbeat_seconds=2, presence_seconds=10, max_workers=128)
    all(value -> value isa Real && !(value isa Bool) && isfinite(value) && value > 0,
        (heartbeat_seconds, presence_seconds)) && presence_seconds > 2heartbeat_seconds &&
        presence_seconds <= 60 || throw(ArgumentError("worker presence requires finite bounded lifetimes"))
    max_workers isa Integer && !(max_workers isa Bool) && 1 <= max_workers <= 4096 ||
        throw(ArgumentError("worker inventory bound must be in 1:4096"))
    return WorkerInventory(store, profiles, string(uuid4()), clock,
        Float64(heartbeat_seconds), Float64(presence_seconds), max_workers,
        Dict{String,Tuple{String,Float64}}(), Dict{String,WorkerPresence}(),
        Dict{String,Set{String}}(), ReentrantLock())
end

function inventory_registration(inventory::WorkerInventory, id::String)
    return lock(inventory.store.lock) do
        registered_worker(inventory.store.db, id)
    end
end

"""
    probe_worker!(inventory, worker_id) -> WorkerProbe

Issue one bounded freshness challenge for a registered, non-pending identity.
Replace its previous outstanding challenge without refreshing presence.
Publication and broker authentication belong to the transport adapter.
"""
function probe_worker!(inventory::WorkerInventory, worker_id::AbstractString)
    id = Protocol.runtime_token(worker_id)
    return lock(inventory.lock) do
        registration = inventory_registration(inventory, id)
        registration.state != :pending || throw(AccessDenied(409, "Worker approval is pending"))
        known = union(keys(inventory.presence), keys(inventory.challenges))
        id in known || length(known) < inventory.max_workers ||
            throw(CapacityUnavailable())
        challenge = string(uuid4())
        inventory.challenges[id] = (challenge, inventory.clock())
        return Protocol.WorkerProbe("2.0", id, inventory.coordinator_id, challenge)
    end
end

"""
    accept_report!(inventory, authenticated_worker_id, report) -> WorkerPresence

Accept an exact outstanding challenge from its broker-authorized identity.
Reject duplicate/delayed reports, old coordinators, retired boots, incompatible
environments and capacity inflation. A changed worker boot fences old assignments
through its changed identity; assignment checks must also consult this inventory.
"""
function accept_report!(inventory::WorkerInventory, authenticated_worker_id::AbstractString,
        report::Protocol.WorkerAnnouncement)
    id = Protocol.runtime_token(authenticated_worker_id)
    Protocol.validate(report)
    report.worker_id == id || throw(AccessDenied(403, "Worker subject and payload identities disagree"))
    return lock(inventory.lock) do
        registration = inventory_registration(inventory, id)
        registration.state != :pending || throw(AccessDenied(403, "Worker is not approved"))
        report.coordinator_id == inventory.coordinator_id ||
            throw(AccessDenied(409, "Worker report belongs to an old coordinator"))
        challenge = get(inventory.challenges, id, nothing)
        timestamp = inventory.clock()
        challenge !== nothing && challenge[1] == report.challenge &&
            0 <= timestamp - challenge[2] <= inventory.presence_seconds ||
            throw(AccessDenied(409, "Worker freshness challenge is missing or expired"))
        report.capacity <= registration.capacity ||
            throw(AccessDenied(409, "Worker exceeds approved capacity"))
        for installed in report.profiles
            installed.profile_id in registration.profiles &&
                haskey(inventory.profiles.definitions, installed.profile_id) ||
                throw(AccessDenied(409, "Worker profile is not approved"))
            approved = inventory.profiles.definitions[installed.profile_id]
            installed.version == string(approved.version) &&
                installed.fingerprint == approved.fingerprint ||
                throw(AccessDenied(409, "Worker environment does not match its approved profile"))
        end
        retired = get!(inventory.retired, id, Set{String}())
        report.boot_id in retired && throw(AccessDenied(409, "Worker boot was retired"))
        previous = get(inventory.presence, id, nothing)
        if previous !== nothing
            if previous.report.boot_id == report.boot_id
                report.sequence > previous.report.sequence ||
                    throw(AccessDenied(409, "Worker report sequence did not advance"))
            else
                length(retired) < 256 ||
                    throw(AccessDenied(409, "Worker incarnation history requires coordinator reconciliation"))
                push!(retired, previous.report.boot_id)
            end
        end
        presence = WorkerPresence(report, Float64(timestamp))
        inventory.presence[id] = presence
        delete!(inventory.challenges, id)
        return presence
    end
end

function presence_state(inventory::WorkerInventory, presence::Union{Nothing,WorkerPresence})
    presence === nothing && return :unknown
    age = inventory.clock() - presence.received_at
    age < 0 && return :unknown
    age > inventory.presence_seconds && return :offline
    return age > 2inventory.heartbeat_seconds ? :stale : :online
end

"""
    worker_inventory(inventory, principal)

Return approval and liveness as separate dimensions. Installed profiles do not
imply prepared executors. Offline rows remain visible for diagnostics and pinned
placement, but cannot supply a live assignment.
"""
function worker_inventory(inventory::WorkerInventory, principal::Principal)
    return lock(inventory.lock) do
        [(registration=worker_payload(registration),
          liveness=String(presence_state(inventory, get(inventory.presence, registration.worker_id, nothing))),
          report=public_worker_report(get(inventory.presence, registration.worker_id, nothing)))
         for registration in list_registrations(inventory.store, principal)]
    end
end

public_worker_report(::Nothing) = nothing
function public_worker_report(presence::WorkerPresence)
    report = presence.report
    return (boot_id=report.boot_id, sequence=report.sequence, capacity=report.capacity,
        profiles=[(profile_id=p.profile_id, version=p.version, fingerprint=p.fingerprint) for p in report.profiles])
end
