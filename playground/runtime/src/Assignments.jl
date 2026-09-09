"""
    AssignmentManager(inventory, applications; limits=AssignmentLimits())

Reserve compatible workers transactionally. The inventory owns current
incarnations; SQLite owns capacity. A reservation remains unusable until the
separate control lifecycle accepts the matching acknowledgement.
"""
struct AssignmentManager{I<:WorkerInventory}
    "Challenged worker presence and its durable approval store."
    inventory::I
    "Immutable application definitions and their permitted roles."
    applications::ApplicationRegistry
    "Global, owner and run allocation limits."
    limits::AssignmentLimits
end
AssignmentManager(inventory::WorkerInventory, applications::ApplicationRegistry;
    limits=AssignmentLimits()) = AssignmentManager(inventory, applications, limits)

function assignment_requirement(manager::AssignmentManager, run::RunRecord, role::String, profile::String)
    definition = get(manager.applications.definitions, run.application, nothing)
    definition !== nothing && definition.version == run.version ||
        throw(AccessDenied(409, "Run application version is no longer registered"))
    requirements = filter(item -> item.role == role, definition.requirements)
    length(requirements) == 1 && profile in only(requirements).profiles ||
        throw(AccessDenied(400, "Profile is not permitted for this application role"))
    haskey(manager.inventory.profiles.definitions, profile) ||
        throw(AccessDenied(409, "Required profile is not installed in the coordinator"))
    return manager.inventory.profiles.definitions[profile]
end

function allocation_candidate(manager::AssignmentManager, occupied, run::RunRecord,
        profile::ProfileDefinition, placement::AbstractPlacement)
    inventory = manager.inventory
    pinned = placement_worker(placement)
    candidates = Tuple{Int,String,WorkerPresence}[]
    for (id, presence) in inventory.presence
        pinned === nothing || id == pinned || continue
        presence_state(inventory, presence) == :online || continue
        registration = registered_worker(inventory.store.db, id)
        registration.state == :approved && profile.id in registration.profiles || continue
        any(p -> p.profile_id == profile.id && p.version == string(profile.version) &&
            p.fingerprint == profile.fingerprint, presence.report.profiles) || continue
        assigned = filter(item -> item.worker_id == id, occupied)
        length(assigned) < min(registration.capacity, presence.report.capacity) || continue
        any(item -> item.placement == "dedicated" && item.run_id != string(run.id), assigned) && continue
        placement isa DedicatedPlacement && any(item -> item.run_id != string(run.id), assigned) && continue
        push!(candidates, (length(assigned), id, presence))
    end
    isempty(candidates) && throw(AccessDenied(409, pinned === nothing ?
        "No compatible approved worker has available capacity" :
        "Pinned worker is unavailable, incompatible or fully allocated"))
    sort!(candidates; by=item -> (item[1], item[2]))
    return first(candidates)[2:3]
end

"""
    reserve_assignment!(manager, principal, run_id, role, profile;
        placement=AutomaticPlacement(), request_id=uuid4()) -> LeaseRecord

Authorize the run and declared role, require fresh approved compatible presence,
and reserve capacity using one SQLite IMMEDIATE transaction. Pinned placement
never falls back. A dedicated reservation excludes other runs while permitting
sibling roles of its own run. It counts toward all limits before a grant is sent.

Idempotent retries return the current record for their original request. Reusing
a request ID with changed inputs fails. Persisted active/unreconciled rows remain
occupied after restart; only explicit lifecycle reconciliation can release them.
"""
function reserve_assignment!(manager::AssignmentManager, principal::Principal,
        run_id::UUID, role::AbstractString, profile::AbstractString;
        placement::AbstractPlacement=AutomaticPlacement(), request_id::UUID=uuid4())
    role_id, profile_id = Protocol.runtime_token(role), Protocol.runtime_token(profile)
    input = (run_id=string(run_id), role=role_id, profile=profile_id,
        placement=placement_kind(placement), worker=placement_worker(placement))
    fingerprint = bytes2hex(sha256(JSON3.write(input)))
    inventory = manager.inventory
    return lock(inventory.lock) do
        transaction(inventory.store) do
            db = inventory.store.db
            run = owned_run(db, principal, run_id)
            prior = sql_rows(db, "SELECT * FROM leases WHERE owner=? AND request_id=?",
                (run.owner, string(request_id)))
            if !isempty(prior)
                row = only(prior)
                row.input_hash == fingerprint ||
                    throw(AccessDenied(409, "Assignment request ID was already used with different inputs"))
                return lease_record(row)
            end
            # Administrators may inspect/control another run, but allocations
            # remain charged to its actual owner rather than the operator.
            run.state in (:reserved, :starting, :running) ||
                throw(AccessDenied(409, "Run is not accepting assignments"))
            definition = assignment_requirement(manager, run, role_id, profile_id)
            occupied = sql_rows(db, "SELECT * FROM leases WHERE state IN $OCCUPIED_LEASE_SQL")
            any(item -> item.run_id == string(run_id) && item.role == role_id, occupied) &&
                throw(AccessDenied(409, "Role already has an unreleased assignment"))
            limits = manager.limits
            length(occupied) < limits.total &&
                count(item -> item.owner == run.owner, occupied) < limits.per_owner &&
                count(item -> item.run_id == string(run_id), occupied) < limits.per_run ||
                throw(CapacityUnavailable())
            worker, presence = allocation_candidate(manager, occupied, run, definition, placement)
            previous = only(sql_rows(db, "SELECT COALESCE(MAX(generation),0) AS generation FROM leases WHERE run_id=? AND role=?",
                (string(run_id), role_id))).generation
            generation = Protocol.runtime_sequence(previous + 1)
            fence = Protocol.AssignmentFence(string(uuid4()), string(run_id), run.owner, role_id,
                worker, presence.report.boot_id, inventory.coordinator_id, profile_id,
                string(definition.version), definition.fingerprint, generation)
            Protocol.validate(fence)
            timestamp = string(now(UTC))
            sql_rows(db, """
                INSERT INTO leases(lease_id,owner,run_id,role,worker_id,worker_boot,coordinator_id,generation,
                    fence_json,placement,state,revision,request_id,input_hash,created_at,updated_at)
                VALUES (?,?,?,?,?,?,?,?,?,?,'reserving',0,?,?,?,?)
            """, (fence.lease_id, run.owner, fence.run_id, role_id, worker, fence.worker_boot,
                fence.coordinator_id, generation, Protocol.encode_message(fence), placement_kind(placement),
                string(request_id), fingerprint, timestamp, timestamp))
            return owned_lease(db, principal, UUID(fence.lease_id))
        end
    end
end
