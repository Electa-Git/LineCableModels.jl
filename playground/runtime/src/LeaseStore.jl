"""
    AbstractPlacement

Choose automatic, pinned-worker or dedicated-run allocation without exposing
broker subjects or launch commands to callers.
"""
abstract type AbstractPlacement end
"""Select any currently compatible worker with available approved capacity."""
struct AutomaticPlacement <: AbstractPlacement end
"""Select one worker; unavailability never falls back to another worker."""
struct PinnedPlacement <: AbstractPlacement
    "Required worker identity."
    worker_id::String
    PinnedPlacement(id::AbstractString) = new(Protocol.runtime_token(id))
end
"""Reserve a worker exclusively for this run, optionally pinning its identity."""
struct DedicatedPlacement <: AbstractPlacement
    "Required worker identity, or nothing for automatic dedicated placement."
    worker_id::Union{Nothing,String}
    DedicatedPlacement(id::Union{Nothing,AbstractString}=nothing) =
        new(id === nothing ? nothing : Protocol.runtime_token(id))
end
placement_kind(::AutomaticPlacement) = "automatic"
placement_kind(::PinnedPlacement) = "pinned"
placement_kind(::DedicatedPlacement) = "dedicated"
placement_worker(::AutomaticPlacement) = nothing
placement_worker(value::Union{PinnedPlacement,DedicatedPlacement}) = value.worker_id

"""
    AssignmentLimits(; total=128, per_owner=8, per_run=4)

Bound simultaneous reserved, active or unreconciled assignments. Counts include
leases awaiting acknowledgement or cleanup; persisted history never frees a
resource merely because a coordinator restarted.
"""
struct AssignmentLimits
    "Maximum occupied assignments across the coordinator."
    total::Int
    "Maximum occupied assignments for one authenticated owner."
    per_owner::Int
    "Maximum occupied assignments for one application run."
    per_run::Int
    function AssignmentLimits(; total=128, per_owner=8, per_run=4)
        all(v -> v isa Integer && !(v isa Bool), (total, per_owner, per_run)) &&
            1 <= per_run <= per_owner <= total <= 4096 ||
            throw(ArgumentError("assignment limits require 1 <= per_run <= per_owner <= total <= 4096"))
        new(total, per_owner, per_run)
    end
end

"""Describe durable allocation bookkeeping, not live executor authority."""
struct LeaseRecord
    "Exact run, role, environment and incarnation fence."
    fence::Protocol.AssignmentFence
    "Automatic, pinned or dedicated placement."
    placement::Symbol
    "Reserving, active, releasing, reconciling, released, expired or failed."
    state::Symbol
    "Last persisted control revision."
    revision::Int
    "Idempotent reservation request UUID."
    request_id::UUID
    "UTC creation timestamp for history, not lease expiry."
    created_at::DateTime
    "UTC transition timestamp for history, not lease expiry."
    updated_at::DateTime
end

const OCCUPIED_LEASE_STATES = ("reserving", "active", "releasing", "reconciling")
const OCCUPIED_LEASE_SQL = "('reserving','active','releasing','reconciling')"

function create_lease_schema!(db)
    sql_rows(db, """
        CREATE TABLE leases (
            lease_id TEXT PRIMARY KEY,
            owner TEXT NOT NULL,
            run_id TEXT NOT NULL REFERENCES runs(run_id),
            role TEXT NOT NULL,
            worker_id TEXT NOT NULL REFERENCES workers(worker_id),
            worker_boot TEXT NOT NULL,
            coordinator_id TEXT NOT NULL,
            generation INTEGER NOT NULL CHECK(generation>0),
            fence_json TEXT NOT NULL,
            placement TEXT NOT NULL CHECK(placement IN ('automatic','pinned','dedicated')),
            state TEXT NOT NULL CHECK(state IN ('reserving','active','releasing','reconciling','released','expired','failed')),
            revision INTEGER NOT NULL CHECK(revision>=0),
            request_id TEXT NOT NULL,
            input_hash TEXT NOT NULL,
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL,
            UNIQUE(owner,request_id),
            UNIQUE(run_id,role,generation)
        )
    """)
    sql_rows(db, "CREATE INDEX leases_worker_state ON leases(worker_id,state)")
    sql_rows(db, """
        CREATE UNIQUE INDEX leases_live_role ON leases(run_id,role)
        WHERE state IN ('reserving','active','releasing','reconciling')
    """)
end

check_lease_schema!(db) = sql_rows(db, """
    SELECT lease_id,owner,run_id,role,worker_id,worker_boot,coordinator_id,generation,
        fence_json,placement,state,revision,request_id,input_hash,created_at,updated_at FROM leases LIMIT 0
""")

function lease_record(row)
    fence = Protocol.decode_runtime_message(Protocol.AssignmentFence, row.fence_json)
    (fence.lease_id, fence.owner, fence.run_id, fence.role, fence.worker_id, fence.worker_boot,
        fence.coordinator_id, fence.generation) ==
        (row.lease_id, row.owner, row.run_id, row.role, row.worker_id, row.worker_boot,
         row.coordinator_id, row.generation) || throw(ArgumentError("inconsistent persisted assignment fence"))
    return LeaseRecord(fence, Symbol(row.placement), Symbol(row.state), row.revision,
        UUID(row.request_id), DateTime(row.created_at), DateTime(row.updated_at))
end

function owned_lease(db, principal::Principal, id::UUID)
    rows = sql_rows(db, "SELECT * FROM leases WHERE lease_id=?", (string(id),))
    isempty(rows) && throw(AccessDenied(404, "Assignment not found"))
    row = only(rows)
    row.owner == principal.id || principal.administrator || throw(AccessDenied(404, "Assignment not found"))
    return lease_record(row)
end

"""
    get_assignment(store, principal, lease_id) -> LeaseRecord

Read owned allocation bookkeeping. Foreign and missing assignments have the same
404 response. An active row alone is not permission to execute a job.
"""
get_assignment(store::RuntimeStore, principal::Principal, id::UUID) =
    lock(() -> owned_lease(store.db, principal, id), store.lock)

"""
    list_assignments(store, principal; run_id=nothing)

List caller-owned allocation history, or all records for an explicit operator.
A supplied run is authorized before querying its assignments.
"""
function list_assignments(store::RuntimeStore, principal::Principal; run_id::Union{Nothing,UUID}=nothing)
    return lock(store.lock) do
        if run_id !== nothing
            owned_run(store.db, principal, run_id)
            return lease_record.(sql_rows(store.db,
                "SELECT * FROM leases WHERE run_id=? ORDER BY created_at,lease_id", (string(run_id),)))
        end
        rows = principal.administrator ?
            sql_rows(store.db, "SELECT * FROM leases ORDER BY created_at,lease_id") :
            sql_rows(store.db, "SELECT * FROM leases WHERE owner=? ORDER BY created_at,lease_id", (principal.id,))
        return lease_record.(rows)
    end
end

