"""
    WorkerTrust(worker_id, credential_ref, profiles; capacity=1)

Bind an operator-provisioned broker identity reference to its permitted profiles
and capacity ceiling. This declaration contains no secret and does not enroll a
worker. Only trusted configuration may construct the gateway's trust catalogue.
"""
struct WorkerTrust
    "Broker-authorized worker subject identity."
    worker_id::String
    "Opaque operator credential reference, never credential contents."
    credential_ref::String
    "Profiles this identity is permitted to announce."
    profiles::Tuple{Vararg{String}}
    "Maximum admitted executor slots."
    capacity::Int
    function WorkerTrust(worker_id, credential_ref, profiles; capacity=1)
        ids = Tuple(Protocol.runtime_token.(profiles))
        !isempty(ids) && length(ids) <= 64 && allunique(ids) ||
            throw(ArgumentError("worker trust requires distinct approved profiles"))
        capacity isa Integer && !(capacity isa Bool) && 1 <= capacity <= 256 ||
            throw(ArgumentError("worker capacity must be an integer in 1:256"))
        return new(Protocol.runtime_token(worker_id), checked_id(credential_ref), ids, capacity)
    end
end

"""
    WorkerRegistration

Return persistent operator approval, not current presence or executor warmth.
A disabled or draining registration rejects new assignments without implicitly
terminating existing work.
"""
struct WorkerRegistration
    "Approved broker subject identity."
    worker_id::String
    "Opaque operator credential reference."
    credential_ref::String
    "Approved profile identifiers."
    profiles::Tuple{Vararg{String}}
    "Approved capacity ceiling."
    capacity::Int
    "Pending, approved, draining or disabled."
    state::Symbol
    "Monotonic registration revision."
    revision::Int
    "UTC creation timestamp."
    created_at::DateTime
    "UTC last-change timestamp."
    updated_at::DateTime
end

function create_worker_schema!(db)
    sql_rows(db, """
        CREATE TABLE workers (
            worker_id TEXT PRIMARY KEY,
            credential_ref TEXT NOT NULL UNIQUE,
            profiles_json TEXT NOT NULL,
            capacity INTEGER NOT NULL CHECK(capacity BETWEEN 1 AND 256),
            state TEXT NOT NULL CHECK(state IN ('pending','approved','draining','disabled')),
            revision INTEGER NOT NULL CHECK(revision > 0),
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL
        )
    """)
    sql_rows(db, """
        CREATE TABLE worker_actions (
            owner TEXT NOT NULL,
            request_id TEXT NOT NULL,
            action TEXT NOT NULL,
            input_hash TEXT NOT NULL,
            response_json TEXT NOT NULL,
            PRIMARY KEY(owner, request_id)
        )
    """)
end

function check_worker_schema!(db)
    sql_rows(db, "SELECT worker_id,credential_ref,profiles_json,capacity,state,revision,created_at,updated_at FROM workers LIMIT 0")
    sql_rows(db, "SELECT owner,request_id,action,input_hash,response_json FROM worker_actions LIMIT 0")
end

worker_record(row) = WorkerRegistration(row.worker_id, row.credential_ref,
    Tuple(String.(JSON3.read(row.profiles_json))), row.capacity, Symbol(row.state),
    row.revision, DateTime(row.created_at), DateTime(row.updated_at))

worker_payload(worker::WorkerRegistration) = (
    worker_id=worker.worker_id, credential_ref=worker.credential_ref,
    profiles=worker.profiles, capacity=worker.capacity, state=String(worker.state),
    revision=worker.revision, created_at=string(worker.created_at), updated_at=string(worker.updated_at),
)

function require_administrator(principal::Principal)
    principal.administrator || throw(AccessDenied(403, "Worker administration is not permitted"))
end

function registered_worker(db, id::AbstractString)
    rows = sql_rows(db, "SELECT * FROM workers WHERE worker_id = ?", (String(id),))
    isempty(rows) && throw(AccessDenied(404, "Worker registration not found"))
    return worker_record(only(rows))
end

"""
    list_registrations(store, principal)

Read approved inventory metadata for an authenticated caller. No report in this
list establishes liveness, preparation or a usable assignment.
"""
function list_registrations(store::RuntimeStore, principal::Principal)
    return lock(store.lock) do
        worker_record.(sql_rows(store.db, "SELECT * FROM workers ORDER BY worker_id"))
    end
end

function registration_action!(action, store::RuntimeStore, principal::Principal,
        request_id::UUID, name::String, inputs)
    require_administrator(principal)
    digest = bytes2hex(sha256(JSON3.write(inputs)))
    return transaction(store) do
        prior = sql_rows(store.db, "SELECT * FROM worker_actions WHERE owner=? AND request_id=?",
            (principal.id, string(request_id)))
        if !isempty(prior)
            row = only(prior)
            row.action == name && row.input_hash == digest ||
                throw(AccessDenied(409, "Request identity already belongs to another control action"))
            value = JSON3.read(row.response_json)
            return WorkerRegistration(value.worker_id, value.credential_ref, Tuple(String.(value.profiles)),
                value.capacity, Symbol(value.state), value.revision,
                DateTime(value.created_at), DateTime(value.updated_at))
        end
        worker = action()
        sql_rows(store.db, """
            INSERT INTO worker_actions(owner,request_id,action,input_hash,response_json)
            VALUES (?,?,?,?,?)
        """, (principal.id, string(request_id), name, digest, JSON3.write(worker_payload(worker))))
        return worker
    end
end

"""
    enroll_worker!(store, principal, trust; request_id=uuid4())

Create a pending registration from an operator-provisioned identity binding.
Only administrators may enroll; discovery cannot approve itself. Duplicate
worker IDs or credential references are rejected unless the request is replayed
with identical inputs.
"""
function enroll_worker!(store::RuntimeStore, principal::Principal, trust::WorkerTrust;
        request_id::UUID=uuid4())
    input = (worker_id=trust.worker_id, credential_ref=trust.credential_ref,
        profiles=trust.profiles, capacity=trust.capacity)
    return registration_action!(store, principal, request_id, "enroll", input) do
        isempty(sql_rows(store.db, "SELECT worker_id FROM workers WHERE worker_id=? OR credential_ref=?",
            (trust.worker_id, trust.credential_ref))) ||
            throw(AccessDenied(409, "Worker identity or credential reference is already registered"))
        timestamp = string(now(UTC))
        sql_rows(store.db, """
            INSERT INTO workers(worker_id,credential_ref,profiles_json,capacity,state,revision,created_at,updated_at)
            VALUES (?,?,?,?,'pending',1,?,?)
        """, (trust.worker_id, trust.credential_ref, JSON3.write(trust.profiles),
            trust.capacity, timestamp, timestamp))
        return registered_worker(store.db, trust.worker_id)
    end
end

"""
    set_registration_state!(store, principal, worker_id, state;
        expected_revision, request_id=uuid4())

Apply explicit approval, drain or disable under optimistic revision checking.
Idempotent retries return the original response; delayed actions cannot undo a
newer operator decision. This action does not kill processes or revoke secrets.
"""
function set_registration_state!(store::RuntimeStore, principal::Principal,
        worker_id::AbstractString, state::Symbol; expected_revision::Integer,
        request_id::UUID=uuid4())
    id = Protocol.runtime_token(worker_id)
    state in (:approved, :draining, :disabled) || throw(ArgumentError("invalid registration action"))
    expected_revision isa Bool && throw(ArgumentError("revision must be an integer"))
    Protocol.runtime_sequence(expected_revision)
    input = (worker_id=id, state=String(state), expected_revision=expected_revision)
    return registration_action!(store, principal, request_id, "state", input) do
        worker = registered_worker(store.db, id)
        worker.revision == expected_revision ||
            throw(AccessDenied(409, "Registration changed; refresh before applying another action"))
        worker.state == :pending && state == :draining &&
            throw(AccessDenied(409, "A pending registration cannot be drained"))
        worker.state == state && return worker
        sql_rows(store.db, "UPDATE workers SET state=?,revision=revision+1,updated_at=? WHERE worker_id=?",
            (String(state), string(now(UTC)), id))
        return registered_worker(store.db, id)
    end
end

"""
    migrate_runtime!(path) -> Union{Nothing,String}

Migrate a stopped earlier runtime database (schemas 1–4) to the current schema. Acquire the
same exclusive supervisor lock, create a private consistent backup before any
schema change, then migrate transactionally. Return its backup path, or nothing
when already current. Never overwrite an existing backup or migrate a live
coordinator's database.
"""
function migrate_runtime!(path::AbstractString)
    absolute = abspath(path)
    Sys.islinux() || throw(ArgumentError("runtime migration requires Linux"))
    isfile(absolute) && !islink(absolute) || throw(ArgumentError("runtime database must exist and not be linked"))
    lockpath = absolute * ".supervisor.lock"
    islink(lockpath) && throw(ArgumentError("supervisor lock cannot be linked"))
    ownership = open(lockpath, "a+")
    if ccall(:flock, Cint, (Cint, Cint), fd(ownership), 6) != 0
        close(ownership)
        throw(ArgumentError("stop the runtime before migrating its database"))
    end
    db = nothing
    try
        db = SQLite.DB(absolute)
        sql_rows(db, "PRAGMA busy_timeout = 5000")
        version = only(sql_rows(db, "PRAGMA user_version")).user_version
        version in (1, 2, 3, 4, RUN_SCHEMA_VERSION) || throw(ArgumentError("unsupported runtime database schema"))
        sql_rows(db, "SELECT run_id,owner,application,application_version,state,request_id,created_at,updated_at,reason FROM runs LIMIT 0")
        if version == RUN_SCHEMA_VERSION
            check_worker_schema!(db)
            check_lease_schema!(db)
            check_job_schema!(db)
            check_job_cancellation_schema!(db)
            return nothing
        end
        version >= 2 && check_worker_schema!(db)
        version >= 3 && check_lease_schema!(db)
        version >= 4 && check_job_schema!(db)
        directory = mktempdir(dirname(absolute); prefix=basename(absolute) * ".schema-$version-backup-", cleanup=false)
        backup = joinpath(directory, "runtime.sqlite")
        # A bound SQLite busy timeout replaces an indefinitely retrying backup
        # loop. INTO writes a new snapshot and does not vacuum the source.
        try
            sql_rows(db, "PRAGMA synchronous = FULL")
            sql_rows(db, "VACUUM INTO ?", (backup,))
            chmod(backup, 0o600)
        catch
            # Only this freshly allocated backup directory is eligible. Once
            # the snapshot is complete, retain it even if migration fails.
            islink(directory) || rm(directory; recursive=true)
            rethrow()
        end
        store = RuntimeStore(db, absolute, ReentrantLock())
        transaction(store) do
            only(sql_rows(db, "PRAGMA user_version")).user_version == version ||
                throw(ArgumentError("database version changed during migration"))
            version == 1 && create_worker_schema!(db)
            version < 3 && create_lease_schema!(db)
            version < 4 && create_job_schema!(db)
            create_job_cancellation_schema!(db)
            sql_rows(db, "PRAGMA user_version = 5")
        end
        return backup
    finally
        db === nothing || close(db)
        close(ownership)
    end
end
