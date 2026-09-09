"""
    CapacityUnavailable()

Report that admission would exceed the configured global or owner limit.
"""
struct CapacityUnavailable <: Exception end
Base.showerror(io::IO, ::CapacityUnavailable) = print(io, "Application capacity is unavailable")

"""
    RunRecord

Return the durable, passive description of an application run. This record does
not prove that its process is live or prepared.
"""
struct RunRecord
    "Globally unique run identity."
    id::UUID
    "Authenticated owner identity."
    owner::String
    "Registered application identity."
    application::String
    "Application contract version at admission."
    version::VersionNumber
    "Durable lifecycle state; not a heartbeat."
    state::Symbol
    "Idempotent launch request identity."
    request_id::UUID
    "UTC creation timestamp."
    created_at::DateTime
    "UTC last-transition timestamp."
    updated_at::DateTime
    "Sanitized lifecycle reason."
    reason::String
end

"""
    RuntimeStore(path)

Open current SQLite run, worker, lease and job bookkeeping. Earlier schemas require the explicit
backed-up migration command; unknown schemas and non-runtime databases are
rejected without modification. Call close when finished.
SQLite transactions serialize admission across connections as well as tasks.
"""
struct RuntimeStore
    "Owned SQLite connection."
    db::SQLite.DB
    "Absolute database path."
    path::String
    "Serialize operations on this connection."
    lock::ReentrantLock
end

"""
    SQLRow

Hold an internal, detached SQL row with column-name property access. Its concrete
type is independent of query columns and cell values, including NULL. Convert
rows to the existing typed domain records before returning them from store APIs.
This avoids compiling new row/container specializations under live control locks
when a previously empty query first returns data.
"""
struct SQLRow
    "Copied column values, independent of the SQLite statement lifetime."
    columns::Dict{Symbol,Any}
end
Base.getproperty(row::SQLRow, name::Symbol) = getfield(row, :columns)[name]
Base.propertynames(row::SQLRow, private::Bool=false) = collect(keys(getfield(row, :columns)))

function sql_rows(db, sql, params=())
    statement = DBInterface.prepare(db, sql)
    try
        result = DBInterface.execute(statement, params)
        try
            rows = SQLRow[]
            for row in result
                columns = Dict{Symbol,Any}()
                for name in propertynames(row)
                    columns[name] = getproperty(row, name)
                end
                push!(rows, SQLRow(columns))
            end
            return rows
        finally
            DBInterface.close!(result)
        end
    finally
        DBInterface.close!(statement)
    end
end

function transaction(f, store::RuntimeStore)
    lock(store.lock) do
        sql_rows(store.db, "BEGIN IMMEDIATE")
        try
            result = f()
            sql_rows(store.db, "COMMIT")
            return result
        catch
            sql_rows(store.db, "ROLLBACK")
            rethrow()
        end
    end
end

const RUN_SCHEMA_VERSION = 5
const ACTIVE_RUN_STATES = (:reserved, :starting, :running, :stopping)
const RUN_TRANSITIONS = Dict(
    :reserved => (:starting, :failed, :stopped),
    :starting => (:running, :failed, :stopping),
    :running => (:stopping, :failed),
    :stopping => (:stopped, :failed),
    :failed => (),
    :stopped => (),
)

function RuntimeStore(path::AbstractString)
    absolute = abspath(path)
    islink(absolute) && throw(ArgumentError("runtime database cannot be a symbolic link"))
    parent = dirname(absolute)
    isdir(parent) || mkpath(parent; mode=0o700)
    existed = isfile(absolute)
    db = SQLite.DB(absolute)
    try
        existed || chmod(absolute, 0o600)
        sql_rows(db, "PRAGMA busy_timeout = 5000")
        sql_rows(db, "PRAGMA foreign_keys = ON")
        store = RuntimeStore(db, absolute, ReentrantLock())
        transaction(store) do
            version = only(sql_rows(db, "PRAGMA user_version")).user_version
            version in (1, 2, 3, 4) && throw(ArgumentError("older runtime schema requires lcm runtime migrate before startup"))
            version in (0, RUN_SCHEMA_VERSION) ||
                throw(ArgumentError("unsupported runtime database schema"))
            if version == 0
                isempty(sql_rows(db, "SELECT name FROM sqlite_master WHERE type='table'")) ||
                    throw(ArgumentError("refusing to initialize a non-runtime database"))
                sql_rows(db, """
                    CREATE TABLE runs (
                        run_id TEXT PRIMARY KEY,
                        owner TEXT NOT NULL,
                        application TEXT NOT NULL,
                        application_version TEXT NOT NULL,
                        state TEXT NOT NULL CHECK(state IN ('reserved','starting','running','stopping','failed','stopped')),
                        request_id TEXT NOT NULL,
                        created_at TEXT NOT NULL,
                        updated_at TEXT NOT NULL,
                        reason TEXT NOT NULL DEFAULT '',
                        UNIQUE(owner, request_id)
                    )
                """)
                sql_rows(db, "CREATE INDEX runs_owner_state ON runs(owner, state)")
                create_worker_schema!(db)
                create_lease_schema!(db)
                create_job_schema!(db)
                create_job_cancellation_schema!(db)
                sql_rows(db, "PRAGMA user_version = 5")
            else
                # Fail early if an unrelated database only claims our version.
                sql_rows(db, "SELECT run_id, owner, application, application_version, state, request_id, created_at, updated_at, reason FROM runs LIMIT 0")
                check_worker_schema!(db)
                check_lease_schema!(db)
                check_job_schema!(db)
                check_job_cancellation_schema!(db)
            end
        end
        sql_rows(db, "PRAGMA journal_mode = WAL")
        return store
    catch
        close(db)
        rethrow()
    end
end

Base.close(store::RuntimeStore) = lock(() -> close(store.db), store.lock)

record(row) = RunRecord(UUID(row.run_id), row.owner, row.application,
    VersionNumber(row.application_version), Symbol(row.state), UUID(row.request_id),
    DateTime(row.created_at), DateTime(row.updated_at), row.reason)

function owned_run(db, principal::Principal, id::UUID)
    rows = sql_rows(db, "SELECT * FROM runs WHERE run_id = ?", (string(id),))
    isempty(rows) && throw(AccessDenied(404, "Run not found"))
    row = only(rows)
    (row.owner == principal.id || principal.administrator) ||
        throw(AccessDenied(404, "Run not found"))
    return record(row)
end

"""
    get_run(store, principal, id) -> RunRecord

Read an owned run. Missing and foreign runs produce the same AccessDenied(404).
An explicitly configured administrator may inspect other owners' runs.
"""
get_run(store::RuntimeStore, principal::Principal, id::UUID) =
    lock(() -> owned_run(store.db, principal, id), store.lock)

"""
    list_runs(store, principal) -> Vector{RunRecord}

List only the caller's runs, or all runs for an explicit administrator.
"""
function list_runs(store::RuntimeStore, principal::Principal)
    lock(store.lock) do
        rows = principal.administrator ?
            sql_rows(store.db, "SELECT * FROM runs ORDER BY created_at, run_id") :
            sql_rows(store.db, "SELECT * FROM runs WHERE owner = ? ORDER BY created_at, run_id", (principal.id,))
        return record.(rows)
    end
end

"""
    reserve_run!(store, principal, application; limits=RunLimits(), request_id=uuid4())

Atomically reserve bounded capacity before starting a process. Repeating the same
owner/request ID returns its original run; reusing it for a different application
or version raises ArgumentError. Exceeding capacity raises CapacityUnavailable.
"""
function reserve_run!(store::RuntimeStore, principal::Principal,
        application::ApplicationDefinition; limits::RunLimits=RunLimits(), request_id::UUID=uuid4())
    transaction(store) do
        previous = sql_rows(store.db, "SELECT * FROM runs WHERE owner = ? AND request_id = ?",
            (principal.id, string(request_id)))
        if !isempty(previous)
            run = record(only(previous))
            run.application == application.id && run.version == application.version ||
                throw(ArgumentError("request ID already belongs to a different application"))
            return run
        end
        counts = only(sql_rows(store.db, """
            SELECT COUNT(*) AS total,
                COALESCE(SUM(CASE WHEN owner = ? THEN 1 ELSE 0 END), 0) AS owned
            FROM runs WHERE state IN ('reserved','starting','running','stopping')
        """, (principal.id,)))
        counts.total < limits.max_runs && counts.owned < limits.max_runs_per_owner ||
            throw(CapacityUnavailable())
        id, timestamp = uuid4(), string(now(UTC))
        sql_rows(store.db, """
            INSERT INTO runs (run_id,owner,application,application_version,state,
                request_id,created_at,updated_at) VALUES (?,?,?,?,'reserved',?,?,?)
        """, (string(id), principal.id, application.id, string(application.version),
            string(request_id), timestamp, timestamp))
        return owned_run(store.db, principal, id)
    end
end

"""
    transition_run!(store, principal, id, state; reason="") -> RunRecord

Apply an allowed lifecycle transition under ownership and a database transaction.
Repeated transitions to the current state are idempotent. A terminal run cannot
be revived; restarting requires a new run identity.
"""
function transition_run!(store::RuntimeStore, principal::Principal, id::UUID,
        state::Symbol; reason::AbstractString="")
    haskey(RUN_TRANSITIONS, state) || throw(ArgumentError("unknown run state"))
    length(reason) <= 240 && !occursin(r"[\r\n\x00]", reason) ||
        throw(ArgumentError("run reason must be one bounded, sanitized line"))
    transaction(store) do
        run = owned_run(store.db, principal, id)
        run.state == state && return run
        state in RUN_TRANSITIONS[run.state] ||
            throw(ArgumentError("invalid run lifecycle transition"))
        sql_rows(store.db, "UPDATE runs SET state=?,updated_at=?,reason=? WHERE run_id=?",
            (String(state), string(now(UTC)), String(reason), string(id)))
        return owned_run(store.db, principal, id)
    end
end
