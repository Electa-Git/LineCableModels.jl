"""
    JobRecord

Retain an authorized submission receipt, not proof of execution or a cached
scientific result. The exact wire request survives HTTP retries and coordinator
restart. JetStream remains the durable result authority.
"""
struct JobRecord
    "Exact server-authored request, deadline and prepared executor identity."
    job::Protocol.AssignedJob
    "Authenticated caller's idempotency key."
    request_id::UUID
    "Hash of the caller's lease, operation and normalized inputs."
    submission_hash::String
    "Queued, submitted, succeeded, failed, uncertain, canceled or revoked."
    state::Symbol
    "UTC receipt creation time."
    created_at::DateTime
    "UTC last bookkeeping change."
    updated_at::DateTime
end

const PENDING_JOB_SQL = "('queued','submitted')"

function create_job_schema!(db)
    sql_rows(db,"""
        CREATE TABLE jobs (
            job_id TEXT PRIMARY KEY,
            owner TEXT NOT NULL,
            run_id TEXT NOT NULL REFERENCES runs(run_id),
            lease_id TEXT NOT NULL REFERENCES leases(lease_id),
            request_id TEXT NOT NULL,
            submission_hash TEXT NOT NULL,
            job_json TEXT NOT NULL,
            state TEXT NOT NULL CHECK(state IN ('queued','submitted','succeeded','failed','uncertain','canceled','revoked')),
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL,
            UNIQUE(owner,request_id)
        )
    """)
    sql_rows(db,"CREATE INDEX jobs_run_created ON jobs(run_id,created_at,job_id)")
    sql_rows(db,"CREATE INDEX jobs_lease ON jobs(lease_id)")
    sql_rows(db,"CREATE UNIQUE INDEX jobs_pending_lease ON jobs(lease_id) WHERE state IN $PENDING_JOB_SQL")
end

check_job_schema!(db)=sql_rows(db,"""
    SELECT job_id,owner,run_id,lease_id,request_id,submission_hash,job_json,state,created_at,updated_at FROM jobs LIMIT 0
""")

function job_record(row)
    job=Protocol.decode_runtime_message(Protocol.AssignedJob,row.job_json)
    (job.request.job_id,job.fence.owner,job.fence.run_id,job.fence.lease_id)==
        (row.job_id,row.owner,row.run_id,row.lease_id) || throw(ArgumentError("Inconsistent persisted job identity"))
    digest=job_submission_hash(UUID(job.fence.lease_id),job.request.operation,job.request.parameters)
    digest==row.submission_hash || throw(ArgumentError("Inconsistent persisted submission inputs"))
    JobRecord(job,UUID(row.request_id),row.submission_hash,Symbol(row.state),DateTime(row.created_at),DateTime(row.updated_at))
end

function job_submission_hash(lease_id::UUID,operation::AbstractString,parameters::AbstractDict)
    normalized=Protocol.normalize_wire(parameters)
    ncodeunits(JSON3.write(normalized))<=65536 || throw(AccessDenied(400,"Scientific inputs exceed their byte limit"))
    Protocol.input_hash("runtime.submission",Dict("lease_id"=>string(lease_id),"operation"=>String(operation),"parameters"=>normalized))
end

function owned_job(db,principal::Principal,id::UUID)
    rows=sql_rows(db,"SELECT * FROM jobs WHERE job_id=?",(string(id),))
    isempty(rows) && throw(AccessDenied(404,"Job not found"))
    row=only(rows)
    row.owner==principal.id || principal.administrator || throw(AccessDenied(404,"Job not found"))
    job_record(row)
end

"""Read an owned job receipt; missing and foreign jobs have the same 404 response."""
get_job(store::RuntimeStore,principal::Principal,id::UUID)=lock(()->owned_job(store.db,principal,id),store.lock)

"""
    prior_job_submission(store, principal, lease_id, operation, parameters, request_id)

Resolve an identical retry before requiring new readiness. Return its original
receipt even after expiry; changed inputs with the same owner/request ID fail.
This read cannot publish, revive authority, or extend the saved deadline.
"""
function prior_job_submission(store::RuntimeStore,principal::Principal,lease_id::UUID,
        operation::AbstractString,parameters::AbstractDict,request_id::UUID)
    digest=job_submission_hash(lease_id,operation,parameters)
    lock(store.lock) do
        lease=owned_lease(store.db,principal,lease_id)
        rows=sql_rows(store.db,"SELECT * FROM jobs WHERE owner=? AND request_id=?",(lease.fence.owner,string(request_id)))
        isempty(rows) && return nothing
        previous=job_record(only(rows))
        previous.submission_hash==digest || throw(AccessDenied(409,"Job request identity was reused with different inputs"))
        previous
    end
end

"""
    reserve_job!(store, principal, job, request_id) -> JobRecord

Persist the exact server-authored submission before attempting broker delivery.
Allow one pending job per lease and at most 256 receipts in that lease's lifetime.
An identical retry returns its original job ID, deadline and prepared target.
The execution coordinator must establish live lease and fresh preparation before
calling this low-level transaction; a durable active row is not that authority.
"""
function reserve_job!(store::RuntimeStore,principal::Principal,job::Protocol.AssignedJob,request_id::UUID)
    Protocol.validate(job)
    fence=job.fence
    digest=job_submission_hash(UUID(fence.lease_id),job.request.operation,job.request.parameters)
    transaction(store) do
        lease=owned_lease(store.db,principal,UUID(fence.lease_id))
        previous=prior_job_submission(store,principal,UUID(fence.lease_id),job.request.operation,job.request.parameters,request_id)
        previous===nothing || return previous
        lease.fence==fence && lease.state==:active || throw(AccessDenied(409,"Job assignment is not active"))
        run=owned_run(store.db,principal,UUID(fence.run_id))
        run.state in (:reserved,:starting,:running) || throw(AccessDenied(409,"Job run is not active"))
        Protocol.parse_utc_timestamp(job.request.deadline)>now(UTC) || throw(AccessDenied(409,"Job deadline expired"))
        counts=only(sql_rows(store.db,"""
            SELECT COUNT(*) AS total,COALESCE(SUM(CASE WHEN state IN $PENDING_JOB_SQL THEN 1 ELSE 0 END),0) AS pending
            FROM jobs WHERE lease_id=?
        """,(fence.lease_id,)))
        counts.total<256 && counts.pending==0 || throw(AccessDenied(409,"Job admission requires an idle assignment with available history"))
        timestamp=string(now(UTC))
        sql_rows(store.db,"""
            INSERT INTO jobs(job_id,owner,run_id,lease_id,request_id,submission_hash,job_json,state,created_at,updated_at)
            VALUES (?,?,?,?,?,?,?,'queued',?,?)
        """,(job.request.job_id,fence.owner,fence.run_id,fence.lease_id,string(request_id),digest,
            Protocol.encode_message(job),timestamp,timestamp))
        owned_job(store.db,principal,UUID(job.request.job_id))
    end
end

"""
    transition_job!(store, principal, job_id, state) -> JobRecord

Advance submission bookkeeping monotonically. Confirmed terminal receipts are immutable.
An uncertain outcome may be resolved by later durable result evidence, never
returned to a queued/submitted state or used as permission to execute again.
Only the execution coordinator may call this after checking delivery/result
evidence; changing this row is not an acknowledgement to the worker or broker.
"""
function transition_job!(store::RuntimeStore,principal::Principal,id::UUID,state::Symbol)
    state in (:submitted,:succeeded,:failed,:uncertain,:canceled,:revoked) || throw(ArgumentError("Invalid job state"))
    transaction(store) do
        previous=owned_job(store.db,principal,id)
        previous.state==state && return previous
        previous.state in (:queued,:submitted) ||
            (previous.state==:uncertain && state in (:succeeded,:failed,:canceled)) ||
            throw(AccessDenied(409,"Job bookkeeping is already terminal"))
        sql_rows(store.db,"UPDATE jobs SET state=?,updated_at=? WHERE job_id=?",(String(state),string(now(UTC)),string(id)))
        owned_job(store.db,principal,id)
    end
end

"""Read at most 256 recent receipts for an authorized run, without retrieving scientific results."""
function list_jobs(store::RuntimeStore,principal::Principal,run_id::UUID)
    lock(store.lock) do
        owned_run(store.db,principal,run_id)
        job_record.(sql_rows(store.db,"SELECT * FROM jobs WHERE run_id=? ORDER BY created_at DESC,job_id DESC LIMIT 256",(string(run_id),)))
    end
end
