"""Retain cancellation intent separately from the eventual durable job outcome."""
struct JobCancellation
    "Exact owned job identity."
    job_id::UUID
    "Original idempotent cancellation request."
    request_id::UUID
    "Whether the exact worker incarnation acknowledged its cancellation tombstone."
    acknowledged::Bool
    "UTC intent creation time, not an execution completion time."
    created_at::DateTime
end

function create_job_cancellation_schema!(db)
    sql_rows(db,"""
        CREATE TABLE job_cancellations (
            job_id TEXT PRIMARY KEY REFERENCES jobs(job_id),
            owner TEXT NOT NULL,
            request_id TEXT NOT NULL,
            acknowledged INTEGER NOT NULL DEFAULT 0 CHECK(acknowledged IN (0,1)),
            created_at TEXT NOT NULL,
            UNIQUE(owner,request_id)
        )
    """)
end
check_job_cancellation_schema!(db)=sql_rows(db,"SELECT job_id,owner,request_id,acknowledged,created_at FROM job_cancellations LIMIT 0")
cancellation_record(row)=JobCancellation(UUID(row.job_id),UUID(row.request_id),row.acknowledged==1,DateTime(row.created_at))

"""Read owned cancellation intent; acknowledgement is not proof of job termination."""
function job_cancellation(store::RuntimeStore,principal::Principal,id::UUID)
    lock(store.lock) do
        owned_job(store.db,principal,id)
        rows=sql_rows(store.db,"SELECT * FROM job_cancellations WHERE job_id=?",(string(id),))
        isempty(rows) ? nothing : cancellation_record(only(rows))
    end
end

"""
    request_job_cancellation!(store, principal, job_id, request_id)

Persist an owned cancellation intent before sending control traffic. Repeated
cancellation of the same job returns its original intent. A reused request ID
cannot target another job. Return nothing for an already terminal receipt;
neither a canceled intent nor an acknowledged tombstone fabricates a result.
"""
function request_job_cancellation!(store::RuntimeStore,principal::Principal,id::UUID,request_id::UUID)
    transaction(store) do
        receipt=owned_job(store.db,principal,id)
        reused=sql_rows(store.db,"SELECT job_id FROM job_cancellations WHERE owner=? AND request_id=?",
            (receipt.job.fence.owner,string(request_id)))
        isempty(reused) || only(reused).job_id==string(id) || throw(AccessDenied(409,"Cancellation identity belongs to another job"))
        previous=job_cancellation(store,principal,id)
        previous===nothing || return previous
        receipt.state in (:queued,:submitted) || return nothing
        sql_rows(store.db,"INSERT INTO job_cancellations(job_id,owner,request_id,created_at) VALUES (?,?,?,?)",
            (string(id),receipt.job.fence.owner,string(request_id),string(now(UTC))))
        job_cancellation(store,principal,id)
    end
end

function acknowledge_job_cancellation!(store::RuntimeStore,principal::Principal,id::UUID,request_id::UUID)
    transaction(store) do
        intent=job_cancellation(store,principal,id)
        intent!==nothing && intent.request_id==request_id || throw(AccessDenied(409,"Job cancellation intent differs"))
        sql_rows(store.db,"UPDATE job_cancellations SET acknowledged=1 WHERE job_id=?",(string(id),))
        job_cancellation(store,principal,id)
    end
end
