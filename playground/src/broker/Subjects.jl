"JetStream stream containing durable calculation jobs."
const JOB_STREAM = "LCM_JOBS"

"JetStream stream containing durable terminal results."
const RESULT_STREAM = "LCM_RESULTS"

function job_subject(operation::AbstractString, priority::AbstractString="normal")
    priority in PRIORITIES || throw(ArgumentError("unsupported job priority: $priority"))
    return "lcm.jobs.v1.$priority.$(operation_subject_token(operation))"
end
result_subject(job_id::AbstractString) = "lcm.results.v1.$job_id"
event_subject(job_id::AbstractString) = "lcm.events.v1.$job_id"
cancel_subject(job_id::AbstractString) = "lcm.control.v1.cancel.$job_id"
heartbeat_subject(worker_id::AbstractString="*") =
    "lcm.workers.v1.heartbeat.$worker_id"
capability_subject(worker_id::AbstractString="*") =
    "lcm.workers.v1.capabilities.$worker_id"
