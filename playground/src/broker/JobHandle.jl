const JOB_STATES = (
    :editable,
    :dirty,
    :submitting,
    :queued,
    :unavailable,
    :running,
    :cancelling,
    :ready,
    :failed,
    :canceled,
    :expired,
)
const JOB_TERMINAL_STATES = (:ready, :failed, :canceled, :expired)

"""
    JobHandle

Hold one browser session's observable calculation lifecycle. A pending request
does not replace `last_successful_result`; the preceding valid result therefore
remains available while newer inputs are queued or running.
"""
mutable struct JobHandle
    "Current request, when one has been submitted."
    request::Union{Nothing,JobRequest}

    "Current job identifier."
    job_id::Observable{Union{Nothing,String}}

    "Lifecycle state."
    state::Observable{Symbol}

    "Completed fraction in the closed interval [0, 1]."
    progress::Observable{Float64}

    "Concise lifecycle stage."
    stage::Observable{String}

    "Messages associated with the current request."
    messages::Observable{Vector{ConsoleEntry}}

    "Most recent successful result, retained while newer work is pending."
    last_successful_result::Observable{Union{Nothing,JobResult}}

    "Input hash awaiting a terminal result."
    pending_input_hash::Observable{Union{Nothing,String}}

    "Whether selectors changed after the current request snapshot."
    dirty_inputs::Observable{Bool}

    "Monotonic selector-edit revision."
    input_revision::Int

    "Selector revision captured by the current request."
    submitted_revision::Int

    "Browser-safe current failure."
    failure::Observable{Union{Nothing,FailureInfo}}

    "Worker currently executing the job."
    worker::Observable{Union{Nothing,String}}

    "Highest lifecycle sequence already applied."
    last_event_sequence::Int

    "Subscription handles owned by the broker client."
    subscriptions::Vector{Any}

    "Lock protecting callbacks delivered by NATS tasks."
    lock::ReentrantLock
end

function JobHandle()
    return JobHandle(
        nothing,
        Observable{Union{Nothing,String}}(nothing),
        Observable(:editable),
        Observable(0.0),
        Observable("Not submitted"),
        Observable(ConsoleEntry[]),
        Observable{Union{Nothing,JobResult}}(nothing),
        Observable{Union{Nothing,String}}(nothing),
        Observable(false),
        0,
        0,
        Observable{Union{Nothing,FailureInfo}}(nothing),
        Observable{Union{Nothing,String}}(nothing),
        0,
        Any[],
        ReentrantLock()
    )
end

function set_job_state!(handle::JobHandle, state::Symbol, stage::AbstractString)
    state in JOB_STATES || throw(ArgumentError("unsupported job state: $state"))
    handle.state[] = state
    handle.stage[] = string(stage)
    return handle
end

function mark_dirty!(handle::JobHandle)
    lock(handle.lock) do
        handle.input_revision += 1
        handle.dirty_inputs[] = true
        handle.failure[] = nothing
        if handle.state[] in (:submitting, :queued, :unavailable, :running, :cancelling)
            handle.stage[] = "Inputs changed · current request still pending"
        else
            handle.pending_input_hash[] = nothing
            set_job_state!(handle, :dirty, "Inputs changed")
        end
    end
    return handle
end

function job_console_entry(event::JobEvent)
    channel = event.event_type == "JobLog" ? :stdout :
              event.event_type == "JobFailed" ? :error : :broker
    message = something(event.message, event.stage, event.event_type)
    return ConsoleEntry(channel, message)
end

function apply_event!(handle::JobHandle, event::JobEvent)
    lock(handle.lock) do
        event.event_sequence <= handle.last_event_sequence && return handle
        handle.last_event_sequence = event.event_sequence
        !isnothing(event.progress) && (handle.progress[] = event.progress)
        !isnothing(event.worker_id) && (handle.worker[] = event.worker_id)
        handle.messages[] = vcat(handle.messages[], [job_console_entry(event)])

        if event.event_type in ("JobClaimed", "JobStarted", "JobProgress", "JobLog")
            set_job_state!(
                handle,
                event.event_type == "JobClaimed" ? :queued : :running,
                something(event.stage, event.message, event.event_type)
            )
        elseif event.event_type == "JobSucceeded"
            handle.progress[] = 1.0
            set_job_state!(handle, :running, "Result stored · loading")
        elseif event.event_type == "JobFailed"
            set_job_state!(handle, :failed, something(event.stage, "Failed"))
        elseif event.event_type == "JobCanceled"
            set_job_state!(handle, :canceled, something(event.stage, "Canceled"))
        elseif event.event_type == "JobExpired"
            set_job_state!(handle, :expired, something(event.stage, "Expired"))
        end
    end
    return handle
end

function apply_result!(handle::JobHandle, result::JobResult)
    lock(handle.lock) do
        if isnothing(result.failure)
            handle.last_successful_result[] = result
            handle.failure[] = nothing
            handle.progress[] = 1.0
            handle.pending_input_hash[] = nothing
            if handle.dirty_inputs[]
                set_job_state!(handle, :dirty, "Result ready · selectors changed")
            else
                set_job_state!(handle, :ready,
                    result.cache_status == "hit" ? "Ready · cache hit" : "Ready")
            end
        else
            handle.failure[] = result.failure
            handle.pending_input_hash[] = nothing
            state = result.failure.category == "canceled" ? :canceled :
                    result.failure.category == "deadline_expired" ? :expired : :failed
            if handle.dirty_inputs[]
                set_job_state!(handle, :dirty,
                    "Previous request $(string(state)) · selectors changed")
            else
                set_job_state!(handle, state, result.failure.message)
            end
        end
    end
    return handle
end

function expire_handle!(handle::JobHandle, job_id::AbstractString)
    lock(handle.lock) do
        handle.job_id[] == job_id || return handle
        handle.state[] in JOB_TERMINAL_STATES && return handle
        failure = FailureInfo(
            "deadline_expired",
            "The calculation deadline expired before a result was received",
            string(uuid4()),
            false
        )
        handle.failure[] = failure
        handle.pending_input_hash[] = nothing
        set_job_state!(handle, :expired, failure.message)
    end
    return handle
end

function submission_failed!(handle::JobHandle, error)
    diagnostic_id = string(uuid4())
    failure = FailureInfo(
        "broker_unavailable",
        sprint(showerror, error),
        diagnostic_id,
        true
    )
    lock(handle.lock) do
        handle.failure[] = failure
        set_job_state!(handle, :unavailable, failure.message)
    end
    return handle
end
