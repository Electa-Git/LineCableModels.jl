"""
    AbstractScientificDriver

Implement the physical-resource boundary for scientific assignments. Required
hooks verify installed profiles and launch-time policy, recover owned resources,
provide an owned ExecutorSupervisor, release exact fences, and close. These
operator-installed hooks cannot replace ScientificResources' lease, preparation,
admission, cancellation or result checks. No default driver bypasses preflight.
"""
abstract type AbstractScientificDriver end
const ScientificProcess = ExecutionCore.ExecutorSupervisor

@required AbstractScientificDriver begin
    installed_profiles(::AbstractScientificDriver)
    recover_owned!(::AbstractScientificDriver)
    verify_executor!(::AbstractScientificDriver, ::ProfileDefinition, ::AssignmentFence)
    executor_for!(::AbstractScientificDriver, ::ProfileDefinition, ::AssignmentFence)
    release_owned!(::AbstractScientificDriver, ::AssignmentFence)
    Base.close(::AbstractScientificDriver)
end

"""
    verify_executor!(driver, profile, fence)

Recheck the approved environment and physical resource policy before execution.
Return nothing only on success. Native source drift and unsupported required
limits must fail here, not be downgraded to a ready label.
"""
function verify_executor! end

"""
    executor_for!(driver, profile, fence) -> ExecutorSupervisor

Create or retrieve the exact fence's physically owned process supervisor.
The driver records partial acquisition so release_owned! can clean it even if
this hook throws. It must not prepare a model or evaluate a scientific operation.
"""
function executor_for! end

"""Retain one exact assignment's execution state, never shared between runs."""
mutable struct ScientificHandle
    "Complete lease authority, including both process incarnations."
    fence::AssignmentFence
    "Fresh resource identity for correlated diagnostics."
    id::UUID
    "Driver-owned process supervisor, absent before resource acquisition."
    supervisor::Union{Nothing,ScientificProcess}
    "At most one admitted task; completed results are retained only until the next request."
    task::Union{Nothing,Task}
    "Cancellation token for that admitted task."
    token::ExecutionCore.CancellationToken
    "Current explicit request identity."
    request_id::String
    "Current normalized passive request identity."
    request_hash::String
    "Preparation, operation, or inspection."
    request_kind::Symbol
    "Idle, starting, preparing, executing, failed or closing."
    phase::Symbol
    "Normalized successful preparation key, or nothing."
    preparation_key::Union{Nothing,String}
    "Process generation which actually produced that preparation."
    prepared_generation::Int
    "Local monotonic expiry; inspection does not extend it."
    prepared_until::Float64
    "Whether the last explicit preparation attempt failed."
    preparation_failed::Bool
    "Monotonic start of the current operation."
    started_at::Float64
    "Monotonic completion time, retained while idle."
    finished_at::Float64
    "Most recent bounded progress fraction."
    progress::Float64
    "Bounded output-line count; raw engine messages are not ordinary diagnostics."
    output_lines::Int
    "Fixed non-sensitive failure code."
    failure::Union{Nothing,String}
    "Whether starts have been permanently revoked for this handle."
    closing::Bool
    "Serialize physical release attempts without blocking other handles."
    cleanup_lock::ReentrantLock
end

function ScientificHandle(fence::AssignmentFence)
    ScientificHandle(fence, uuid4(), nothing, nothing, ExecutionCore.CancellationToken(),
        "", "", :prepare, :idle, nothing, 0, 0.0, false, 0.0, 0.0, 0.0, 0, nothing, false, ReentrantLock())
end
Base.show(io::IO, handle::ScientificHandle) = print(io, "ScientificHandle(", handle.id, ", <owned>)")

"""Retain bounded canceled job identities until their exact lease loses authority."""
struct JobCancellations
    "Full assignment identity, not only its UUID."
    fence::AssignmentFence
    "At most 256 canceled job UUIDs; no input values or results."
    ids::Set{String}
end

"""
    ScientificResources(driver)

Enforce one independently supervised scientific process per exact assignment.
Construction is passive. The agent binds its live ledger before recovery; only
explicit prepare_assigned! or execute_assigned! calls admit work. Different
assignments proceed independently. One pending request is allowed per assignment;
concurrent matching preparations share it, while different work is rejected.

The physical driver must enforce its approved resource policy. This owner does
not turn a native process into a container sandbox and never advertises a profile
for which the driver has not completed its installed-profile checks.
"""
mutable struct ScientificResources{D<:AbstractScientificDriver} <: AbstractAgentResources
    "Approved physical-resource adapter; not application-supplied callbacks."
    driver::D
    "Verified scientific definitions, not warm executors."
    profiles::ProfileRegistry
    "One agent incarnation's authoritative lease ledger."
    ledger::Union{Nothing,AgentLeaseLedger}
    "Occupied resources keyed by lease UUID, including unresolved cleanup."
    handles::Dict{String,ScientificHandle}
    "Cancellation tombstones independent of executor preparation or replacement."
    canceled_jobs::Dict{String,JobCancellations}
    "Whether startup recovery completed."
    recovered::Bool
    "Whether all new work is forbidden."
    closed::Bool
    "Serialize short admission/status updates, never numerical work."
    lock::ReentrantLock
    "Join whole-owner recovery and shutdown."
    lifecycle_lock::ReentrantLock
end

function ScientificResources(driver::AbstractScientificDriver)
    RequiredInterfaces.check_interface_implemented(AbstractScientificDriver, typeof(driver)) === true ||
        throw(ArgumentError("scientific driver does not implement required ownership hooks"))
    profiles = installed_profiles(driver)
    profiles isa ProfileRegistry && all(p -> p.kind == :scientific, values(profiles.definitions)) ||
        throw(ArgumentError("scientific driver must return verified scientific profiles"))
    return ScientificResources(driver, profiles, nothing, Dict{String,ScientificHandle}(),Dict{String,JobCancellations}(),
        false, false, ReentrantLock(), ReentrantLock())
end
installed_profiles(resources::ScientificResources) = resources.profiles

function bind_agent!(resources::ScientificResources, ledger::AgentLeaseLedger)
    lock(resources.lock) do
        resources.closed && throw(ArgumentError("scientific resources are closed"))
        resources.ledger === nothing || resources.ledger === ledger ||
            throw(ArgumentError("scientific resources cannot change agent incarnation"))
        for (id, profile) in resources.profiles.definitions
            get(ledger.profiles.definitions, id, nothing) === profile ||
                throw(ArgumentError("scientific resource definitions do not match lease authority"))
        end
        resources.ledger = ledger
    end
    return nothing
end

function recover_owned!(resources::ScientificResources)
    lock(resources.lifecycle_lock) do
        resources.closed && throw(ArgumentError("scientific resources are closed"))
        resources.ledger === nothing && throw(ArgumentError("scientific resources lack lease authority"))
        resources.recovered && return nothing
        recover_owned!(resources.driver) === nothing || throw(ArgumentError("scientific resource recovery failed"))
        resources.recovered = true
    end
    return nothing
end

function scientific_authority(resources::ScientificResources, fence::AssignmentFence)
    resources.closed && throw(AccessDenied(409, "Scientific resources are closed"))
    resources.recovered || throw(AccessDenied(409, "Scientific resources are not recovered"))
    ledger = resources.ledger
    ledger !== nothing && agent_lease_usable(ledger, fence) ||
        throw(AccessDenied(409, "Scientific assignment is no longer usable"))
    profile = get(resources.profiles.definitions, fence.profile_id, nothing)
    profile !== nothing && profile.fingerprint == fence.fingerprint &&
        string(profile.version) == fence.profile_version ||
        throw(AccessDenied(409, "Scientific profile does not match its assignment"))
    return profile
end

function prune_job_cancellations!(resources::ScientificResources)
    for (id,entry) in collect(resources.canceled_jobs)
        agent_lease_usable(resources.ledger,entry.fence) || delete!(resources.canceled_jobs,id)
    end
    return nothing
end

function check_job_cancellation(resources::ScientificResources,fence::AssignmentFence,job_id::String)
    lock(resources.lock) do
        scientific_authority(resources,fence)
        entry=get(resources.canceled_jobs,fence.lease_id,nothing)
        entry===nothing && return nothing
        entry.fence==fence || throw(AccessDenied(409,"Job cancellation fence differs"))
        job_id in entry.ids && throw(ExecutionCore.OperationCanceled())
        return nothing
    end
end

function invalidate_preparation!(handle::ScientificHandle)
    handle.preparation_key = nothing
    handle.prepared_generation = 0
    handle.prepared_until = 0.0
    return nothing
end

function preparation_current(resources::ScientificResources, handle::ScientificHandle)
    supervisor = handle.supervisor
    current = handle.preparation_key !== nothing && supervisor !== nothing &&
        supervisor.process !== nothing && process_running(supervisor.process) &&
        supervisor.generation == handle.prepared_generation &&
        resources.ledger.clock() < handle.prepared_until
    current || invalidate_preparation!(handle)
    return current
end

function retain_preparation!(resources, handle, status, key, sampled_at)
    status isa AbstractDict && get(status, "preparation_input_hash", nothing) == key &&
        get(status, "ready", nothing) isa Bool || throw(ArgumentError("invalid executor preparation status"))
    remaining = get(status, "remaining_seconds", nothing)
    remaining isa Real && !(remaining isa Bool) && isfinite(remaining) && 0 <= remaining <= 86400 ||
        throw(ArgumentError("invalid executor preparation lifetime"))
    lock(resources.lock) do
        if status["ready"] && remaining > 0 && !handle.closing
            handle.preparation_key = key
            handle.prepared_generation = handle.supervisor.generation
            # Subtract the whole request latency conservatively; receiving a
            # reply must not extend the child's original model expiry.
            handle.prepared_until = sampled_at + remaining
        else
            invalidate_preparation!(handle)
        end
    end
    return nothing
end

function query_preparation!(resources, handle, context, key)
    sampled_at = resources.ledger.clock()
    status = ExecutionCore.inspect_preparation!(handle.supervisor, context, key)
    retain_preparation!(resources, handle, status, key, sampled_at)
    return status
end

function scientific_context(resources, handle, deadline)
    progress = (fraction, _, _) -> lock(resources.lock) do
        fraction isa Real && isfinite(fraction) && 0 <= fraction <= 1 ||
            throw(ArgumentError("invalid executor progress"))
        handle.progress = Float64(fraction)
    end
    output = _ -> lock(resources.lock) do
        handle.output_lines = min(handle.output_lines + 1, 1_000_000)
    end
    return ExecutionCore.ExecutionContext(handle.request_id, handle.token, progress, output,
        String[], ExecutionCore.PreparedResourceCache(), deadline)
end

function scientific_failure(error)
    error isa ExecutionCore.OperationCanceled && return "canceled"
    error isa ExecutionCore.OperationDeadlineExpired && return "deadline"
    error isa AccessDenied && return "lease-unavailable"
    error isa ExecutionCore.PermanentOperationError && return "operation-rejected"
    return "executor-failed"
end

"""Retain one completed operation's value and child-reported schema independently of later inspections."""
struct ScientificOutput
    "Normalized scientific result."
    value::Dict{String,Any}
    "Schema version from the actual child operation registry."
    schema_version::String
    "Bounded scientific warnings returned by that operation."
    warnings::Vector{String}
end

function run_scientific_request!(resources, handle, profile, kind, parameters, operation, deadline)
    context = scientific_context(resources, handle, deadline)
    watching = Ref(true)
    guard = @async while watching[]
        if resources.closed || handle.closing || !agent_lease_usable(resources.ledger, handle.fence)
            ExecutionCore.cancel!(context.token)
        end
        sleep(0.025)
    end
    try
        scientific_authority(resources, handle.fence)
        verify_executor!(resources.driver, profile, handle.fence) === nothing ||
            throw(ArgumentError("scientific driver preflight failed"))
        ExecutionCore.check_canceled(context)
        if handle.supervisor === nothing
            acquired = executor_for!(resources.driver, profile, handle.fence)
            acquired isa ScientificProcess || throw(ArgumentError("scientific driver returned no owned process supervisor"))
            handle.supervisor = acquired
            handle.id = uuid4()
        end
        scientific_authority(resources, handle.fence)
        ExecutionCore.check_canceled(context)
        result = if kind == :prepare
            lock(resources.lock) do
                handle.phase = :preparing
                invalidate_preparation!(handle)
            end
            prepared = ExecutionCore.prepare_supervised!(handle.supervisor, context, parameters;
                timeout_seconds=profile.budget.prepare_seconds)
            key = get(prepared, "preparation_input_hash", nothing)
            key isa String && occursin(r"^[a-f0-9]{64}$", key) ||
                throw(ArgumentError("executor did not return preparation identity"))
            inspected = query_preparation!(resources, handle, context, key)
            get(inspected, "ready", false) && preparation_current(resources, handle) ||
                throw(ArgumentError("executor did not retain its prepared resources"))
            prepared
        else
            key = handle.preparation_key
            key !== nothing && preparation_current(resources, handle) ||
                throw(AccessDenied(409, "Explicit preparation is required"))
            inspected = query_preparation!(resources, handle, context, key)
            get(inspected, "ready", false) || throw(AccessDenied(409, "Scientific preparation expired"))
            if kind == :inspect
                inspected
            else
                lock(resources.lock) do
                    handle.phase = :executing
                end
                spec = ExecutionCore.OperationSpec(operation, identity, (_,_) -> nothing;
                    timeout_seconds=profile.budget.job_seconds, execution_mode=:supervised)
                schema = Ref("")
                value = ExecutionCore.execute_supervised!(handle.supervisor, spec, context, parameters; preparation_key=key,
                    on_result_schema=version->(schema[]=version))
                value isa Dict{String,Any} || throw(ArgumentError("Scientific operation must return an object"))
                inspected = query_preparation!(resources, handle, context, key)
                ScientificOutput(value,schema[],copy(context.warnings))
            end
        end
        ExecutionCore.check_canceled(context)
        # A queued old result cannot restore state after lease/source replacement.
        scientific_authority(resources, handle.fence)
        verify_executor!(resources.driver, profile, handle.fence) === nothing ||
            throw(ArgumentError("scientific driver post-execution verification failed"))
        lock(resources.lock) do
            handle.closing && throw(AccessDenied(409, "Scientific assignment is closing"))
            handle.phase = :idle
            handle.progress = 1.0
            handle.failure = nothing
        end
        return result
    catch error
        permanent = kind == :operation && error isa ExecutionCore.PermanentOperationError &&
            !ExecutionCore.iscanceled(context.token) && agent_lease_usable(resources.ledger, handle.fence)
        if permanent && error.category == "preparation_expired"
            lock(() -> invalidate_preparation!(handle), resources.lock)
        end
        if !permanent
            lock(resources.lock) do
                invalidate_preparation!(handle)
            end
            released = try
                handle.supervisor === nothing || ExecutionCore.stop_executor!(handle.supervisor)
                release_owned!(resources.driver, handle.fence) === true
            catch
                false
            end
            lock(resources.lock) do
                if released
                    handle.supervisor = nothing
                else
                    handle.closing = true
                end
            end
        end
        lock(resources.lock) do
            handle.phase = handle.closing ? :closing : :failed
            handle.failure = handle.closing ? "cleanup-unresolved" : scientific_failure(error)
            handle.preparation_failed = kind == :prepare
        end
    finally
        watching[] = false
        wait(guard)
        lock(resources.lock) do
            handle.finished_at = resources.ledger.clock()
        end
    end
    # Throw outside the catch block, after the original exception context has
    # unwound. Failed public tasks retain only the fixed, non-sensitive outcome.
    throw(AccessDenied(409, "Scientific request failed; inspect its owned status"))
end

function admit_scientific_request!(resources::ScientificResources, fence::AssignmentFence,
        kind::Symbol, parameters::AbstractDict; request_id::String, operation="", deadline=nothing, request_digest=nothing)
    kind in (:prepare, :operation, :inspect) || throw(ArgumentError("unsupported scientific request"))
    Protocol.validate(fence)
    Protocol.runtime_uuid(request_id)
    normalized = Protocol.normalize_wire(parameters)
    ncodeunits(JSON3.write(normalized)) <= 256 * 1024 ||
        throw(ArgumentError("scientific input exceeds its byte limit"))
    hash = something(request_digest, Protocol.input_hash(kind == :operation ? operation : "runtime.$kind", normalized))
    return lock(resources.lock) do
        profile = scientific_authority(resources, fence)
        kind == :operation && check_job_cancellation(resources,fence,request_id)
        kind == :operation && !(operation in profile.operations) &&
            throw(AccessDenied(403, "Operation is not permitted by this profile"))
        handle = get(resources.handles, fence.lease_id, nothing)
        if handle === nothing
            kind == :prepare || throw(AccessDenied(409, "Explicit preparation is required"))
            length(resources.handles) < resources.ledger.capacity ||
                throw(AccessDenied(409, "Scientific resource capacity is occupied"))
            handle = ScientificHandle(fence)
            resources.handles[fence.lease_id] = handle
        end
        handle.fence == fence && !handle.closing || throw(AccessDenied(409, "Scientific resource fence mismatch"))
        handle.request_id == request_id && (handle.request_hash != hash || handle.request_kind != kind) &&
            throw(AccessDenied(409, "Scientific request identity was reused"))
        if handle.task !== nothing && !istaskdone(handle.task)
            kind == :prepare && handle.request_kind == :prepare && handle.request_hash == hash && return handle.task
            handle.request_id == request_id && handle.request_hash == hash && return handle.task
            throw(AccessDenied(409, "Scientific assignment is busy"))
        end
        kind == :operation && handle.task !== nothing && handle.request_id == request_id && return handle.task
        if kind != :prepare
            preparation_current(resources, handle) || throw(AccessDenied(409, "Explicit preparation is required"))
        end
        handle.token = ExecutionCore.CancellationToken()
        handle.request_id, handle.request_hash, handle.request_kind = request_id, hash, kind
        handle.phase, handle.failure, handle.progress = :starting, nothing, 0.0
        handle.started_at = resources.ledger.clock()
        handle.finished_at = 0.0
        handle.preparation_failed = false
        kind == :prepare && invalidate_preparation!(handle)
        handle.output_lines = 0
        handle.task = Task(() -> run_scientific_request!(resources, handle, profile, kind,
            normalized, operation, deadline))
        schedule(handle.task)
        handle.task
    end
end

"""
    prepare_assigned!(resources, fence, parameters; request_id=uuid4()) -> Task

Admit explicit representative preparation under the current live lease. Matching
pending preparations share one task; different work is rejected, not queued
without a bound. No model is prepared merely because a page or lease exists.
"""
function prepare_assigned!(resources::ScientificResources, fence::AssignmentFence,
        parameters::AbstractDict; request_id=uuid4())
    admit_scientific_request!(resources, fence, :prepare, parameters; request_id=string(request_id))
end

"""
    execute_assigned!(resources, job) -> Task

Admit an allowlisted v2 job only after its exact process has prepared successfully.
Check its absolute deadline and lease continuously, and recheck authority before
returning its result. Durable result-before-ack remains the broker owner's duty;
this function does not claim exactly-once effects or publish a result.
The task returns ScientificOutput, keeping value, schema and warnings together.
"""
function execute_assigned!(resources::ScientificResources, job::Protocol.AssignedJob)
    Protocol.validate(job)
    deadline = Protocol.parse_utc_timestamp(job.request.deadline)
    deadline > now(UTC) || throw(AccessDenied(409, "Scientific request deadline expired"))
    job.request.engine_constraint === nothing || throw(AccessDenied(400,
        "Assigned jobs use the approved profile version and fingerprint, not a legacy engine constraint"))
    request_digest = Protocol.input_hash("runtime.assigned_request", JSON3.read(JSON3.write(job)))
    lock(resources.lock) do
        check_job_cancellation(resources,job.fence,job.request.job_id)
        prepared_execution(resources,job.fence) == job.execution ||
            throw(AccessDenied(409,"Prepared executor or model changed before job admission"))
        admit_scientific_request!(resources, job.fence, :operation, job.request.parameters;
            request_id=job.request.job_id, operation=job.request.operation, deadline, request_digest)
    end
end

"""
    prepared_execution(resources, fence) -> Protocol.PreparedExecution

Read the exact currently retained child/model identity under live lease authority.
This local snapshot is not a substitute for a fresh remote status inspection.
Throw AccessDenied when no usable preparation exists; never start a child here.
"""
function prepared_execution(resources::ScientificResources,fence::AssignmentFence)
    lock(resources.lock) do
        scientific_authority(resources,fence)
        handle=get(resources.handles,fence.lease_id,nothing)
        handle!==nothing && handle.fence==fence && !handle.closing && preparation_current(resources,handle) ||
            throw(AccessDenied(409,"Explicit preparation is required"))
        Protocol.PreparedExecution(string(handle.id),handle.prepared_generation,handle.preparation_key)
    end
end

"""
    refresh_preparation!(resources, fence) -> Task

Explicitly inspect retained preparation in the current child without recomputing
it or extending its idle TTL. Like a job, inspection cannot race another admitted
request or cross an expired fence.
"""
function refresh_preparation!(resources::ScientificResources, fence::AssignmentFence)
    admit_scientific_request!(resources, fence, :inspect, Dict{String,Any}(); request_id=string(uuid4()))
end

"""
    cancel_assigned!(resources, fence, request_id) -> Bool

Cancel only the exact currently pending request. A delayed cancellation cannot
stop its successor. Cancellation does not release the lease or replay any input.
"""
function cancel_assigned!(resources::ScientificResources, fence::AssignmentFence, request_id::String)
    lock(resources.lock) do
        scientific_authority(resources, fence)
        handle = get(resources.handles, fence.lease_id, nothing)
        handle !== nothing && handle.fence == fence && !handle.closing &&
            handle.request_id == request_id && handle.task !== nothing && !istaskdone(handle.task) || return false
        ExecutionCore.cancel!(handle.token)
        return true
    end
end

"""
    cancel_job!(resources, fence, job_id) -> Bool

Record cancellation before a queued job can be admitted, and cancel that exact
job if it is already running. Return true only when a current operation received
the interrupt token; false still means future admission is fenced. A previously
stored result wins over a late cancellation. Preparation and unrelated jobs are
untouched. Tombstones survive executor cleanup/repreparation under the same live
lease; lease loss or worker restart invalidates that job's authority instead.
"""
function cancel_job!(resources::ScientificResources,fence::AssignmentFence,job_id::String)
    Protocol.runtime_uuid(job_id)
    lock(resources.lock) do
        scientific_authority(resources,fence)
        prune_job_cancellations!(resources)
        entry=get(resources.canceled_jobs,fence.lease_id,nothing)
        if entry===nothing
            length(resources.canceled_jobs)<resources.ledger.capacity || throw(AccessDenied(409,"Cancellation capacity is occupied"))
            entry=JobCancellations(fence,Set{String}())
            resources.canceled_jobs[fence.lease_id]=entry
        end
        entry.fence==fence || throw(AccessDenied(409,"Job cancellation fence differs"))
        job_id in entry.ids || length(entry.ids)<256 || throw(AccessDenied(409,"Job cancellation history requires a new assignment"))
        push!(entry.ids,job_id)
        handle=get(resources.handles,fence.lease_id,nothing)
        handle!==nothing && handle.request_kind==:operation && handle.request_id==job_id || return false
        cancel_assigned!(resources,fence,job_id)
    end
end

"""
    scientific_status(resources, fence)

Return a local, non-sensitive status snapshot. Readiness requires live authority,
a live matching process generation and unexpired retained evidence. Use
refresh_preparation! before publishing a fresh remote preparation report; a
snapshot alone cannot prove the child still retains an evictable model.
"""
function scientific_status(resources::ScientificResources, fence::AssignmentFence)
    lock(resources.lock) do
        scientific_authority(resources, fence)
        handle = get(resources.handles, fence.lease_id, nothing)
        handle === nothing && return (phase=:idle, preparation=:cold, executor_id=nothing,
            generation=0, request_id=nothing, preparation_key=nothing, progress=0.0,
            elapsed_seconds=0.0, output_lines=0, failure=nothing)
        handle.fence == fence || throw(AccessDenied(409, "Scientific resource fence mismatch"))
        ready = preparation_current(resources, handle)
        (phase=handle.phase, preparation=handle.request_kind == :prepare && handle.phase in (:starting, :preparing) ? :preparing : ready ? :ready :
            handle.preparation_failed ? :failed : :cold, executor_id=string(handle.id),
            generation=handle.supervisor === nothing ? 0 : handle.supervisor.generation,
            request_id=handle.request_id, preparation_key=handle.preparation_key,
            progress=handle.progress, elapsed_seconds=max(0.0,
                (handle.task !== nothing && !istaskdone(handle.task) ? resources.ledger.clock() : handle.finished_at)-handle.started_at),
            output_lines=handle.output_lines, failure=handle.failure)
    end
end

function release_owned!(resources::ScientificResources, fence::AssignmentFence)
    handle = lock(resources.lock) do
        prune_job_cancellations!(resources)
        found = get(resources.handles, fence.lease_id, nothing)
        found === nothing && return nothing
        found.fence == fence || throw(AccessDenied(409, "Scientific cleanup fence mismatch"))
        found.closing = true
        found.phase = :closing
        invalidate_preparation!(found)
        ExecutionCore.cancel!(found.token)
        found
    end
    # Include partial acquisitions not yet represented by a completed supervisor.
    handle === nothing && return release_owned!(resources.driver, fence) === true
    return lock(handle.cleanup_lock) do
        handle.supervisor === nothing || ExecutionCore.stop_executor!(handle.supervisor)
        if handle.task !== nothing
            timedwait(() -> istaskdone(handle.task), 10; pollint=0.025)
            istaskdone(handle.task) || return false
        end
        release_owned!(resources.driver, fence) === true || return false
        lock(resources.lock) do
            get(resources.handles, fence.lease_id, nothing) === handle && delete!(resources.handles, fence.lease_id)
        end
        true
    end
end

function Base.close(resources::ScientificResources)
    lock(resources.lifecycle_lock) do
        lock(resources.lock) do
            resources.closed = true
            foreach(handle -> ExecutionCore.cancel!(handle.token), values(resources.handles))
        end
        fences = lock(() -> [h.fence for h in values(resources.handles)], resources.lock)
        tasks = [@async release_owned!(resources, fence) for fence in fences]
        successes = [try fetch(task) === true catch; false end for task in tasks]
        all(successes) || throw(ArgumentError("scientific resource cleanup remains unresolved"))
        close(resources.driver)
        lock(()->empty!(resources.canceled_jobs),resources.lock)
    end
    return nothing
end
