"""
    control_capabilities(service)

Report configured control services, not worker readiness or isolation evidence.
Both public capability discovery and private inventory use this declaration.
Submission still requires the job owner's live lease and fresh prepared target.
Private terminal discovery reports the configured relay only. Attachment still
requires owner authorization, a live lease and a verified container profile.
"""
control_capabilities(::Nothing) =
    (worker_control=false, preparation_control=false, assigned_execution=false, private_terminal=false)
control_capabilities(service::ControlService) = (
    worker_control=!service.closed,
    preparation_control=!service.closed && !service.science.closed,
    assigned_execution=!service.closed && !service.science.closed && !service.jobs.closed,
    private_terminal=!service.closed && !service.terminals.closed,
)

function control_snapshot(::Nothing, principal::Principal)
    return (; schema_version=1, enabled=false, control_capabilities(nothing)..., broker="disabled", administrator=principal.administrator,
        workers=NamedTuple[], profiles=NamedTuple[], provisioned=NamedTuple[])
end

"""
    control_snapshot(service, principal)

Return shared inventory dimensions for the panel and owned widgets. Environment
paths, connection URLs, credentials and other owners' assignment identities are
excluded. Profile advertisement is explicitly not preparation evidence.
"""
function control_snapshot(service::ControlService, principal::Principal)
    occupied = lock(service.inventory.store.lock) do
        sql_rows(service.inventory.store.db,
            "SELECT worker_id,COUNT(*) AS occupied FROM leases WHERE state IN $OCCUPIED_LEASE_SQL GROUP BY worker_id")
    end
    counts = Dict{String,Int}(row.worker_id => row.occupied for row in occupied)
    workers = [merge(item, (occupied=get(counts, item.registration.worker_id, 0),
        preparation="unknown")) for item in worker_inventory(service.inventory, principal)]
    profiles = [(id=p.id, version=string(p.version), kind=String(p.kind), isolation=String(p.isolation),
        fingerprint=p.fingerprint, operations=p.operations, preparation=p.preparation)
        for p in sort!(collect(values(service.config.profiles.definitions)); by=p -> p.id)]
    provisioned = principal.administrator ?
        [(worker_id=w.worker_id, credential_ref=w.credential_ref, profiles=w.profiles, capacity=w.capacity)
            for w in sort!(collect(values(service.config.workers)); by=w -> w.worker_id)] : NamedTuple[]
    return (; schema_version=1, enabled=!service.closed, control_capabilities(service)...,
        broker=String(lock(() -> service.state, service.lock)), administrator=principal.administrator,
        workers=workers, profiles=profiles, provisioned=provisioned)
end

function assignment_payload(service::ControlService, principal::Principal, lease::LeaseRecord)
    fence = lease.fence
    usable = assignment_usable(service.coordinator,principal,UUID(fence.lease_id))
    preparation = usable ? try remote_scientific_status(service.science,principal,UUID(fence.lease_id)).preparation catch; "unknown" end : "unknown"
    return (id=fence.lease_id, run_id=fence.run_id, role=fence.role, profile=fence.profile_id,
        worker_id=fence.worker_id, worker_boot=fence.worker_boot, generation=fence.generation,
        placement=String(lease.placement), state=String(lease.state), revision=lease.revision,
        usable, preparation)
end

function requested_fields(data, names)
    Set(keys(data)) == Set(names) || throw(AccessDenied(400, "Unexpected or missing control fields"))
    return data
end

function job_payload(service::ControlService,principal::Principal,receipt::JobRecord)
    job=receipt.job
    cancellation=job_cancellation(service.inventory.store,principal,UUID(job.request.job_id))
    return (id=job.request.job_id,request_id=string(receipt.request_id),run_id=job.fence.run_id,
        lease_id=job.fence.lease_id,role=job.fence.role,operation=job.request.operation,
        input_hash=job.request.input_hash,worker_id=job.fence.worker_id,worker_boot=job.fence.worker_boot,
        generation=job.fence.generation,executor_id=job.execution.executor_id,
        executor_generation=job.execution.executor_generation,preparation_key=job.execution.preparation_key,
        submitted_at=job.request.submitted_at,deadline=job.request.deadline,state=String(receipt.state),
        current_assignment=current_job_assignment(service.jobs,receipt),channel=String(service.jobs.state),
        cancel_requested=cancellation!==nothing,
        cancel_acknowledged=cancellation!==nothing && cancellation.acknowledged)
end

function requested_placement(data)
    data isa AbstractDict || throw(AccessDenied(400, "Expected a placement object"))
    mode = get(data, "mode", nothing)
    if mode == "automatic"
        requested_fields(data, ("mode",))
        return AutomaticPlacement()
    elseif mode == "pinned"
        requested_fields(data, ("mode", "worker_id"))
        data["worker_id"] isa AbstractString || throw(AccessDenied(400, "Invalid worker identity"))
        return PinnedPlacement(data["worker_id"])
    elseif mode == "dedicated"
        requested_fields(data, ("mode", "worker_id"))
        id = data["worker_id"]
        id === nothing || id isa AbstractString || throw(AccessDenied(400, "Invalid worker identity"))
        return DedicatedPlacement(id)
    end
    throw(AccessDenied(400, "Unsupported placement mode"))
end

function control_request(service::Union{Nothing,ControlService}, supervisor::UIHostSupervisor,
        principal::Principal, stream, path::String)
    method = stream.message.method
    if path == "/runtime/api/control" && method == "GET"
        gateway_response(stream, 200, JSON3.write(control_snapshot(service, principal)))
        return true
    end
    is_control = startswith(path, "/runtime/api/control/") || path == "/runtime/api/workers" ||
        startswith(path, "/runtime/api/workers/") || startswith(path, "/runtime/api/assignments/") ||
        startswith(path,"/runtime/api/jobs/") ||
        occursin(r"^/runtime/api/runs/[a-f0-9-]{36}/(?:assignments|jobs)$", path)
    is_control || return false
    service === nothing && throw(AccessDenied(503, "Worker control is not configured"))
    store = service.inventory.store
    submission_route=match(r"^/runtime/api/assignments/([a-f0-9-]{36})/jobs$",path)
    if submission_route!==nothing
        id=requested_uuid(submission_route[1])
        get_assignment(store,principal,id) # ownership before parsing private input
        method=="POST" || throw(AccessDenied(405,"Method is not supported"))
        data=requested_fields(read_gateway_json(stream),("operation","parameters","request_id"))
        data["operation"] isa AbstractString && data["parameters"] isa AbstractDict ||
            throw(AccessDenied(400,"Expected registered operation and input object"))
        receipt=submit_job!(service.jobs,principal,id,data["operation"],data["parameters"];
            request_id=requested_uuid(data["request_id"]))
        gateway_response(stream,202,JSON3.write(job_payload(service,principal,receipt)))
        return true
    end
    jobs_route=match(r"^/runtime/api/runs/([a-f0-9-]{36})/jobs$",path)
    if jobs_route!==nothing
        id=requested_uuid(jobs_route[1])
        get_run(store,principal,id)
        method=="GET" || throw(AccessDenied(405,"Method is not supported"))
        gateway_response(stream,200,JSON3.write([job_payload(service,principal,r) for r in list_jobs(store,principal,id)]))
        return true
    end
    job_route=match(r"^/runtime/api/jobs/([a-f0-9-]{36})(?:/(result|cancel|artifact))?$",path)
    if job_route!==nothing
        id=requested_uuid(job_route[1])
        receipt=get_job(store,principal,id)
        action=job_route[2]
        if action=="cancel" && method=="POST"
            data=requested_fields(read_gateway_json(stream),("request_id",))
            receipt=cancel_job!(service.jobs,principal,id;request_id=requested_uuid(data["request_id"]))
            gateway_response(stream,202,JSON3.write(job_payload(service,principal,receipt)))
        elseif action=="artifact" && method in ("GET","HEAD")
            bytes=owned_job_artifact(service.jobs,principal,id)
            gateway_response(stream,200,bytes;headers=["Content-Disposition"=>"attachment; filename=\"result.json\""])
        elseif action=="result" && method=="GET"
            result=owned_job_result(service.jobs,principal,id)
            gateway_response(stream,200,JSON3.write((schema_version=1,
                job=job_payload(service,principal,get_job(store,principal,id)),result)))
        elseif action===nothing && method=="GET"
            gateway_response(stream,200,JSON3.write(job_payload(service,principal,receipt)))
        else
            throw(AccessDenied(405,"Method is not supported"))
        end
        return true
    end
    if path == "/runtime/api/control/events" && method == "GET"
        query = URIs.queryparampairs(URIs.URI(stream.message.target))
        length(query) <= 2 && allunique(first.(query)) &&
            all(p -> first(p) in ("after", "epoch"), query) ||
            throw(AccessDenied(400, "Invalid event cursor"))
        parameters = Dict(query)
        after = tryparse(Int, get(parameters, "after", "0"))
        after !== nothing && after >= 0 || throw(AccessDenied(400, "Invalid event cursor"))
        epoch = get(parameters, "epoch", nothing)
        gateway_response(stream, 200, JSON3.write(control_events(service.events, principal; after, epoch)))
        return true
    elseif path == "/runtime/api/workers"
        if method == "GET"
            gateway_response(stream, 200, JSON3.write(control_snapshot(service, principal).workers))
        elseif method == "POST"
            require_administrator(principal)
            data = requested_fields(read_gateway_json(stream), ("worker_id", "request_id"))
            data["worker_id"] isa AbstractString || throw(AccessDenied(400, "Invalid worker identity"))
            trust = get(service.config.workers, data["worker_id"], nothing)
            trust === nothing && throw(AccessDenied(400, "Worker identity was not provisioned by the operator"))
            worker = enroll_worker!(store, principal, trust; request_id=requested_uuid(data["request_id"]))
            record_event!(service.events, :worker_enrolled; worker_id=worker.worker_id)
            gateway_response(stream, 200, JSON3.write(worker_payload(worker)))
        else
            throw(AccessDenied(405, "Method is not supported"))
        end
        return true
    end
    worker_route = match(r"^/runtime/api/workers/([a-z0-9][a-z0-9_-]{0,63})$", path)
    if worker_route !== nothing && method == "PATCH"
        require_administrator(principal)
        data = requested_fields(read_gateway_json(stream), ("state", "expected_revision", "request_id"))
        data["state"] in ("approved", "draining", "disabled") ||
            throw(AccessDenied(400, "Unsupported registration state"))
        revision = data["expected_revision"]
        revision isa Integer && !(revision isa Bool) || throw(AccessDenied(400, "Invalid registration revision"))
        haskey(service.config.workers, worker_route[1]) ||
            throw(AccessDenied(409, "Registration no longer has a provisioned trust binding"))
        worker = set_registration_state!(store, principal, worker_route[1], Symbol(data["state"]);
            expected_revision=revision, request_id=requested_uuid(data["request_id"]))
        record_event!(service.events, :worker_registration_changed; worker_id=worker.worker_id)
        gateway_response(stream, 200, JSON3.write(worker_payload(worker)))
        return true
    end
    assignments = match(r"^/runtime/api/runs/([a-f0-9-]{36})/assignments$", path)
    if assignments !== nothing
        run_id = requested_uuid(assignments[1])
        get_run(store, principal, run_id) # authorize before inputs or availability
        if method == "GET"
            leases = list_assignments(store, principal; run_id)
            gateway_response(stream, 200, JSON3.write([assignment_payload(service, principal, l) for l in leases]))
        elseif method == "POST"
            data = requested_fields(read_gateway_json(stream), ("role", "profile", "placement", "request_id"))
            all(k -> data[k] isa AbstractString, ("role", "profile")) ||
                throw(AccessDenied(400, "Invalid role or profile"))
            payload = lock(service.coordinator.lock) do
                lease = reserve_assignment!(service.coordinator.assignments, principal, run_id,
                    data["role"], data["profile"]; placement=requested_placement(data["placement"]),
                    request_id=requested_uuid(data["request_id"]))
                if lease.state == :reserving && lease.revision == 0
                    grant_assignment!(service.coordinator, principal, UUID(lease.fence.lease_id))
                    record_event!(service.events, :assignment_reserved; fence=lease.fence)
                end
                current = get_assignment(store, principal, UUID(lease.fence.lease_id))
                assignment_payload(service, principal, current)
            end
            gateway_response(stream, 202, JSON3.write(payload))
        else
            throw(AccessDenied(405, "Method is not supported"))
        end
        return true
    end
    lease_route = match(r"^/runtime/api/assignments/([a-f0-9-]{36})$", path)
    science_route = match(r"^/runtime/api/assignments/([a-f0-9-]{36})/science$",path)
    if science_route !== nothing
        id = requested_uuid(science_route[1])
        get_assignment(store,principal,id) # authorize before reading private inputs
        if method == "GET"
            try
                request_scientific!(service.science,principal,id,"status")
            catch error
                error isa BrokerUnavailable || rethrow()
            end
        elseif method == "POST"
            data = read_gateway_json(stream)
            action = get(data,"action",nothing)
            if action == "prepare"
                requested_fields(data,("action","parameters","request_id"))
                data["parameters"] isa AbstractDict || throw(AccessDenied(400,"Expected preparation input object"))
                request_scientific!(service.science,principal,id,action;parameters=data["parameters"],request_id=requested_uuid(data["request_id"]))
            elseif action == "cancel"
                requested_fields(data,("action","target_id","request_id"))
                request_scientific!(service.science,principal,id,action;target_id=string(requested_uuid(data["target_id"])),request_id=requested_uuid(data["request_id"]))
            else
                throw(AccessDenied(400,"Unsupported scientific control action"))
            end
        else
            throw(AccessDenied(405,"Method is not supported"))
        end
        gateway_response(stream,method=="GET" ? 200 : 202,JSON3.write(remote_scientific_status(service.science,principal,id)))
        return true
    end
    if lease_route !== nothing
        id = requested_uuid(lease_route[1])
        lease = get_assignment(store, principal, id)
        if method == "DELETE"
            data = requested_fields(read_gateway_json(stream), ("request_id",))
            release_assignment!(service.coordinator, principal, id; request_id=requested_uuid(data["request_id"]))
            record_event!(service.events, :assignment_releasing; fence=lease.fence)
        elseif method != "GET"
            throw(AccessDenied(405, "Method is not supported"))
        end
        gateway_response(stream, 200, JSON3.write(assignment_payload(service, principal,
            get_assignment(store, principal, id))))
        return true
    end
    throw(AccessDenied(404, "Worker control route not found"))
end
