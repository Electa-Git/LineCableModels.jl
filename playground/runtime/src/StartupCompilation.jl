# Compile the installed orchestration specializations before announcing live
# control availability. These are compiler requests, not invocations: no lease,
# job, prepared model or child process is created by this startup step. In
# particular, scientific packages remain exclusively in the executor process.
function compile_runtime_paths(service::ControlService)
    compile_terminal_wire_paths()
    coordinator, jobs, science = typeof(service.coordinator), typeof(service.jobs), typeof(service.science)
    assignments = typeof(service.coordinator.assignments)
    # The first reservation must not compile its placement branches while
    # holding the inventory lock: that prevents fresh reports being accepted.
    # Compile every supported placement without allocating a synthetic lease.
    for placement in (AutomaticPlacement, PinnedPlacement, DedicatedPlacement)
        keywords = NamedTuple{(:placement, :request_id),Tuple{placement,UUID}}
        precompile(Core.kwcall, (keywords, typeof(reserve_assignment!), assignments,
            Principal, UUID, String, String))
        precompile(allocation_candidate, (assignments, Vector{SQLRow}, RunRecord,
            ProfileDefinition, placement))
    end
    for control in (grant_assignment!, renew_assignment!, release_assignment!)
        precompile(control, (coordinator, Principal, UUID))
        precompile(Core.kwcall, (NamedTuple{(:request_id,),Tuple{UUID}},
            typeof(control), coordinator, Principal, UUID))
    end
    precompile(accept_lease_ack!, (coordinator, String, Protocol.LeaseAcknowledgement))
    precompile(accept_control_record!, (typeof(service), ControlEnvelope{Protocol.LeaseAcknowledgement}))
    precompile(request_scientific!, (science, Principal, UUID, String))
    precompile(submit_job!, (jobs, Principal, UUID, String, Dict{String,Any}))
    precompile(Core.kwcall,(NamedTuple{(:request_id,),Tuple{UUID}},typeof(submit_job!),
        jobs,Principal,UUID,String,Dict{String,Any}))
    precompile(reconcile_job!, (jobs, UUID))
    precompile(owned_job_result, (jobs, Principal, UUID))
    terminal = typeof(service.terminals)
    for session in (Nothing,String), writer in (Nothing,String)
        keywords = NamedTuple{(:session_id,:writer_id,:input_sequence,:after,:columns,:rows,:bytes,:request_id,:retry),
            Tuple{session,writer,Int,Int,Int,Int,Vector{UInt8},UUID,Bool}}
        precompile(Core.kwcall,(keywords,typeof(request_terminal!),terminal,Principal,UUID,String))
    end
    precompile(accept_terminal_report!, (terminal,String,Protocol.TerminalReport))
    return nothing
end

function compile_runtime_paths(agent::AgentService)
    compile_terminal_wire_paths()
    precompile(receive_agent_record!, (typeof(agent), Protocol.LeaseControl))
    precompile(close_agent_service!, (typeof(agent),))
    compile_scientific_paths(scientific_resources(agent.resources))
    compile_terminal_paths(terminal_resources(agent.resources))
    if agent.terminals !== nothing
        precompile(receive_terminal_command!, (typeof(agent.terminals), Protocol.TerminalCommand))
        precompile(perform_terminal_action, (typeof(agent.terminals), Protocol.TerminalCommand))
    end
    if agent.jobs !== nothing
        jobs = typeof(agent.jobs)
        connection = BrokerJobs{WorkerIdentity}
        precompile(accept_assigned_delivery!, (jobs,AssignedDelivery))
        precompile(admit_delivered_job!, (jobs,AgentJobFlight))
        precompile(execute_delivered_job!, (jobs,AgentJobFlight,connection))
        precompile(run_delivered_job!, (jobs,AgentJobFlight,connection))
        precompile(tick_jobs!, (jobs,))
    end
    return nothing
end

function compile_terminal_wire_paths()
    for record in (Protocol.TerminalCommand,Protocol.TerminalReport)
        precompile(Protocol.encode_message,(record,))
        precompile(Protocol.decode_terminal_message,(Type{record},String))
        precompile(Protocol.decode_terminal_message,(Type{record},Vector{UInt8}))
    end
    return nothing
end

compile_terminal_paths(::Nothing) = nothing
function compile_terminal_paths(resources::TerminalResources)
    resource = typeof(resources)
    precompile(run_terminal_session!, (resource, TerminalSession, ProfileDefinition, Int, Int))
    for action in (:open,:status,:read,:input,:resize,:keepalive,:disconnect,:stop,:restart)
        precompile(terminal_action!, (Val{action}, resource, Protocol.TerminalCommand))
    end
    return nothing
end

compile_scientific_paths(::Nothing) = nothing
function compile_scientific_paths(resources::ScientificResources)
    resource = typeof(resources)
    for deadline in (Nothing, DateTime)
        precompile(run_scientific_request!, (resource, ScientificHandle, ProfileDefinition,
            Symbol, Dict{String,Any}, String, deadline))
    end
    precompile(execute_assigned!, (resource, Protocol.AssignedJob))
    return nothing
end
