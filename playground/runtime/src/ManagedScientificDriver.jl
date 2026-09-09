"""Retain exact partial acquisition and process ownership for one managed lease."""
mutable struct ManagedScientificHandle
    "Complete assignment identity."
    fence::AssignmentFence
    "Bound physical receipt; the journal also retains earlier failed intents."
    receipt::Union{Nothing,ResourceReceipt}
    "Owned attached CLI/process framing, not another executor shared across runs."
    supervisor::Union{Nothing,ScientificProcess,TerminalProcess}
    "Whether new process creation has been revoked."
    closing::Bool
    "Serialize acquisition/release for this lease only."
    lock::ReentrantLock
end

"""
    ManagedScientificDriver(config_file, config)

Open the managed agent's private physical owner and inspect its approved local
native and container profiles. No numerical process is started. Unsupported/missing profiles
are recorded as unavailable; installed_profiles returns only the verified subset.
The actual service's post-stop recovery must be present before opening ownership.

ScientificResources retains lease authority and preparation orchestration. This
driver owns stopped acquisition, the attached executor command and exact physical
retirement. Source/image installation and arbitrary container flags are not hooks.
"""
mutable struct ManagedScientificDriver <: AbstractScientificDriver
    "Operator configuration, not browser state."
    config::AgentConfig
    "Verified service incarnation and fixed recovery target."
    manager::ManagedAgentIdentity
    "Only locally verified scientific definitions."
    profiles::ProfileRegistry
    "Fixed per-profile unavailability codes."
    unavailable::Dict{String,Symbol}
    "Shared bounded owner for finite engine/manager commands."
    runner::CommandRunner
    "Exclusive resource journal, retained through unresolved cleanup."
    journal::ResourceJournal
    "Inspected local engine; nothing when it is unavailable."
    host::Union{Nothing,ContainerHostCheck}
    "Inspected engine scope, not a context-name guess."
    scope::Union{Nothing,String}
    "Inspected local native-service prerequisites, absent when unused/unavailable."
    native_host::Union{Nothing,NativeHostCheck}
    "Per-lease handles, including partial acquisitions."
    handles::Dict{String,ManagedScientificHandle}
    "Whether old ownership was reconciled."
    recovered::Bool
    "Whether new work is permanently forbidden."
    closed::Bool
    "Whether physical teardown and local handle closure completed."
    cleanup_complete::Bool
    "Short handle/admission updates, never numerical work."
    lock::ReentrantLock
    "Serialize whole-owner recovery/closure."
    lifecycle_lock::ReentrantLock
end
Base.show(io::IO, driver::ManagedScientificDriver) = print(io,"ManagedScientificDriver(",driver.config.worker_id,", <owned>)")
Base.show(io::IO, ::MIME"text/plain", driver::ManagedScientificDriver) = show(io,driver)

function ManagedScientificDriver(config_file::AbstractString,config::AgentConfig)
    runner = CommandRunner()
    journal = nothing
    try
        manager = verify_managed_agent(runner,config_file,config)
        journal = ResourceJournal(manager.journal_root,config.worker_id;capacity=config.capacity)
        needs_container = any(p->p.isolation==:container,values(config.profiles.definitions))
        needs_native = any(p->p.kind==:scientific && p.isolation==:trusted_process,values(config.profiles.definitions))
        host = needs_container ? try check_container_host(runner;requested=string(config.container_runtime)) catch; nothing end : nothing
        native_host = needs_native ? try check_native_host(runner) catch; nothing end : nothing
        scope = host === nothing ? nothing : try container_scope(runner,host) catch; nothing end
        profiles,unavailable = ProfileRegistry(),Dict{String,Symbol}()
        for profile in values(config.profiles.definitions)
            failure = if profile.isolation == :trusted_process
                if native_host === nothing
                    :native_host_unavailable
                elseif !isempty(native_host.failures)
                    first(native_host.failures)
                else
                    try
                        verify_native_environment(profile)
                        ExecutionCore.ExecutorLimits(profile.budget.cpus,profile.budget.memory_bytes,
                            profile.budget.pids,profile.budget.scratch_bytes)
                        nothing
                    catch
                        :native_environment_unavailable
                    end
                end
            elseif host === nothing || scope === nothing
                :container_engine_unavailable
            elseif !isempty(host.failures)
                first(host.failures)
            else
                try
                    object = container_inspection(runner,host,["image","inspect","--format","{{json .}}",profile.environment];
                        invoke=args->run_owned_command!(runner,setenv(Cmd(args),container_command_environment())))
                    verify_container_image(profile,object)
                    nothing
                catch error
                    error isa ExecutionCore.IsolationError ? error.code : :image_unavailable
                end
            end
            failure === nothing ? register!(profiles,profile) : (unavailable[profile.id] = failure)
        end
        return ManagedScientificDriver(config,manager,profiles,unavailable,runner,journal,host,scope,native_host,
            Dict{String,ManagedScientificHandle}(),false,false,false,ReentrantLock(),ReentrantLock())
    catch
        journal === nothing || close(journal)
        close(runner)
        rethrow()
    end
end
installed_profiles(driver::ManagedScientificDriver) = driver.profiles

function recover_owned!(driver::ManagedScientificDriver)
    lock(driver.lifecycle_lock) do
        driver.closed && throw(ArgumentError("scientific driver is closed"))
        driver.recovered && return nothing
        isempty(driver.handles) || throw(ArgumentError("scientific driver has unreconciled live handles"))
        verify_managed_agent(driver.runner,driver.manager.config_file,driver.config;previous=driver.manager)
        recover_agent_journal!(driver.journal,driver.runner)
        driver.recovered = true
    end
    return nothing
end

function managed_driver_profile(driver,profile,fence)
    driver.closed && throw(ArgumentError("scientific driver is closed"))
    driver.recovered || throw(ArgumentError("scientific driver is not recovered"))
    get(driver.profiles.definitions,profile.id,nothing) === profile &&
        fence.worker_id == driver.config.worker_id && (fence.profile_id,fence.profile_version,fence.fingerprint) ==
        (profile.id,string(profile.version),profile.fingerprint) || throw(ArgumentError("scientific assignment is not approved"))
    return nothing
end

function scientific_driver_profile(driver,profile,fence)
    profile.kind==:scientific || throw(ArgumentError("scientific profile required"))
    return managed_driver_profile(driver,profile,fence)
end

function verified_driver_host(driver)
    driver.host !== nothing && driver.scope !== nothing || throw(CommandFailure(:container_engine_unavailable))
    verify_managed_agent(driver.runner,driver.manager.config_file,driver.config;previous=driver.manager)
    host,scope = recheck_container_host(driver.runner,driver.host)
    container_require(isempty(host.failures),:host_prerequisites_unavailable)
    container_require(scope == driver.scope,:container_scope_changed)
    return host
end

function verify_executor!(driver::ManagedScientificDriver,profile::ProfileDefinition,fence::AssignmentFence)
    scientific_driver_profile(driver,profile,fence)
    verify_scientific_backend!(Val(profile.isolation),driver,profile,fence)
end

function verify_scientific_backend!(::Val{:container},driver,profile,fence)
    host = verified_driver_host(driver)
    invoke = args->run_owned_command!(driver.runner,setenv(Cmd(args),container_command_environment()))
    image = verify_container_image(profile,container_inspection(driver.runner,host,
        ["image","inspect","--format","{{json .}}",profile.environment];invoke))
    handle = lock(()->get(driver.handles,fence.lease_id,nothing),driver.lock)
    if handle !== nothing
        lock(handle.lock) do
            handle.fence == fence && !handle.closing || throw(ArgumentError("scientific handle no longer accepts work"))
            if handle.receipt !== nothing
                policy = ContainerPolicy(profile,handle.receipt)
                object = container_inspection(driver.runner,host,["container","inspect","--format","{{json .}}",handle.receipt.physical_id];invoke)
                process = handle.supervisor === nothing ? nothing : handle.supervisor.process
                running = process !== nothing && process_running(process)
                verify_created_container(policy,object,image;running)
            end
        end
    end
    return nothing
end

function managed_handle!(driver::ManagedScientificDriver,fence::AssignmentFence)
    return lock(driver.lock) do
        driver.closed && throw(ArgumentError("managed resource driver is closed"))
        found = get(driver.handles,fence.lease_id,nothing)
        if found === nothing
            length(driver.handles) < driver.config.capacity || throw(CapacityUnavailable())
            found = ManagedScientificHandle(fence,nothing,nothing,false,ReentrantLock())
            driver.handles[fence.lease_id] = found
        end
        found.fence == fence || throw(ArgumentError("managed handle fence differs"))
        found
    end
end

function executor_for!(driver::ManagedScientificDriver,profile::ProfileDefinition,fence::AssignmentFence)
    scientific_driver_profile(driver,profile,fence)
    handle = managed_handle!(driver,fence)
    return lock(handle.lock) do
        !handle.closing && !driver.closed || throw(ArgumentError("scientific handle is closing"))
        handle.supervisor === nothing || return handle.supervisor::ScientificProcess
        acquire_scientific_backend!(Val(profile.isolation),driver,handle,profile,fence)
        return handle.supervisor
    end
end

stop_managed_process!(process::ScientificProcess) = ExecutionCore.stop_executor!(process)
stop_managed_process!(process::TerminalProcess) = close(process)


function verified_native_driver_host(driver)
    driver.native_host !== nothing || throw(CommandFailure(:native_host_unavailable))
    verify_managed_agent(driver.runner,driver.manager.config_file,driver.config;previous=driver.manager)
    host = check_native_host(driver.runner)
    isempty(host.failures) || throw(CommandFailure(:host_prerequisites_unavailable))
    previous = driver.native_host
    (host.scope,host.command,host.uid,host.gid) == (previous.scope,previous.command,previous.uid,previous.gid) ||
        throw(CommandFailure(:native_scope_changed))
    return host
end

function verify_scientific_backend!(::Val{:trusted_process},driver,profile,fence)
    host = verified_native_driver_host(driver)
    environment = verify_native_environment(profile)
    handle = lock(()->get(driver.handles,fence.lease_id,nothing),driver.lock)
    if handle !== nothing
        lock(handle.lock) do
            handle.fence == fence && !handle.closing || throw(ArgumentError("scientific handle no longer accepts work"))
            if handle.receipt !== nothing && handle.supervisor !== nothing && handle.supervisor.process !== nothing
                handle.receipt = verify_native_service!(driver.journal,driver.runner,NativePolicy(profile,handle.receipt),
                    driver.manager,host,environment)
            end
        end
    end
    return nothing
end

function acquire_scientific_backend!(::Val{:container},driver,handle,profile,fence)
    host = verified_driver_host(driver)
    handle.receipt = create_owned_container!(driver.journal,driver.runner,host,profile,fence)
    command = setenv(Cmd([host.command;"container";"start";"--attach";"--interactive";"--detach-keys=";handle.receipt.physical_id]),
        container_command_environment())
    handle.supervisor = ScientificProcess(;command,discard_stderr=true,startup_timeout_seconds=120)
    return nothing
end

function acquire_scientific_backend!(::Val{:trusted_process},driver,handle,profile,fence)
    host = verified_native_driver_host(driver)
    environment = verify_native_environment(profile)
    # The intent precedes even construction of the command. Every later failure
    # remains covered by lease release and the same agent post-stop journal.
    handle.receipt = reserve_resource!(driver.journal,fence,:native,host.scope)
    command = native_launch_command(NativePolicy(profile,handle.receipt),driver.manager,host,environment)
    handle.supervisor = ScientificProcess(;command,discard_stderr=true,startup_timeout_seconds=120)
    return nothing
end

function release_owned!(driver::ManagedScientificDriver,fence::AssignmentFence)
    driver.cleanup_complete && return true
    handle = lock(driver.lock) do
        found = get(driver.handles,fence.lease_id,nothing)
        if found !== nothing
            found.fence == fence || throw(ArgumentError("scientific cleanup fence differs"))
            found.closing = true
        end
        found
    end
    function release()
        supervisor = handle === nothing ? nothing : handle.supervisor
        process_complete = supervisor === nothing || try
            stop_managed_process!(supervisor)
            true
        catch
            false
        end
        receipts = filter(r->r.fence.lease_id == fence.lease_id,resource_receipts(driver.journal))
        all(r->r.fence == fence,receipts) || throw(ArgumentError("scientific receipt fence differs"))
        completed = true
        for receipt in receipts
            try
                if receipt.backend == :native
                    completed &= remove_owned_native!(driver.journal,driver.runner,receipt)
                else
                    host = check_container_host(driver.runner;requested=string(receipt.backend))
                    completed &= remove_owned_container!(driver.journal,driver.runner,host,receipt)
                end
            catch
                completed = false
            end
        end
        # Physical removal is still attempted if an attached CLI/pipe could not
        # retire. Its original handle remains owned until a joined retry passes.
        if !process_complete && completed
            process_complete = try
                stop_managed_process!(supervisor)
                true
            catch
                false
            end
        end
        completed &= process_complete
        if completed
            lock(driver.lock) do
                get(driver.handles,fence.lease_id,nothing) === handle && delete!(driver.handles,fence.lease_id)
            end
        end
        return completed
    end
    return handle === nothing ? release() : lock(release,handle.lock)
end

function Base.close(driver::ManagedScientificDriver)
    lock(driver.lifecycle_lock) do
        driver.cleanup_complete && return nothing
        lock(driver.lock) do
            driver.closed = true
        end
        fences = lock(()->[h.fence for h in values(driver.handles)],driver.lock)
        completed = true
        for fence in fences
            completed &= try release_owned!(driver,fence) catch; false end
        end
        # Includes startup leftovers and acquisition failures without a handle.
        try recover_agent_journal!(driver.journal,driver.runner) catch; completed = false end
        completed && isempty(driver.handles) && isempty(resource_receipts(driver.journal)) ||
            throw(ArgumentError("scientific driver cleanup remains unresolved"))
        close(driver.runner)
        close(driver.journal)
        driver.cleanup_complete = true
    end
    return nothing
end

"""Compatibility name; native and container profiles share one journal owner."""
const ContainerScientificDriver = ManagedScientificDriver
const ContainerScientificHandle = ManagedScientificHandle
const container_driver_profile = scientific_driver_profile

"""Shared physical owner; the earlier scientific-only name remains compatible."""
const ManagedResourceDriver = ManagedScientificDriver
const ManagedProcessHandle = ManagedScientificHandle
