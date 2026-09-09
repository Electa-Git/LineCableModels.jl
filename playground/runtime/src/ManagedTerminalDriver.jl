"""
    AbstractTerminalDriver

Provide verified terminal profiles, retained PTY acquisition, fixed guarded
startup, exact physical retirement and whole-driver cleanup. These operator hooks
do not authorize a writer, renew a lease or evaluate user input in the agent.
The terminal state owner must enforce current authority around every operation.
"""
abstract type AbstractTerminalDriver end

@required AbstractTerminalDriver begin
    installed_profiles(::AbstractTerminalDriver)
    recover_owned!(::AbstractTerminalDriver)
    verify_terminal!(::AbstractTerminalDriver,::ProfileDefinition,::AssignmentFence)
    terminal_for!(::AbstractTerminalDriver,::ProfileDefinition,::AssignmentFence)
    start_owned_terminal!(::AbstractTerminalDriver,::ProfileDefinition,::AssignmentFence,::Int,::Int)
    terminal_ready_marker(::AbstractTerminalDriver,::ProfileDefinition,::AssignmentFence)
    release_owned!(::AbstractTerminalDriver,::AssignmentFence)
    Base.close(::AbstractTerminalDriver)
end

"""
    verify_terminal!(driver, profile, fence) -> Nothing

Recheck the approved terminal image and created container in the same managed
agent/engine scope. Missing kernel prerequisites are errors, never native fallback.
"""
function verify_terminal! end

"""
    terminal_for!(driver, profile, fence) -> TerminalProcess

Retain a passive PTY handle and stopped-container acquisition under the exact
fence. Repeated calls retrieve that handle; they never start or replace a process.
Partial acquisition remains covered by the physical owner's normal cleanup.
"""
function terminal_for! end

"""
    start_owned_terminal!(driver, profile, fence, columns, rows) -> TerminalProcess

Start the approved container's fixed, guarded Julia REPL using the retained PTY.
This grants no writer authority and does not claim the REPL has finished startup.
Repeated calls never restart an exited/failed process. New sessions require an
explicit higher-level replacement and fresh stream identity.
"""
function start_owned_terminal! end

"""
    terminal_ready_marker(driver, profile, fence) -> Vector{UInt8}

Return the exact acquisition's startup marker. The session owner must observe it
before admitting input and recheck physical policy afterwards. A marker alone is
not authority, an operation result, or evidence that a process remains alive.
"""
function terminal_ready_marker end

function terminal_ready_marker(driver::ManagedResourceDriver,profile::ProfileDefinition,fence::AssignmentFence)
    terminal_driver_profile(driver,profile,fence)
    handle=lock(()->get(driver.handles,fence.lease_id,nothing),driver.lock)
    handle !== nothing || throw(ArgumentError("terminal has not been acquired"))
    return lock(handle.lock) do
        handle.fence==fence && !handle.closing && handle.receipt !== nothing ||
            throw(ArgumentError("terminal acquisition is not current"))
        collect(codeunits("\x1elcm-terminal-ready:"*string(handle.receipt.id)*"\x1f"))
    end
end

function terminal_driver_profile(driver,profile,fence)
    profile.kind==:terminal && profile.isolation==:container ||
        throw(ArgumentError("private terminal requires an approved container profile"))
    return managed_driver_profile(driver,profile,fence)
end

function verify_terminal!(driver::ManagedResourceDriver,profile::ProfileDefinition,fence::AssignmentFence)
    terminal_driver_profile(driver,profile,fence)
    return verify_scientific_backend!(Val(:container),driver,profile,fence)
end

function terminal_for!(driver::ManagedResourceDriver,profile::ProfileDefinition,fence::AssignmentFence)
    terminal_driver_profile(driver,profile,fence)
    handle=managed_handle!(driver,fence)
    return lock(handle.lock) do
        !handle.closing && !driver.closed || throw(ArgumentError("terminal handle is closing"))
        handle.supervisor === nothing || return handle.supervisor::TerminalProcess
        host=verified_driver_host(driver)
        # Shared journal/create policy already validates source, scope, labels,
        # cgroup limits, mounts and image side effects before binding a full ID.
        handle.receipt=create_owned_container!(driver.journal,driver.runner,host,profile,fence)
        handle.supervisor=TerminalProcess()
        return handle.supervisor
    end
end

function terminal_attach_command(host::ContainerHostCheck,receipt::ResourceReceipt)
    receipt.backend==host.engine.name && receipt.physical_id !== nothing &&
        occursin(r"^[a-f0-9]{64}$",receipt.physical_id) && !isempty(host.command) ||
        throw(ArgumentError("terminal attachment requires an exact owned container ID"))
    return setenv(Cmd([host.command;"container";"start";"--attach";"--interactive";
        "--detach-keys=";receipt.physical_id]),container_command_environment())
end

function start_owned_terminal!(driver::ManagedResourceDriver,profile::ProfileDefinition,
        fence::AssignmentFence,columns::Int,rows::Int)
    terminal_size(columns,rows)
    terminal_driver_profile(driver,profile,fence)
    handle=managed_handle!(driver,fence)
    return lock(handle.lock) do
        !handle.closing && !driver.closed || throw(ArgumentError("terminal handle is closing"))
        process=terminal_for!(driver,profile,fence)
        process.started_at==0 || return process
        verify_terminal!(driver,profile,fence)
        host=verified_driver_host(driver)
        # PTY ownership precedes the CLI start; any subsequent failure is still
        # retained in this shared handle and the durable physical receipt.
        start_terminal!(process,terminal_attach_command(host,handle.receipt);columns,rows)
        return process
    end
end

"""Join admitted operations before collecting a borrowed partition's handles."""
mutable struct ManagedViewLifecycle
    "Short admission lock and wakeup for the final in-flight operation."
    condition::Threads.Condition
    "Operations admitted before closure, including partial acquisitions."
    active::Int
    "Whether further operations are permanently forbidden."
    closed::Bool
    "Serialize repeatable partition cleanup without holding the admission lock."
    cleanup_lock::ReentrantLock
end
ManagedViewLifecycle() = ManagedViewLifecycle(Threads.Condition(),0,false,ReentrantLock())

function with_managed_view(action,view)
    lifecycle=view.lifecycle
    lock(lifecycle.condition) do
        lifecycle.closed && throw(ArgumentError("managed resource partition is closed"))
        lifecycle.active+=1
    end
    try
        return action()
    finally
        lock(lifecycle.condition) do
            lifecycle.active-=1
            lifecycle.active==0 && notify(lifecycle.condition;all=true)
        end
    end
end

"""
    ManagedTerminalView(parent)

Borrow the terminal subset of one managed physical owner. Closing this view
releases only terminal handles and never closes the shared journal/command runner.
The root resource owner must close the parent after all partitions have retired.
"""
mutable struct ManagedTerminalView <: AbstractTerminalDriver
    "The same physical owner used by the scientific partition."
    parent::ManagedResourceDriver
    "Only the parent's verified terminal definitions."
    profiles::ProfileRegistry
    "Join admitted operations before partition cleanup."
    lifecycle::ManagedViewLifecycle
end

function managed_profiles(parent,kind)
    profiles=ProfileRegistry()
    for profile in values(installed_profiles(parent).definitions)
        profile.kind==kind && register!(profiles,profile)
    end
    return profiles
end
ManagedTerminalView(parent::ManagedResourceDriver) = ManagedTerminalView(parent,managed_profiles(parent,:terminal),ManagedViewLifecycle())
installed_profiles(view::ManagedTerminalView) = view.profiles
Base.show(io::IO, ::ManagedTerminalView) = print(io,"ManagedTerminalView(<borrowed>)")

function managed_view_profile(view,profile)
    get(view.profiles.definitions,profile.id,nothing) === profile ||
        throw(ArgumentError("managed resource view does not admit this profile"))
    return nothing
end
function recover_owned!(view::ManagedTerminalView)
    return with_managed_view(()->recover_owned!(view.parent),view)
end
function verify_terminal!(view::ManagedTerminalView,profile::ProfileDefinition,fence::AssignmentFence)
    return with_managed_view(view) do
        managed_view_profile(view,profile)
        verify_terminal!(view.parent,profile,fence)
    end
end
function terminal_for!(view::ManagedTerminalView,profile::ProfileDefinition,fence::AssignmentFence)
    return with_managed_view(view) do
        managed_view_profile(view,profile)
        terminal_for!(view.parent,profile,fence)
    end
end
function start_owned_terminal!(view::ManagedTerminalView,profile::ProfileDefinition,fence::AssignmentFence,columns::Int,rows::Int)
    return with_managed_view(view) do
        managed_view_profile(view,profile)
        start_owned_terminal!(view.parent,profile,fence,columns,rows)
    end
end
function terminal_ready_marker(view::ManagedTerminalView,profile::ProfileDefinition,fence::AssignmentFence)
    return with_managed_view(view) do
        managed_view_profile(view,profile)
        terminal_ready_marker(view.parent,profile,fence)
    end
end

function release_managed_view!(view,fence)
    profile=get(view.profiles.definitions,fence.profile_id,nothing)
    profile !== nothing && string(profile.version)==fence.profile_version && profile.fingerprint==fence.fingerprint ||
        throw(ArgumentError("cleanup fence belongs to another resource partition"))
    return release_owned!(view.parent,fence)
end
release_owned!(view::ManagedTerminalView,fence::AssignmentFence) = release_managed_view!(view,fence)

function close_managed_view!(view)
    lock(view.lifecycle.cleanup_lock) do
        lock(view.lifecycle.condition) do
            view.lifecycle.closed=true
            while view.lifecycle.active>0
                wait(view.lifecycle.condition)
            end
        end
        fences=lock(()->[h.fence for h in values(view.parent.handles)
            if haskey(view.profiles.definitions,h.fence.profile_id)],view.parent.lock)
        completed=true
        for fence in fences
            completed &= try release_owned!(view,fence) === true catch; false end
        end
        completed || throw(ArgumentError("managed partition cleanup remains unresolved"))
    end
    return nothing
end
Base.close(view::ManagedTerminalView) = close_managed_view!(view)

"""
    ManagedScientificView(parent)

Borrow only verified scientific definitions from the shared physical owner.
Partition cleanup never closes terminal handles or the parent journal. Existing
ScientificResources retains all scientific lease/preparation/operation checks.
"""
mutable struct ManagedScientificView <: AbstractScientificDriver
    "Shared physical owner, not a separately opened journal."
    parent::ManagedResourceDriver
    "Only verified scientific definitions."
    profiles::ProfileRegistry
    "Join admitted operations before partition cleanup."
    lifecycle::ManagedViewLifecycle
end
ManagedScientificView(parent::ManagedResourceDriver) = ManagedScientificView(parent,managed_profiles(parent,:scientific),ManagedViewLifecycle())
installed_profiles(view::ManagedScientificView) = view.profiles
Base.show(io::IO, ::ManagedScientificView) = print(io,"ManagedScientificView(<borrowed>)")
function recover_owned!(view::ManagedScientificView)
    return with_managed_view(()->recover_owned!(view.parent),view)
end
function verify_executor!(view::ManagedScientificView,profile::ProfileDefinition,fence::AssignmentFence)
    return with_managed_view(view) do
        managed_view_profile(view,profile)
        verify_executor!(view.parent,profile,fence)
    end
end
function executor_for!(view::ManagedScientificView,profile::ProfileDefinition,fence::AssignmentFence)
    return with_managed_view(view) do
        managed_view_profile(view,profile)
        executor_for!(view.parent,profile,fence)
    end
end
release_owned!(view::ManagedScientificView,fence::AssignmentFence) = release_managed_view!(view,fence)
Base.close(view::ManagedScientificView) = close_managed_view!(view)
