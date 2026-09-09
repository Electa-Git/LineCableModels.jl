"""
    ManagedAgentIdentity

Bind an agent to its actual user-systemd service incarnation and fixed post-stop
recovery command. This is process-lifetime evidence, not a resource-limit check,
broker identity, lease grant or executor readiness report.
"""
struct ManagedAgentIdentity
    "Generated unit name for this provisioned worker."
    unit::String
    "Systemd's exact invocation UUID, not a saved process ID."
    invocation::String
    "Observed service control group containing this agent."
    cgroup::String
    "Resolved operator configuration path used by the fixed start command."
    config_file::String
    "Fixed post-stop resource journal, independent of later config-file edits."
    journal_root::String
    "Provisioned worker identity retained by post-stop recovery."
    worker_id::String
end
Base.show(io::IO, identity::ManagedAgentIdentity) = print(io,"ManagedAgentIdentity(",identity.unit,", <private>)")
Base.show(io::IO, ::MIME"text/plain", identity::ManagedAgentIdentity) = show(io,identity)

managed_require(condition) = condition === true ? nothing : throw(CommandFailure(:managed_agent_unverified))
agent_unit_name(worker::AbstractString) = "lcm-agent-" * Protocol.runtime_token(worker) * ".service"
agent_journal_root(config::AgentConfig) = joinpath(config.scratch_root,"resources")

function agent_julia_command(action::String,arguments::Vector{String})
    executable = realpath(joinpath(Sys.BINDIR,Base.julia_exename()))
    project = realpath(joinpath(@__DIR__,".."))
    return [executable,"--startup-file=no","--compiled-modules=existing","--threads=2",
        "--project=" * project,"-m","LineCableModelsRuntime","runtime",action,arguments...]
end
agent_start_command(path::String) = agent_julia_command("start-agent",["--config",path])
agent_recovery_command(config::AgentConfig) = agent_julia_command("recover-agent",
    ["--journal",agent_journal_root(config),"--worker",config.worker_id])

function systemd_quote(argument::AbstractString)
    ncodeunits(argument) <= 4096 && !occursin(r"[\x00-\x1f\x7f]",argument) ||
        throw(ArgumentError("service command argument is invalid"))
    return "\"" * replace(argument,"\\"=>"\\\\","\""=>"\\\"","%"=>"%%", "\$"=>"\$\$") * "\""
end
systemd_command(arguments::Vector{String}) = join(systemd_quote.(arguments)," ")

"""
    agent_service_unit(config_file) -> String

Render an operator-installable user service without installing or starting it.
The unit runs the lightweight agent and invokes exact journal/worker recovery
after normal exit, crash or forced termination. Recovery arguments capture the
original journal location; editing the config later cannot redirect cleanup.
Both shutdown phases have finite timeouts. Automatic service restart is disabled.
"""
function agent_service_unit(config_file::AbstractString)
    path = realpath(config_file)
    config = read_agent_config(path)
    return """
    [Unit]
    Description=LCM runtime agent $(config.worker_id)

    [Service]
    Type=exec
    ExitType=main
    RemainAfterExit=no
    ExecStart=$(systemd_command(agent_start_command(path)))
    ExecStopPost=$(systemd_command(agent_recovery_command(config)))
    KillMode=control-group
    KillSignal=SIGINT
    SendSIGKILL=yes
    FinalKillSignal=SIGKILL
    TimeoutStartSec=120
    TimeoutStopSec=60
    TimeoutStopFailureMode=terminate
    Restart=no
    UMask=0077
    StandardInput=null
    StandardOutput=journal
    StandardError=journal

    [Install]
    WantedBy=default.target
    """
end

function systemd_bus_command(runner,arguments;address=nothing,invoke=arguments->run_owned_command!(runner,
        setenv(Cmd(arguments),container_command_environment())),which=Sys.which)
    executable = which("busctl")
    managed_require(executable !== nothing)
    endpoint = address === nothing ? "--user" : "--address=" * address
    result = invoke([[executable,endpoint,"--no-pager","--timeout=5","--json=short"];arguments])
    managed_require(result isa CommandResult && result.exitcode == 0 && ncodeunits(result.output) <= 256*1024)
    lines = split(strip(result.output),'\n';keepempty=false)
    managed_require(1 <= length(lines) <= 32)
    records = [try JSON3.read(line) catch; nothing end for line in lines]
    managed_require(all(r->r isa JSON3.Object && haskey(r,:type) && haskey(r,:data),records))
    return records
end

function systemd_properties(runner,path,interface,properties;options...)
    records = systemd_bus_command(runner,["get-property","org.freedesktop.systemd1",path,
        "org.freedesktop.systemd1." * interface,properties...];options...)
    managed_require(length(records) == length(properties))
    return Dict(key=>record for (key,record) in zip(properties,records))
end

# Observe only the already-verified service incarnation. Loss/replacement is an
# error; only its explicit deactivating state requests cooperative root cleanup.
function managed_agent_stopping(runner,object,identity::ManagedAgentIdentity;options...)
    properties = systemd_properties(runner,object,"Unit",["Id","InvocationID","ActiveState"];options...)
    invocation = systemd_property(properties,"InvocationID","ay")
    managed_require(systemd_property(properties,"Id","s")==identity.unit &&
        invocation isa AbstractVector && length(invocation)==16 &&
        all(v->v isa Integer && !(v isa Bool) && 0<=v<=255,invocation) &&
        bytes2hex(UInt8.(invocation))==identity.invocation)
    state = systemd_property(properties,"ActiveState","s")
    managed_require(state in ("active","deactivating"))
    return state == "deactivating"
end

function systemd_property(properties,key,type)
    record = get(properties,key,nothing)
    managed_require(record isa JSON3.Object && get(record,:type,nothing) == type)
    return record.data
end

function systemd_exec_matches(record,expected)
    record isa AbstractVector && length(record) == 1 || return false
    entry = only(record)
    entry isa AbstractVector && length(entry) == 10 || return false
    return entry[1] == first(expected) && entry[2] == expected && entry[3] === false
end

function managed_cgroup(text,expected)
    ncodeunits(text) <= 4096 || return false
    return strip(text) == "0::" * expected
end

"""
    verify_managed_agent(runner, config_file, config; previous=nothing)

Inspect typed service-manager properties and this process's cgroup. Require the
exact worker unit, current PID, nonzero invocation identity, active service,
fixed Julia start/recovery commands, complete-group termination and finite stop
timeouts. A flag or environment variable is not accepted as proof. With previous,
also require the same service incarnation and journal identity on recheck.
Return ManagedAgentIdentity, or a fixed CommandFailure without private bus output.
"""
function verify_managed_agent(runner::CommandRunner,config_file::AbstractString,config::AgentConfig;
        previous::Union{Nothing,ManagedAgentIdentity}=nothing,pid=getpid(),
        read_cgroup=()->read("/proc/self/cgroup",String),options...)
    path = realpath(config_file)
    unit = agent_unit_name(config.worker_id)
    records = systemd_bus_command(runner,["call","org.freedesktop.systemd1","/org/freedesktop/systemd1",
        "org.freedesktop.systemd1.Manager","GetUnit","s",unit];options...)
    managed_require(length(records) == 1 && records[1].type == "o" &&
        records[1].data isa AbstractVector && length(records[1].data) == 1)
    object_path = only(records[1].data)
    managed_require(object_path isa String && occursin(r"^/org/freedesktop/systemd1/unit/[A-Za-z0-9_]+$",object_path))
    properties = systemd_properties(runner,object_path,"Unit",["Id","InvocationID","ActiveState"];options...)
    managed_require(systemd_property(properties,"Id","s") == unit &&
        systemd_property(properties,"ActiveState","s") == "active")
    invocation = systemd_property(properties,"InvocationID","ay")
    managed_require(invocation isa AbstractVector && length(invocation) == 16 &&
        all(v->v isa Integer && !(v isa Bool) && 0 <= v <= 255,invocation) && any(!iszero,invocation))
    service = systemd_properties(runner,object_path,"Service",["MainPID","Type","ExitType","ControlGroup",
        "KillMode","KillSignal","SendSIGKILL","FinalKillSignal","TimeoutStopUSec","Restart",
        "RemainAfterExit","TimeoutStopFailureMode","ExecStart","ExecStop","ExecStopPost"];options...)
    main_pid = systemd_property(service,"MainPID","u")
    managed_require(main_pid isa Integer && !(main_pid isa Bool) && main_pid == pid && pid > 0 &&
        systemd_property(service,"Type","s") == "exec" && systemd_property(service,"ExitType","s") == "main" &&
        systemd_property(service,"KillMode","s") == "control-group" && systemd_property(service,"KillSignal","i") == 2 &&
        systemd_property(service,"SendSIGKILL","b") === true && systemd_property(service,"FinalKillSignal","i") == 9 &&
        systemd_property(service,"TimeoutStopUSec","t") == 60_000_000 && systemd_property(service,"Restart","s") == "no" &&
        systemd_property(service,"RemainAfterExit","b") === false &&
        systemd_property(service,"TimeoutStopFailureMode","s") == "terminate" &&
        isempty(systemd_property(service,"ExecStop","a(sasbttttuii)")))
    cgroup = systemd_property(service,"ControlGroup","s")
    managed_require(cgroup isa String && startswith(cgroup,"/user.slice/") && endswith(cgroup,"/" * unit) &&
        !occursin(r"[\\\x00-\x20]",cgroup) && managed_cgroup(read_cgroup(),cgroup))
    managed_require(systemd_exec_matches(systemd_property(service,"ExecStart","a(sasbttttuii)"),agent_start_command(path)))
    managed_require(systemd_exec_matches(systemd_property(service,"ExecStopPost","a(sasbttttuii)"),agent_recovery_command(config)))
    identity = ManagedAgentIdentity(unit,bytes2hex(UInt8.(invocation)),cgroup,path,agent_journal_root(config),config.worker_id)
    previous === nothing || managed_require(all(k->getfield(identity,k) == getfield(previous,k),fieldnames(ManagedAgentIdentity)))
    return identity
end

"""
    recover_agent_resources!(journal_root, worker_id)

Recover the explicit worker's persisted resources without broker credentials or a
possibly edited agent configuration. No absent journal is created. A live owner
holds the kernel journal lock and prevents this cleanup from racing its launches.
Attempt each recorded backend, retain unresolved receipts, and never infer absence
from a missing engine or service manager.
"""
function recover_agent_resources!(journal_root::AbstractString,worker_id::AbstractString)
    worker = Protocol.runtime_token(worker_id)
    root = abspath(journal_root)
    root in ("/",homedir(),pwd(),dirname(pwd()),tempdir()) && throw(ArgumentError("recovery requires a dedicated journal"))
    ispath(root) || return nothing
    isfile(joinpath(root,"owner.json")) || throw(ArgumentError("resource journal ownership is unavailable"))
    journal = ResourceJournal(root,worker;capacity=256)
    runner = CommandRunner()
    try
        recover_agent_journal!(journal,runner)
    finally
        close(runner)
        close(journal)
    end
    return nothing
end

function recover_agent_journal!(journal,runner)
    completed = true
    for backend in unique(r.backend for r in resource_receipts(journal))
        if backend in (:podman,:docker)
            try
                host = check_container_host(runner;requested=string(backend))
                recover_containers!(journal,runner,host)
            catch
                completed = false
            end
        elseif backend == :native
            try recover_native!(journal,runner) catch; completed = false end
        else
            completed = false
        end
    end
    completed && isempty(resource_receipts(journal)) || throw(ArgumentError("owned agent recovery remains unresolved"))
    return nothing
end
