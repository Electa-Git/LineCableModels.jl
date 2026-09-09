"""
    NativeHostCheck

Record local user-manager prerequisites without starting a service. Passing this
inspection is not effective-limit or preparation evidence: the fixed child guard
must still pass before a scientific package is imported.
"""
struct NativeHostCheck
    "Verified local machine/user scope."
    scope::String
    "Resolved local transient-service command."
    command::String
    "Non-root service user."
    uid::Int
    "Non-root primary group."
    gid::Int
    "Fixed prerequisite failure codes."
    failures::Tuple{Vararg{Symbol}}
end
Base.show(io::IO, h::NativeHostCheck) = print(io,"NativeHostCheck(prerequisites=",
    isempty(h.failures) ? "present" : "unavailable",")")

function native_host_failures(uid,read_file,disk_stat)
    failures = Symbol[]
    if disk_stat("/sys/fs/cgroup").ftype != 0x63677270
        push!(failures,:cgroup_v2_required)
    else
        root = "/sys/fs/cgroup/user.slice/user-$uid.slice/user@$uid.service"
        controllers = split(read_file(root * "/cgroup.controllers"))
        for name in ("cpu","memory","pids")
            name in controllers || push!(failures,Symbol(name * "_controller_missing"))
        end
    end
    count = tryparse(Int,strip(read_file("/proc/sys/user/max_user_namespaces")))
    count !== nothing && count > 0 || push!(failures,:user_namespace_unavailable)
    return Tuple(failures)
end

"""
    check_native_host(runner) -> NativeHostCheck

Inspect the fixed local user bus and delegated cgroup-v2 CPU, memory and task
controllers. Require non-root identity and available user namespaces. This is
read-only: no unit, scratch mount or scientific process is created. Missing
controllers disable native launch, not cleanup or public pages.
"""
function check_native_host(runner::CommandRunner;uid=native_uid(),gid=Int(ccall(:getegid,Cuint,())),
        which=Sys.which,read_file=ExecutionCore.isolation_kernel_text,disk_stat=Base.Filesystem.diskstat,
        scope_options...)
    Sys.islinux() || throw(CommandFailure(:linux_required))
    uid > 0 && gid > 0 || throw(CommandFailure(:native_unprivileged_identity_required))
    scope = native_scope(runner;uid,scope_options...)
    command = which("systemd-run")
    command !== nothing || throw(CommandFailure(:native_launcher_unavailable))
    return NativeHostCheck(scope,realpath(command),uid,gid,native_host_failures(uid,read_file,disk_stat))
end

"""
    NativePolicy(profile, receipt)

Derive a fixed native scientific policy from an approved profile and durable
receipt. Construction does not start a service, attest host support or prepare a
model. Terminal profiles cannot use this trusted-code backend.
"""
struct NativePolicy
    "Approved trusted scientific profile."
    profile::ProfileDefinition
    "Exact acquisition intent, recorded before launch."
    receipt::ResourceReceipt
    "The same finite resource dimensions used by container executors."
    limits::ExecutionCore.ExecutorLimits
    function NativePolicy(profile::ProfileDefinition,receipt::ResourceReceipt)
        profile.kind == :scientific && profile.isolation == :trusted_process && receipt.backend == :native ||
            throw(ArgumentError("native scientific profile and receipt required"))
        fence = receipt.fence
        (fence.profile_id,fence.profile_version,fence.fingerprint) ==
            (profile.id,string(profile.version),profile.fingerprint) || throw(ArgumentError("native profile receipt differs"))
        b = profile.budget
        new(profile,receipt,ExecutionCore.ExecutorLimits(b.cpus,b.memory_bytes,b.pids,b.scratch_bytes))
    end
end
Base.show(io::IO,p::NativePolicy) = print(io,"NativePolicy(",p.profile.id,", <owned>)")
Base.show(io::IO,::MIME"text/plain",p::NativePolicy) = show(io,p)

function native_path(path)
    path isa AbstractString && isabspath(path) && normpath(path) == path && ncodeunits(path)<=2048 &&
        !occursin(r"[\x00-\x1f\x7f\$%:]",path) || throw(ArgumentError("unsupported native launch path"))
    return String(path)
end

native_cpu_quota(policy::NativePolicy) = floor(Int,policy.limits.cpus * ExecutionCore.CONTAINER_CPU_PERIOD)

function native_service_properties(policy::NativePolicy,manager::ManagedAgentIdentity)
    manager.worker_id == policy.receipt.fence.worker_id && manager.unit == agent_unit_name(manager.worker_id) ||
        throw(ArgumentError("native manager and receipt differ"))
    b = policy.limits
    scratch = div(b.scratch_bytes-ExecutionCore.CONTAINER_SHM_BYTES,4096)*4096
    return ["Description=" * native_resource_description(policy.receipt),
        "BindsTo=" * manager.unit,"After=" * manager.unit,"Type=exec","ExitType=main","RemainAfterExit=no",
        "Restart=no","KillMode=control-group","KillSignal=SIGINT","SendSIGKILL=yes","FinalKillSignal=SIGKILL",
        "TimeoutStartSec=120","TimeoutStopSec=2","TimeoutStopFailureMode=terminate",
        "CPUAccounting=yes","CPUQuota=$(native_cpu_quota(policy)/1000)%","CPUQuotaPeriodSec=100ms",
        "MemoryAccounting=yes","MemoryMax=$(b.memory_bytes)","MemorySwapMax=0","TasksAccounting=yes","TasksMax=$(b.pids)",
        "PrivateUsers=yes","PrivateDevices=yes","PrivateNetwork=yes","PrivateIPC=yes",
        "ProtectSystem=strict","ProtectControlGroups=yes","ProtectKernelTunables=yes",
        "ReadOnlyPaths=/dev","NoNewPrivileges=yes","CapabilityBoundingSet=","AmbientCapabilities=",
        "TemporaryFileSystem=/tmp:rw,nosuid,nodev,noexec,size=$scratch,mode=1777 /dev/shm:rw,nosuid,nodev,noexec,size=$(ExecutionCore.CONTAINER_SHM_BYTES),mode=1777",
        "WorkingDirectory=/tmp","UMask=0077","LimitNOFILE=1024","LimitCORE=0","LimitMSGQUEUE=0"]
end

function native_exec_arguments(policy::NativePolicy,host::NativeHostCheck,environment::EnvironmentFingerprint;
        project=realpath(policy.profile.environment),julia=realpath(joinpath(Sys.BINDIR,Base.julia_exename())),
        depots=filter(isdir,DEPOT_PATH),env_executable=Sys.which("env"),
        guard=realpath(joinpath(@__DIR__,"..","..","worker","core","src","native-guard.jl")),
        entry=realpath(joinpath(@__DIR__,"..","..","worker","core","src","native-scientific.jl")))
    environment.digest == policy.profile.fingerprint &&
        occursin(r"^[A-Za-z][A-Za-z0-9_]{0,127}$",environment.package) ||
        throw(ArgumentError("native source identity differs"))
    host.scope == policy.receipt.scope && isempty(host.failures) || throw(CommandFailure(:host_prerequisites_unavailable))
    identity = ExecutionCore.NativeIdentity(host.uid,host.gid,native_cgroup(policy.receipt,host.uid))
    paths = native_path.([project,julia,guard,entry])
    # /tmp and /dev/shm are replaced with empty quota-limited mounts. Source,
    # binary and installed depots must remain visible outside these locations.
    visible(path) = !(path == "/tmp" || startswith(path,"/tmp/") || path == "/dev/shm" || startswith(path,"/dev/shm/"))
    all(visible,paths) || throw(ArgumentError("native source is hidden by private scratch"))
    1 <= length(depots) <= 16 || throw(ArgumentError("native execution requires installed depots"))
    approved_depots = native_path.(String.(depots))
    all(visible,approved_depots) || throw(ArgumentError("native depot is hidden by private scratch"))
    executable = native_path(env_executable)
    b = policy.limits
    values = Dict("HOME"=>"/tmp","TMPDIR"=>"/tmp","LANG"=>"C.UTF-8","LC_ALL"=>"C.UTF-8",
        "PATH"=>"/usr/local/bin:/usr/bin:/bin","JULIA_LOAD_PATH"=>"@:@stdlib",
        "JULIA_DEPOT_PATH"=>join(["/tmp/depot";approved_depots],':'),"JULIA_NUM_THREADS"=>"1",
        "OPENBLAS_NUM_THREADS"=>"1","JULIA_PKG_PRECOMPILE_AUTO"=>"0","JULIA_PKG_OFFLINE"=>"true",
        "LCM_NATIVE_CPUS"=>string(b.cpus),"LCM_NATIVE_MEMORY_BYTES"=>string(b.memory_bytes),
        "LCM_NATIVE_PIDS"=>string(b.pids),"LCM_NATIVE_SCRATCH_BYTES"=>string(b.scratch_bytes),
        "LCM_NATIVE_UID"=>string(identity.uid),"LCM_NATIVE_GID"=>string(identity.gid),"LCM_NATIVE_CGROUP"=>identity.cgroup)
    # env -i prevents inherited user-manager credentials, proxy settings, Julia
    # startup paths and DBus addresses from reaching the scientific process.
    return [executable;"-i";[key * "=" * value for (key,value) in sort!(collect(values);by=first)];
        paths[2];"--startup-file=no";"--history-file=no";"--compiled-modules=existing";
        "--threads=1";"--project=" * paths[1];"--load=" * paths[3];paths[4];string(environment.uuid);environment.package]
end

"""
    native_launch_command(policy, manager, host, environment) -> Cmd

Build one fixed transient-service command with bounded private scratch and shared
executor quotas. Attach framed standard IO through systemd-run; never allocate a
PTY here. The scientific package's zero-argument main runs only after the kernel
guard succeeds. The caller must own the receipt before starting this command and
retain it through failed launch, process death and physical cleanup.
"""
function native_launch_command(policy::NativePolicy,manager::ManagedAgentIdentity,host::NativeHostCheck,
        environment::EnvironmentFingerprint;options...)
    arguments = [host.command,"--user","--quiet","--pipe","--wait","--unit=" * resource_name(policy.receipt),"--slice=app.slice"]
    append!(arguments,["--property=" * value for value in native_service_properties(policy,manager)])
    push!(arguments,"--")
    append!(arguments,native_exec_arguments(policy,host,environment;options...))
    values = container_command_environment()
    values["DBUS_SESSION_BUS_ADDRESS"] = native_bus_address(host.uid)
    return setenv(Cmd(arguments),values)
end

"""
    verify_native_service!(journal, runner, policy, manager, host, environment)

Bind a started service to its durable invocation and check the fixed scientific
command, agent lifetime dependency and configured limits. This supplements the
child's mandatory effective-kernel guard; configuration alone cannot admit work.
Return the bound receipt. Missing/replaced units and policy drift fail closed.
"""
function verify_native_service!(journal::ResourceJournal,runner::CommandRunner,policy::NativePolicy,
        manager::ManagedAgentIdentity,host::NativeHostCheck,environment::EnvironmentFingerprint;
        command_options=(;),bus_options=(;))
    receipt = policy.receipt
    identity = inspect_native_unit(runner,receipt;uid=host.uid,bus_options...)
    native_require(identity !== nothing && identity.invocation !== nothing && identity.pid>0 && identity.state==:active)
    receipt = bind_resource!(journal,receipt,identity.invocation)
    path = native_unit_listing(runner,receipt;uid=host.uid,bus_options...)
    unit = native_properties(runner,path,"Unit",["BindsTo","After"];uid=host.uid,bus_options...)
    for key in ("BindsTo","After")
        names = systemd_property(unit,key,"as")
        native_require(names isa AbstractVector && manager.unit in names)
    end
    specifications = [("CPUQuotaPerSecUSec","t",10*native_cpu_quota(policy)),
        ("CPUQuotaPeriodUSec","t",100_000),("MemoryMax","t",policy.limits.memory_bytes),("MemorySwapMax","t",0),
        ("TasksMax","t",policy.limits.pids),("NoNewPrivileges","b",true),("PrivateUsers","b",true),
        ("PrivateDevices","b",true),("PrivateNetwork","b",true),("PrivateIPC","b",true),
        ("ProtectSystem","s","strict"),("ProtectControlGroups","b",true),("ProtectKernelTunables","b",true),
        ("WorkingDirectory","s","/tmp"),("CapabilityBoundingSet","t",0),("AmbientCapabilities","t",0)]
    properties = native_properties(runner,path,"Service",[first.(specifications);"ExecStart";"ExecStartPre";"ExecStartPost"];
        uid=host.uid,bus_options...)
    for (key,type,expected) in specifications
        native_require(container_same(systemd_property(properties,key,type),expected))
    end
    expected = native_exec_arguments(policy,host,environment;command_options...)
    native_require(systemd_exec_matches(systemd_property(properties,"ExecStart","a(sasbttttuii)"),expected))
    for key in ("ExecStartPre","ExecStartPost")
        commands = systemd_property(properties,key,"a(sasbttttuii)")
        native_require(commands isa AbstractVector && isempty(commands))
    end
    after = inspect_native_unit(runner,receipt;uid=host.uid,bus_options...)
    native_require(after !== nothing && after.invocation==identity.invocation && after.pid==identity.pid && after.state==:active)
    return receipt
end
