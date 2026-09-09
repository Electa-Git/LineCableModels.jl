"""
    ContainerPolicy(profile::ProfileDefinition, receipt::ResourceReceipt)

Derive one immutable launch policy from an approved container profile and its
durable ownership receipt. Scientific and terminal processes use identical kernel
limits; only the fixed Julia entry command and terminal allocation differ.
Construction performs no engine action and proves neither host support nor
preparation. No browser-supplied command, mount, environment or image is accepted.
"""
struct ContainerPolicy
    "Approved digest-pinned environment and operation contract."
    profile::ProfileDefinition
    "Exact acquisition identity, recorded before create."
    receipt::ResourceReceipt
    "Mandatory effective kernel limits, also checked inside Julia."
    limits::ExecutionCore.ContainerLimits
    function ContainerPolicy(profile::ProfileDefinition, receipt::ResourceReceipt)
        profile.isolation == :container || throw(ArgumentError("container profile required"))
        receipt.backend in (:docker,:podman) || throw(ArgumentError("container receipt required"))
        fence = receipt.fence
        (fence.profile_id, fence.profile_version, fence.fingerprint) ==
            (profile.id, string(profile.version), profile.fingerprint) ||
            throw(ArgumentError("container profile does not match its receipt"))
        b = profile.budget
        limits = ExecutionCore.ContainerLimits(b.cpus,b.memory_bytes,b.pids,b.scratch_bytes)
        new(profile,receipt,limits)
    end
end
Base.show(io::IO, p::ContainerPolicy) = print(io,"ContainerPolicy(",p.profile.id,", ",p.profile.kind,", <owned>)")

container_require(condition, code::Symbol) = condition === true ? nothing : throw(ExecutionCore.IsolationError(code))
container_empty(value) = value === nothing || (value isa Union{AbstractVector,AbstractDict,JSON3.Object} && isempty(value))
container_json(value) = value isa Union{AbstractDict,JSON3.Object}
container_vector(value) = value isa AbstractVector && all(v -> v isa AbstractString, value)
container_same(actual, expected) = expected isa Bool ? actual === expected :
    expected isa Integer ? (actual isa Integer && !(actual isa Bool) && actual == expected) : actual == expected

function container_environment(policy::ContainerPolicy)
    environment=container_environment(policy.profile)
    policy.profile.kind==:terminal && (environment["LCM_TERMINAL_READY"]=string(policy.receipt.id))
    return environment
end
function container_environment(profile::ProfileDefinition)
    b = profile.budget
    limits = ExecutionCore.ContainerLimits(b.cpus,b.memory_bytes,b.pids,b.scratch_bytes)
    values = Dict("HOME"=>"/tmp", "TMPDIR"=>"/tmp", "LANG"=>"C.UTF-8", "LC_ALL"=>"C.UTF-8", "HOSTNAME"=>"lcm-executor",
        "PATH"=>"/usr/local/julia/bin:/usr/local/bin:/usr/bin:/bin", "JULIA_DEPOT_PATH"=>"/tmp/depot:/opt/lcm/depot",
        "JULIA_LOAD_PATH"=>"@:@stdlib", "JULIA_NUM_THREADS"=>"1", "OPENBLAS_NUM_THREADS"=>"1",
        "JULIA_PKG_PRECOMPILE_AUTO"=>"0", "JULIA_PKG_OFFLINE"=>"true",
        "LCM_CONTAINER_CPUS"=>string(limits.cpus), "LCM_CONTAINER_MEMORY_BYTES"=>string(limits.memory_bytes),
        "LCM_CONTAINER_PIDS"=>string(limits.pids), "LCM_CONTAINER_SCRATCH_BYTES"=>string(limits.scratch_bytes))
    profile.kind == :terminal && (values["TERM"] = "xterm-256color")
    return values
end

function container_julia_arguments(policy::ContainerPolicy)
    arguments = ["--startup-file=no", "--history-file=no", "--compiled-modules=existing",
        "--project=/opt/lcm/source/playground/worker/profiles/active", "--load=/opt/lcm/container-guard.jl"]
    append!(arguments, policy.profile.kind == :terminal ? ["--interactive", "--color=yes"] : ["/opt/lcm/scientific.jl"])
    return arguments
end

container_tmpfs_bytes(policy::ContainerPolicy) =
    div(policy.limits.scratch_bytes - ExecutionCore.CONTAINER_SHM_BYTES,4096)*4096
container_tmpfs_options(policy::ContainerPolicy) =
    "rw,nosuid,nodev,noexec,size=$(container_tmpfs_bytes(policy)),mode=1777" *
        (policy.receipt.backend == :podman ? ",notmpcopyup" : "")

"""
    container_create_arguments(policy::ContainerPolicy) -> Vector{String}

Render the shared create-only policy as operator CLI arguments. This pure function
does not authorize launching: create_owned_container! additionally requires a
passing host, a verified installed image, fresh scope and a current journal entry.
Docker's private PID/UTS namespaces use empty mode values; Podman accepts private.
The tmpfs allowance rounds down to pages and includes the separate shared memory.
"""
function container_create_arguments(policy::ContainerPolicy)
    limits = policy.limits
    kind = policy.receipt.backend
    arguments = ["container", "create", "--name", resource_name(policy.receipt),
        "--pull=never", "--read-only", "--network=none", "--ipc=private", "--cgroupns=private",
        kind == :podman ? "--pid=private" : "--pid=", kind == :podman ? "--uts=private" : "--uts=",
        "--user=1000:1000", "--cap-drop=ALL", "--security-opt=no-new-privileges",
        "--log-driver=none", "--restart=no", "--workdir=/tmp", "--hostname=lcm-executor", "--interactive",
        "--cpu-period=$(ExecutionCore.CONTAINER_CPU_PERIOD)",
        "--cpu-quota=$(floor(Int,limits.cpus*ExecutionCore.CONTAINER_CPU_PERIOD))",
        "--memory=$(limits.memory_bytes)", "--memory-swap=$(limits.memory_bytes)",
        "--pids-limit=$(limits.pids)", "--shm-size=$(ExecutionCore.CONTAINER_SHM_BYTES)",
        "--ulimit=core=0:0", "--ulimit=msgqueue=0:0", "--ulimit=nofile=1024:1024",
        "--tmpfs", "/tmp:" * container_tmpfs_options(policy), "--entrypoint=/usr/local/julia/bin/julia"]
    policy.profile.kind == :terminal && push!(arguments,"--tty")
    kind == :podman && append!(arguments,["--read-only-tmpfs=false","--image-volume=ignore","--http-proxy=false"])
    for (key,value) in sort!(collect(resource_labels(policy.receipt));by=first)
        append!(arguments,["--label",key * "=" * value])
    end
    for (key,value) in sort!(collect(container_environment(policy));by=first)
        append!(arguments,["--env",key * "=" * value])
    end
    push!(arguments,policy.profile.environment)
    append!(arguments,container_julia_arguments(policy))
    return arguments
end

function container_env_map(value)
    container_require(container_vector(value) && length(value) <= 64,:environment_unverified)
    result = Dict{String,String}()
    for entry in value
        parts = split(entry,'=';limit=2)
        container_require(length(parts) == 2 && (occursin(r"^[A-Z][A-Z0-9_]*$",parts[1]) || parts[1] == "container") &&
            !haskey(result,parts[1]) && ncodeunits(entry) <= 4096 && !occursin(r"[\x00\r\n]",entry),:environment_unverified)
        result[parts[1]] = parts[2]
    end
    return result
end

verify_container_image(policy::ContainerPolicy, object) = verify_container_image(policy.profile,object)
function verify_container_image(profile::ProfileDefinition, object)
    container_require(container_json(object),:image_inspection_unverified)
    digests, config = get(object,:RepoDigests,nothing), get(object,:Config,nothing)
    container_require(container_vector(digests) && profile.environment in digests,:image_digest_unverified)
    id = get(object,:Id,nothing)
    container_require(id isa AbstractString && occursin(r"^(?:sha256:)?[a-f0-9]{64}$",id),:image_identity_unverified)
    # Podman 4.x omits the algorithm prefix on its full content ID. The exact
    # approved RepoDigest above remains mandatory; tags and short IDs fail.
    startswith(id,"sha256:") || (id = "sha256:" * id)
    container_require(get(object,:Os,nothing) == "linux" && container_json(config),:image_configuration_unverified)
    labels = get(config,:Labels,nothing)
    container_require(container_json(labels) && get(labels,Symbol("org.linecablemodels.runtime-image"),nothing) == "1" &&
        get(labels,Symbol("org.linecablemodels.profile"),nothing) == profile.id &&
        get(labels,Symbol("org.linecablemodels.profile-kind"),nothing) == string(profile.kind),:image_contract_unverified)
    for key in (:Volumes,:OnBuild,:Healthcheck)
        container_require(container_empty(get(config,key,nothing)),:image_side_effects_unapproved)
    end
    environment = container_env_map(get(config,:Env,nothing))
    # The official Julia base image carries version/path metadata. Application
    # environment must be explicit; wildcard/proxy/credential inheritance fails.
    allowed = union(keys(container_environment(profile)),("JULIA_PATH","JULIA_VERSION","JULIA_GPG"))
    container_require(all(in(allowed),keys(environment)),:image_environment_unapproved)
    return (; id=String(id), environment)
end

function verify_created_container(policy::ContainerPolicy, object, image; running=false)
    container_require(container_json(object),:container_inspection_unverified)
    config, host, state = get(object,:Config,nothing),get(object,:HostConfig,nothing),get(object,:State,nothing)
    container_require(all(container_json,(config,host,state)),:container_inspection_unverified)
    name = get(object,:Name,""); startswith(name,"/") && (name = name[2:end])
    container_require(matches_resource(policy.receipt,policy.receipt.scope,get(object,:Id,nothing),name,
        get(config,:Labels,nothing)),:container_ownership_unverified)
    observed_image = get(object,:Image,"")
    container_require(observed_image == image.id || "sha256:" * observed_image == image.id,:container_image_changed)
    container_require(get(state,:Running,nothing) === running && (running || get(state,:Pid,nothing) === 0),:container_state_unverified)
    for (key,value) in ((:User,"1000:1000"),(:WorkingDir,"/tmp"),(:Hostname,"lcm-executor"),(:OpenStdin,true),
            (:Tty,policy.profile.kind == :terminal),(:Cmd,container_julia_arguments(policy)))
        container_require(container_same(get(config,key,nothing),value),:container_entry_unverified)
    end
    # Podman 4.x reports a single entrypoint as a string; newer versions and
    # Docker use an argv array. Accept only the same exact fixed executable,
    # never shell splitting or a command string containing extra arguments.
    entry = get(config,:Entrypoint,nothing)
    container_require(entry == ["/usr/local/julia/bin/julia"] ||
        (policy.receipt.backend == :podman && entry == "/usr/local/julia/bin/julia"),:container_entry_unverified)
    expected_env = merge(image.environment,container_environment(policy))
    env = container_env_map(get(config,:Env,nothing))
    # Podman adds only this fixed engine marker; no inherited proxy/host secrets.
    if haskey(env,"container")
        container_require(policy.receipt.backend == :podman && pop!(env,"container") == "podman",:environment_unverified)
    end
    container_require(env == expected_env,:container_environment_changed)
    for key in (:Volumes,:Healthcheck)
        container_require(container_empty(get(config,key,nothing)),:container_mounts_unapproved)
    end
    b = policy.limits
    for (key,value) in ((:ReadonlyRootfs,true),(:Privileged,false),(:NetworkMode,"none"),(:IpcMode,"private"),
            (policy.receipt.backend == :podman ? :CgroupMode : :CgroupnsMode,"private"),
            (:Memory,b.memory_bytes),(:MemorySwap,b.memory_bytes),(:PidsLimit,b.pids),
            (:CpuPeriod,ExecutionCore.CONTAINER_CPU_PERIOD),(:CpuQuota,floor(Int,b.cpus*ExecutionCore.CONTAINER_CPU_PERIOD)),
            (:ShmSize,ExecutionCore.CONTAINER_SHM_BYTES))
        container_require(container_same(get(host,key,nothing),value),:container_limits_changed)
    end
    for key in (:PidMode,:UTSMode)
        container_require(get(host,key,nothing) in (policy.receipt.backend == :podman ? ("private","") : ("",)),:container_namespace_unverified)
    end
    for key in (:Binds,:VolumesFrom,:Devices,:DeviceRequests,:DeviceCgroupRules,:CapAdd,:GroupAdd,:PortBindings,:Links)
        container_require(container_empty(get(host,key,nothing)),:container_privileges_unapproved)
    end
    drops = get(host,:CapDrop,nothing)
    if policy.receipt.backend == :podman
        # Podman expands ALL into its configured defaults and exposes the
        # resulting OCI capability sets separately (JSON null means empty).
        container_require(container_vector(drops) && 1 <= length(drops) <= 64 && allunique(drops) &&
            all(c->occursin(r"^CAP_[A-Z_]+$",c),drops) &&
            all(k->haskey(object,k) && container_empty(object[k]),(:EffectiveCaps,:BoundingCaps)),:container_capabilities_unverified)
    else
        container_require(container_vector(drops) && lowercase.(drops) == ["all"],:container_capabilities_unverified)
    end
    security = get(host,:SecurityOpt,nothing)
    container_require(container_vector(security) && length(security) == 1 &&
        only(security) in ("no-new-privileges","no-new-privileges=true"),:container_security_unverified)
    logs,restart = get(host,:LogConfig,nothing),get(host,:RestartPolicy,nothing)
    container_require(container_json(logs) && get(logs,:Type,nothing) == "none" &&
        container_json(restart) && get(restart,:Name,nothing) == "no",:container_background_work_unapproved)
    tmpfs = get(host,:Tmpfs,nothing)
    container_require(container_json(tmpfs) && length(tmpfs) == 1,:container_scratch_unverified)
    options = get(tmpfs,Symbol("/tmp"),nothing)
    container_require(options isa AbstractString,:container_scratch_unverified)
    items = split(options,',')
    expected = Set(split(container_tmpfs_options(policy),','))
    if policy.receipt.backend == :podman
        # Podman consumes notmpcopyup during specification generation; inspect
        # omits it and, crucially, no longer reports its default tmpcopyup option.
        "notmpcopyup" in items || delete!(expected,"notmpcopyup")
        "rprivate" in items && push!(expected,"rprivate")
    end
    container_require(allunique(items) && Set(items) == expected,:container_scratch_unverified)
    mounts = get(object,:Mounts,nothing)
    container_require(mounts isa AbstractVector && length(mounts) <= 1 && all(m -> container_json(m) &&
        get(m,:Type,nothing) == "tmpfs" && get(m,:Destination,nothing) == "/tmp",mounts),:container_mounts_unapproved)
    ulimits = get(host,:Ulimits,nothing)
    container_require(ulimits isa AbstractVector,:container_rlimits_unverified)
    if policy.receipt.backend == :podman && length(ulimits) == 4
        # Podman materializes a host-derived finite NPROC ceiling at start.
        # It cannot relax the separately verified per-container cgroup PID cap.
        extra = filter(u->container_json(u) && get(u,:Name,nothing) == "RLIMIT_NPROC",ulimits)
        container_require(length(extra) == 1,:container_rlimits_unverified)
        soft,hard = get(only(extra),:Soft,nothing),get(only(extra),:Hard,nothing)
        container_require(all(v->v isa Integer && !(v isa Bool),(soft,hard)) &&
            0 < soft <= hard <= typemax(Int32),:container_rlimits_unverified)
    else
        container_require(length(ulimits) == 3,:container_rlimits_unverified)
    end
    for (key,limit) in (("core",0),("msgqueue",0),("nofile",1024))
        name = policy.receipt.backend == :podman ? "RLIMIT_" * uppercase(key) : key
        matches = filter(u->container_json(u) && get(u,:Name,nothing) == name,ulimits)
        container_require(length(matches) == 1 && get(only(matches),:Soft,nothing) == limit &&
            get(only(matches),:Hard,nothing) == limit,:container_rlimits_unverified)
    end
    return nothing
end

function container_inspection(runner,host,arguments;invoke)
    result = scoped_container_command(runner,host,arguments;invoke)
    container_require(result.exitcode == 0 && ncodeunits(result.output) <= 1024^2,:container_inspection_unavailable)
    object = try JSON3.read(result.output) catch; nothing end
    container_require(container_json(object),:container_inspection_unverified)
    return object
end

"""
    create_owned_container!(journal, runner, host, profile, fence) -> ResourceReceipt

Verify host prerequisites, record intent, verify the already-installed approved
image, create a stopped container, and bind its full ID after checking its policy.
Never pull or start an image. Failure after intent retains its receipt so the
physical owner can retire partial acquisition before releasing capacity. The
entry guard must subsequently verify actual kernel limits inside the started
process; this function is not a readiness or effective-isolation attestation.
"""
function create_owned_container!(journal::ResourceJournal,runner::CommandRunner,host::ContainerHostCheck,
        profile::ProfileDefinition,fence::AssignmentFence;
        probe=arguments->container_probe(runner,arguments),
        invoke=arguments->run_owned_command!(runner,setenv(Cmd(arguments),container_command_environment())),
        scope_options...)
    container_require(isempty(host.failures),:host_prerequisites_unavailable)
    scope = container_scope(runner,host;probe,scope_options...)
    receipt = reserve_resource!(journal,fence,host.engine.name,scope)
    policy = ContainerPolicy(profile,receipt)
    image = verify_container_image(policy,container_inspection(runner,host,
        ["image","inspect","--format","{{json .}}",profile.environment];invoke))
    container_require(container_scope(runner,host;probe,scope_options...) == scope,:container_scope_changed)
    # An existing intent is not permission to repeat creation or restart a process.
    container_require(receipt.physical_id === nothing && container_absent(runner,host,receipt;invoke),:container_already_acquired)
    result = scoped_container_command(runner,host,container_create_arguments(policy);invoke)
    container_require(result.exitcode == 0,:container_create_failed)
    id = strip(result.output)
    container_require(occursin(r"^[a-f0-9]{64}$",id),:container_identity_unverified)
    object = container_inspection(runner,host,["container","inspect","--format","{{json .}}",id];invoke)
    container_require(container_scope(runner,host;probe,scope_options...) == scope,:container_scope_changed)
    verify_created_container(policy,object,image)
    return bind_resource!(journal,receipt,id)
end
