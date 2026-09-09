"""
    ExecutorLimits(cpus, memory_bytes, pids, scratch_bytes)

Describe the mandatory kernel limits for one approved executor. CPU is measured
in logical CPU units; memory and writable scratch are bytes. PID count includes
threads. Runtime operation deadlines and lease authority remain separate.
"""
struct ExecutorLimits
    "Maximum CPU quota in logical CPU units."
    cpus::Float64
    "Maximum cgroup memory in bytes, with no additional swap."
    memory_bytes::Int
    "Maximum cgroup tasks, including threads."
    pids::Int
    "Combined maximum writable tmpfs capacity in bytes."
    scratch_bytes::Int
    function ExecutorLimits(cpus, memory_bytes, pids, scratch_bytes)
        cpus isa Real && !(cpus isa Bool) && isfinite(cpus) && 0.01 <= cpus <= 256 ||
            throw(ArgumentError("executor CPU quota must be in [0.01, 256]"))
        all(v -> v isa Integer && !(v isa Bool) && 0 < v <= typemax(Int),
            (memory_bytes, pids, scratch_bytes)) || throw(ArgumentError("executor byte and task limits must be positive integers"))
        pids <= 65536 || throw(ArgumentError("executor task limit exceeds 65536"))
        scratch_bytes >= CONTAINER_SHM_BYTES + 4096 ||
            throw(ArgumentError("executor scratch must reserve shared memory and at least one page"))
        return new(cpus, memory_bytes, pids, scratch_bytes)
    end
end

"""Compatibility name for the shared native/container executor limits."""
const ContainerLimits = ExecutorLimits

"""Fixed unprivileged container UID and GID for approved runtime images."""
const CONTAINER_USER_ID = 1000
"""Shared-memory tmpfs allowance, in bytes, included in the scratch budget."""
const CONTAINER_SHM_BYTES = 65536
"""CPU quota period in microseconds, shared by both container adapters."""
const CONTAINER_CPU_PERIOD = 100000

"""Report only a fixed isolation failure code, without private kernel paths/data."""
struct IsolationError <: Exception
    "Finite prerequisite or effective-policy failure."
    code::Symbol
end
Base.showerror(io::IO, error::IsolationError) = print(io, "Executor isolation check failed: ", error.code)

function isolation_kernel_text(path::String)
    maximum = path == "/proc/self/mountinfo" ? 256 * 1024 : 64 * 1024
    return open(path, "r") do io
        bytes = read(io, maximum + 1)
        length(bytes) <= maximum || throw(IsolationError(:kernel_evidence_oversized))
        String(bytes)
    end
end

function isolation_integer(text, code::Symbol; zero=false)
    value = tryparse(Int, strip(text))
    value !== nothing && (zero ? value >= 0 : value > 0) || throw(IsolationError(code))
    return value
end

function isolation_status(text::AbstractString)
    ncodeunits(text) <= 64 * 1024 || throw(IsolationError(:kernel_evidence_oversized))
    values = Dict{String,String}()
    for line in split(text, '\n'; keepempty=false)
        parts = split(line, ':'; limit=2)
        length(parts) == 2 || throw(IsolationError(:process_status_unverified))
        key = strip(parts[1])
        haskey(values, key) && throw(IsolationError(:process_status_unverified))
        values[key] = strip(parts[2])
    end
    return values
end

function isolation_process(status)
    for key in ("Uid", "Gid")
        values = split(get(status, key, ""))
        length(values) == 4 && all(==(string(CONTAINER_USER_ID)), values) ||
            throw(IsolationError(:unprivileged_identity_unverified))
    end
    haskey(status,"Groups") || throw(IsolationError(:supplementary_groups_unverified))
    groups = split(status["Groups"])
    all(==(string(CONTAINER_USER_ID)), groups) || throw(IsolationError(:supplementary_groups_unverified))
    for key in ("CapInh", "CapPrm", "CapEff", "CapBnd", "CapAmb")
        value = tryparse(UInt64, get(status, key, ""); base=16)
        value === UInt64(0) || throw(IsolationError(:capabilities_not_dropped))
    end
    get(status, "NoNewPrivs", "") == "1" || throw(IsolationError(:new_privileges_not_blocked))
    get(status, "Seccomp", "") == "2" || throw(IsolationError(:seccomp_not_enforced))
    return nothing
end

function isolation_mounts(text)
    ncodeunits(text) <= 256 * 1024 || throw(IsolationError(:kernel_evidence_oversized))
    lines = split(strip(text), '\n'; keepempty=false)
    1 <= length(lines) <= 512 || throw(IsolationError(:mount_inventory_unverified))
    records = NamedTuple[]
    seen = Set{String}()
    for line in lines
        sides = split(line, " - "; limit=2)
        length(sides) == 2 || throw(IsolationError(:mount_inventory_unverified))
        fields, tail = split(sides[1]), split(sides[2])
        length(fields) >= 6 && length(tail) == 3 || throw(IsolationError(:mount_inventory_unverified))
        path = String(fields[5])
        # Approved mount points are fixed paths; encoded/control-containing
        # alternatives are not silently normalized into an allowed destination.
        startswith(path, "/") && !occursin(r"[\\\x00-\x20]", path) ||
            throw(IsolationError(:mount_inventory_unverified))
        path in seen && throw(IsolationError(:stacked_mounts_unverified))
        push!(seen, path)
        options = Set(String.(split(fields[6], ',')))
        super_options = Set(String.(split(tail[3], ',')))
        readonly = "ro" in options || "ro" in super_options
        push!(records, (; path, filesystem=String(tail[1]), options, readonly))
    end
    return records
end

function isolation_root_owned(path, file_stat; directory=false, device=false)
    metadata = file_stat(path)
    metadata.uid == 0 && (device || metadata.mode & 0o022 == 0) ||
        throw(IsolationError(:writable_runtime_mount))
    if directory
        isdir(metadata) || throw(IsolationError(:mount_type_unverified))
    elseif device
        ischardev(metadata) || throw(IsolationError(:mount_type_unverified))
    else
        isfile(metadata) || throw(IsolationError(:mount_type_unverified))
    end
    return nothing
end

function isolation_scratch_size(path, disk_stat)
    stats = disk_stat(path)
    stats.ftype == 0x01021994 || throw(IsolationError(:scratch_not_tmpfs))
    size = Int128(stats.bsize) * Int128(stats.blocks)
    0 < size <= typemax(Int) || throw(IsolationError(:scratch_limit_unverified))
    return Int(size)
end

function isolation_standard_device(path, number, mount, file_stat)
    metadata = file_stat(path)
    # OCI permits the six default character devices to be bind-mounted.
    # Rootless runc exposes their unmapped host-root owner as 65534. This
    # exception never applies to ordinary files, directories or other devices.
    mount.filesystem in ("tmpfs", "devtmpfs") && "nosuid" in mount.options &&
        ischardev(metadata) && metadata.rdev == number &&
        (metadata.uid == 0 || (mount.filesystem == "devtmpfs" && metadata.uid == 65534)) ||
        throw(IsolationError(:device_mount_unverified))
    return nothing
end

function isolation_mount_policy(mounts, limits::ContainerLimits, file_stat, disk_stat)
    bypath = Dict(m.path=>m for m in mounts)
    for path in ("/", "/proc", "/sys", "/sys/fs/cgroup", "/tmp", "/dev", "/dev/pts")
        haskey(bypath, path) || throw(IsolationError(:required_mount_missing))
    end
    bypath["/"].readonly || throw(IsolationError(:root_filesystem_writable))
    bypath["/sys"].readonly || throw(IsolationError(:sys_filesystem_writable))
    bypath["/sys/fs/cgroup"].filesystem == "cgroup2" && bypath["/sys/fs/cgroup"].readonly ||
        throw(IsolationError(:cgroup_mount_not_private_readonly))
    bypath["/proc"].filesystem == "proc" || throw(IsolationError(:proc_mount_unverified))
    total = 0
    for mount in mounts
        path, fs = mount.path, mount.filesystem
        if path in ("/tmp", "/dev/shm")
            fs == "tmpfs" || throw(IsolationError(:scratch_not_tmpfs))
            if !mount.readonly
                all(in(mount.options), ("nosuid", "nodev", "noexec")) ||
                    throw(IsolationError(:scratch_mount_options_unverified))
                size = isolation_scratch_size(path, disk_stat)
                size <= limits.scratch_bytes - total || throw(IsolationError(:scratch_limit_exceeded))
                path == "/dev/shm" && size > CONTAINER_SHM_BYTES &&
                    throw(IsolationError(:shared_memory_limit_exceeded))
                total += size
            end
        elseif path == "/dev"
            fs == "tmpfs" || throw(IsolationError(:device_mount_unverified))
            isolation_root_owned(path, file_stat; directory=true)
        elseif path == "/dev/pts"
            fs == "devpts" || throw(IsolationError(:terminal_mount_unverified))
        elseif path in ("/dev/null", "/dev/zero", "/dev/full", "/dev/random", "/dev/urandom", "/dev/tty")
            # Fixed Linux major/minor pairs from the OCI default-device set.
            numbers = ("/dev/null"=>0x0103, "/dev/zero"=>0x0105, "/dev/full"=>0x0107,
                "/dev/random"=>0x0108, "/dev/urandom"=>0x0109, "/dev/tty"=>0x0500)
            number = last(only(filter(pair->first(pair) == path, numbers)))
            isolation_standard_device(path, number, mount, file_stat)
        elseif path in ("/dev/console", "/dev/ptmx")
            fs == "devpts" && ischardev(file_stat(path)) ||
                throw(IsolationError(:terminal_mount_unverified))
        elseif path == "/dev/mqueue"
            fs == "mqueue" || throw(IsolationError(:message_queue_mount_unverified))
        elseif path in ("/etc/hosts", "/etc/hostname", "/etc/resolv.conf")
            isolation_root_owned(path, file_stat)
        elseif path in ("/etc/passwd", "/etc/group", "/run/.containerenv")
            # Podman supplies container-local identity metadata for the fixed
            # numeric user. No writable file or entire /run mount is allowed.
            mount.readonly || throw(IsolationError(:writable_runtime_mount))
            isolation_root_owned(path, file_stat)
        elseif path in ("/", "/proc", "/sys", "/sys/fs/cgroup")
            nothing
        elseif startswith(path, "/proc/") || startswith(path, "/sys/")
            if !mount.readonly
                # OCI masks selected proc pseudo-files with /dev/null. No writable
                # ordinary files or unbounded tmpfs are accepted through this path.
                fs in ("tmpfs", "devtmpfs") || throw(IsolationError(:writable_kernel_mount))
                isolation_standard_device(path, 0x0103, mount, file_stat)
                ischardev(file_stat("/dev/null")) && file_stat("/dev/null").rdev == 0x0103 ||
                    throw(IsolationError(:kernel_mask_unverified))
            end
        else
            throw(IsolationError(:unapproved_mount))
        end
    end
    !bypath["/tmp"].readonly && total > 0 || throw(IsolationError(:writable_scratch_missing))
    return total
end

function isolation_rlimits(text)
    ncodeunits(text) <= 64 * 1024 || throw(IsolationError(:kernel_evidence_oversized))
    lines = split(text, '\n')
    for (label, maximum) in (("Max open files", 1024), ("Max msgqueue size", 0), ("Max core file size", 0))
        matches = filter(line -> startswith(line, label * " "), lines)
        length(matches) == 1 || throw(IsolationError(:process_limits_unverified))
        parts = split(strip(only(matches)[length(label)+1:end]))
        length(parts) >= 2 || throw(IsolationError(:process_limits_unverified))
        soft = isolation_integer(parts[1], :process_limits_unverified; zero=true)
        hard = isolation_integer(parts[2], :process_limits_unverified; zero=true)
        soft <= hard <= maximum || throw(IsolationError(:process_limits_unverified))
    end
    return nothing
end

function isolation_cgroup_limits(limits::ExecutorLimits, root::String, read_file)
    quota = split(strip(read_file(root * "/cpu.max")))
    length(quota) == 2 || throw(IsolationError(:cpu_limit_unverified))
    cpu_quota = isolation_integer(quota[1], :cpu_limit_unverified)
    cpu_period = isolation_integer(quota[2], :cpu_limit_unverified)
    1000 <= cpu_period <= 1000000 || throw(IsolationError(:cpu_limit_unverified))
    cpu_quota <= floor(Int, limits.cpus * cpu_period) || throw(IsolationError(:cpu_limit_exceeded))
    memory = isolation_integer(read_file(root * "/memory.max"), :memory_limit_unverified)
    memory <= limits.memory_bytes || throw(IsolationError(:memory_limit_exceeded))
    swap = isolation_integer(read_file(root * "/memory.swap.max"), :swap_limit_unverified; zero=true)
    swap == 0 || throw(IsolationError(:additional_swap_allowed))
    pids = isolation_integer(read_file(root * "/pids.max"), :pids_limit_unverified)
    pids <= limits.pids || throw(IsolationError(:pids_limit_exceeded))
    return (; cpus=cpu_quota / cpu_period, memory_bytes=memory, pids)
end

function isolation_network(read_file)
    devices = split(read_file("/proc/net/dev"), '\n')
    interfaces = String[]
    for line in devices
        occursin(':', line) || continue
        push!(interfaces, strip(first(split(line, ':'; limit=2))))
    end
    interfaces == ["lo"] || throw(IsolationError(:network_not_isolated))
    return nothing
end

function check_container_isolation(limits::ContainerLimits, read_file, file_stat, disk_stat)
    # Both the current Julia process and namespace init must run under the
    # approved unprivileged identity, capability set and seccomp policy.
    isolation_process(isolation_status(read_file("/proc/self/status")))
    isolation_process(isolation_status(read_file("/proc/1/status")))
    strip(read_file("/proc/self/cgroup")) == "0::/" ||
        throw(IsolationError(:private_cgroup_namespace_unverified))
    observed = isolation_cgroup_limits(limits, "/sys/fs/cgroup", read_file)
    isolation_network(read_file)
    isolation_rlimits(read_file("/proc/self/limits"))
    mounts = isolation_mounts(read_file("/proc/self/mountinfo"))
    scratch = isolation_mount_policy(mounts, limits, file_stat, disk_stat)
    return (; schema_version=1, observed..., scratch_bytes=scratch, uid=CONTAINER_USER_ID)
end

"""
    verify_container_isolation(limits::ContainerLimits)

Verify effective kernel evidence inside an approved Linux container before
including a scientific profile or entering the Julia REPL. Require finite CPU,
memory, zero additional swap, task and writable-tmpfs limits; unprivileged UID/GID;
dropped capabilities; no-new-privileges; seccomp; private read-only cgroups; isolated
networking; bounded descriptor/message-queue/core-file limits; and only approved
mount locations. Return a bounded report of observed limits, not preparation or
lease readiness. Missing/unsupported evidence fails closed with an IsolationError.

The agent must also verify the created container's image, namespaces, configuration
and ownership. This check does not authorize arbitrary images, make untrusted code
host-safe, or replace the agent's lease expiry and physical cleanup responsibility.
"""
function verify_container_isolation(limits::ContainerLimits)
    failure = nothing
    report = nothing
    try
        Sys.islinux() || throw(IsolationError(:linux_required))
        report = check_container_isolation(limits, isolation_kernel_text, stat, Base.Filesystem.diskstat)
    catch error
        failure = error isa IsolationError ? error.code : :kernel_evidence_unavailable
    end
    failure === nothing || throw(IsolationError(failure))
    return report
end
