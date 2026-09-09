"""
    NativeIdentity(uid, gid, cgroup)

Declare the unprivileged identity and exact generated user-service control group
expected by a native scientific entry guard. Native execution is trusted code,
not a private terminal or an arbitrary-code sandbox. This value conveys no lease
authority or preparation state.
"""
struct NativeIdentity
    "Expected non-root user ID."
    uid::Int
    "Expected non-root primary group ID."
    gid::Int
    "Absolute generated executor group in the host cgroup-v2 hierarchy."
    cgroup::String
    function NativeIdentity(uid, gid, cgroup)
        all(v->v isa Integer && !(v isa Bool) && 0 < v < typemax(UInt32), (uid,gid)) ||
            throw(ArgumentError("native execution requires a non-root user and group"))
        prefix = "/user.slice/user-$uid.slice/user@$uid.service/app.slice/"
        cgroup isa AbstractString && startswith(cgroup,prefix) &&
            occursin(r"^lcm-exec-[a-f0-9-]{36}\.service$",chopprefix(cgroup,prefix)) ||
            throw(ArgumentError("native executor control group is invalid"))
        token = chopprefix(chopsuffix(cgroup,".service"),prefix * "lcm-exec-")
        try UUID(token) catch; throw(ArgumentError("native executor control group is invalid")) end
        new(uid,gid,String(cgroup))
    end
end
Base.show(io::IO, ::NativeIdentity) = print(io,"NativeIdentity(<private>)")
Base.show(io::IO, ::MIME"text/plain", value::NativeIdentity) = show(io,value)

function isolation_native_process(status, identity::NativeIdentity)
    for (key,id) in (("Uid",identity.uid),("Gid",identity.gid))
        values = split(get(status,key,""))
        length(values)==4 && all(==(string(id)),values) ||
            throw(IsolationError(:unprivileged_identity_unverified))
    end
    for key in ("CapInh","CapPrm","CapEff","CapBnd","CapAmb")
        tryparse(UInt64,get(status,key,"");base=16) === UInt64(0) ||
            throw(IsolationError(:capabilities_not_dropped))
    end
    get(status,"NoNewPrivs","") == "1" || throw(IsolationError(:new_privileges_not_blocked))
    return nothing
end

function isolation_native_mounts(mounts, limits::ExecutorLimits, disk_stat)
    bypath = Dict(m.path=>m for m in mounts)
    for path in ("/","/proc","/sys","/sys/fs/cgroup","/dev","/tmp","/dev/shm")
        haskey(bypath,path) || throw(IsolationError(:required_mount_missing))
    end
    bypath["/proc"].filesystem == "proc" || throw(IsolationError(:proc_mount_unverified))
    bypath["/sys/fs/cgroup"].filesystem == "cgroup2" || throw(IsolationError(:cgroup_mount_unverified))
    for path in ("/","/sys","/sys/fs/cgroup","/dev")
        bypath[path].readonly || throw(IsolationError(:native_filesystem_writable))
    end
    total = 0
    for mount in mounts
        path = mount.path
        if path in ("/tmp","/dev/shm")
            mount.filesystem == "tmpfs" && !mount.readonly &&
                all(in(mount.options),("nosuid","nodev","noexec")) ||
                throw(IsolationError(:scratch_mount_options_unverified))
            size = isolation_scratch_size(path,disk_stat)
            size <= limits.scratch_bytes - total || throw(IsolationError(:scratch_limit_exceeded))
            path == "/dev/shm" && size > CONTAINER_SHM_BYTES &&
                throw(IsolationError(:shared_memory_limit_exceeded))
            total += size
        elseif !mount.readonly
            # Proc and devpts are kernel interfaces, not writable disk scratch.
            # Native profiles remain trusted: this is not a host PID sandbox.
            ((path == "/proc" || startswith(path,"/proc/")) && mount.filesystem == "proc") ||
                (path == "/dev/pts" && mount.filesystem == "devpts") ||
                throw(IsolationError(:native_filesystem_writable))
        end
    end
    return total
end

function check_native_isolation(limits::ExecutorLimits, identity::NativeIdentity, read_file, disk_stat)
    isolation_native_process(isolation_status(read_file("/proc/self/status")),identity)
    strip(read_file("/proc/self/cgroup")) == "0::" * identity.cgroup ||
        throw(IsolationError(:native_cgroup_identity_unverified))
    disk_stat("/sys/fs/cgroup").ftype == 0x63677270 || throw(IsolationError(:cgroup_v2_required))
    observed = isolation_cgroup_limits(limits,"/sys/fs/cgroup" * identity.cgroup,read_file)
    isolation_rlimits(read_file("/proc/self/limits"))
    isolation_network(read_file)
    scratch = isolation_native_mounts(isolation_mounts(read_file("/proc/self/mountinfo")),limits,disk_stat)
    return (;schema_version=1,observed...,scratch_bytes=scratch,uid=identity.uid)
end

"""
    verify_native_isolation(limits::ExecutorLimits, identity::NativeIdentity)

Check effective Linux kernel limits before importing a trusted scientific profile.
Require the exact non-root service identity, finite CPU/memory/task quotas, zero
additional swap, read-only ordinary filesystems, bounded private tmpfs scratch,
no effective capabilities, no new privileges, isolated networking and bounded
descriptor/message-queue/core-file limits. Native profiles retain host visibility;
this does not admit arbitrary terminal code or establish readiness.

Return observed limits on success. Missing, excessive or unsupported evidence
throws an IsolationError with a fixed non-sensitive code. Service ownership,
source verification, lease expiry and physical cleanup remain agent duties.
"""
function verify_native_isolation(limits::ExecutorLimits, identity::NativeIdentity)
    failure, report = nothing,nothing
    try
        Sys.islinux() || throw(IsolationError(:linux_required))
        report = check_native_isolation(limits,identity,isolation_kernel_text,Base.Filesystem.diskstat)
    catch error
        failure = error isa IsolationError ? error.code : :kernel_evidence_unavailable
    end
    failure === nothing || throw(IsolationError(failure))
    return report
end
