"""
    NativeUnitIdentity

Retain one observed, receipt-matched native service state. This is cleanup
identity, not a lease, resource-limit attestation or preparation report. A verified
never-started service has no invocation. The PID
is diagnostic only; cleanup never signals a saved PID.
"""
struct NativeUnitIdentity
    "Generated resource service name."
    unit::String
    "Systemd invocation identity; nothing only for a verified never-started unit."
    invocation::Union{Nothing,String}
    "Expected kernel control group under the local user manager."
    cgroup::String
    "Current service main PID, not a signaling target."
    pid::Int
    "Observed unit state."
    state::Symbol
end
Base.show(io::IO, value::NativeUnitIdentity) = print(io,"NativeUnitIdentity(",value.unit,", <owned>)")
Base.show(io::IO, ::MIME"text/plain", value::NativeUnitIdentity) = show(io,value)

native_require(condition) = condition === true ? nothing : throw(CommandFailure(:native_identity_unverified))
native_uid() = Int(ccall(:geteuid,Cuint,()))
native_bus_address(uid) = "unix:path=/run/user/$uid/bus"
native_cgroup(receipt,uid) = "/user.slice/user-$uid.slice/user@$uid.service/app.slice/" * resource_name(receipt)

function native_bus_command(runner,arguments;uid=native_uid(),options...)
    systemd_bus_command(runner,arguments;address=native_bus_address(uid),options...)
end

function native_bus_result(runner,arguments,type;options...)
    records = native_bus_command(runner,arguments;options...)
    native_require(length(records) == 1 && only(records).type == type &&
        only(records).data isa AbstractVector && length(only(records).data) == 1)
    return only(only(records).data)
end

function native_local_bus(uid)
    path = "/run/user/$uid/bus"
    native_require(realpath(dirname(path)) == dirname(path))
    metadata = lstat(path)
    native_require(metadata.uid == uid && metadata.mode & 0o170000 == 0o140000)
    return nothing
end

"""
    native_scope(runner) -> String

Verify a local user-manager endpoint and hash its machine/user namespace. The
socket is an owned Unix socket at the fixed runtime path; inherited remote bus
addresses are ignored. Typed bus machine and service-owner identities must match
the local machine and effective UID. Missing controllers do not prevent cleanup.
The scope survives manager/reboot replacement; individual resources additionally
require their exact invocation identity and kernel-group absence.
"""
function native_scope(runner::CommandRunner;uid=native_uid(),
        machine_id=strip(read("/etc/machine-id",String)),socket_check=native_local_bus,options...)
    native_require(uid isa Integer && !(uid isa Bool) && uid > 0 &&
        machine_id isa AbstractString && occursin(r"^[a-f0-9]{32}$",machine_id))
    socket_check(uid) === nothing || throw(CommandFailure(:native_identity_unverified))
    bus_machine = native_bus_result(runner,["call","org.freedesktop.DBus","/org/freedesktop/DBus",
        "org.freedesktop.DBus.Peer","GetMachineId"],"s";uid,options...)
    owner = native_bus_result(runner,["call","org.freedesktop.DBus","/org/freedesktop/DBus",
        "org.freedesktop.DBus","GetConnectionUnixUser","s","org.freedesktop.systemd1"],"u";uid,options...)
    native_require(bus_machine == machine_id && owner isa Integer && !(owner isa Bool) && owner == uid)
    return bytes2hex(sha256(JSON3.write(("systemd-user-native-v1",machine_id,string(uid)))))
end

"""
    native_resource_description(receipt) -> String

Render the fixed native-unit ownership marker from every shared receipt label.
This marker is immutable launch metadata, not a browser label or credential.
"""
function native_resource_description(receipt::ResourceReceipt)
    native_require(receipt.backend == :native)
    labels = sort!(collect(resource_labels(receipt));by=first)
    return "LCM native executor v1 " * bytes2hex(sha256(JSON3.write(labels)))
end

function native_unit_listing(runner,receipt;options...)
    values = native_bus_result(runner,["call","org.freedesktop.systemd1","/org/freedesktop/systemd1",
        "org.freedesktop.systemd1.Manager","ListUnitsByPatterns","asas","0","1",resource_name(receipt)],
        "a(ssssssouso)";options...)
    native_require(values isa AbstractVector && length(values) <= 1)
    isempty(values) && return nothing
    item = only(values)
    native_require(item isa AbstractVector && length(item) == 10 && item[1] == resource_name(receipt) &&
        all(v->v isa String,item[[1,2,3,4,5,6,7,9,10]]) && item[3] == "loaded" &&
        occursin(r"^/org/freedesktop/systemd1/unit/[A-Za-z0-9_]+$",item[7]))
    return String(item[7])
end

function native_properties(runner,path,interface,properties;uid=native_uid(),options...)
    systemd_properties(runner,path,interface,properties;address=native_bus_address(uid),options...)
end

function native_invocation(properties)
    value = systemd_property(properties,"InvocationID","ay")
    native_require(value isa AbstractVector && length(value) in (0,16) &&
        all(v->v isa Integer && !(v isa Bool) && 0 <= v <= 255,value))
    return any(!iszero,value) ? bytes2hex(UInt8.(value)) : nothing
end

"""
    inspect_native_unit(runner, receipt) -> Union{Nothing,NativeUnitIdentity}

Inspect the exact generated unit without loading or starting it. Require transient
ownership metadata, invocation identity, fixed cleanup policy and the expected
user control group. Return nothing only after a successful empty unit inventory;
callers must also check pending jobs and the kernel group before proving absence.
"""
function inspect_native_unit(runner::CommandRunner,receipt::ResourceReceipt;uid=native_uid(),options...)
    native_require(receipt.backend == :native)
    path = native_unit_listing(runner,receipt;uid,options...)
    path === nothing && return nothing
    fields = ["Id","Description","InvocationID","Transient","ActiveState"]
    before = native_properties(runner,path,"Unit",fields;uid,options...)
    invocation = native_invocation(before)
    native_require(systemd_property(before,"Id","s") == resource_name(receipt) &&
        systemd_property(before,"Description","s") == native_resource_description(receipt) &&
        systemd_property(before,"Transient","b") === true &&
        (receipt.physical_id === nothing || receipt.physical_id == invocation))
    service = native_properties(runner,path,"Service",["MainPID","ControlGroup","Slice","Type","ExitType",
        "RemainAfterExit","Restart","KillMode","KillSignal","SendSIGKILL","FinalKillSignal",
        "TimeoutStopUSec","TimeoutStopFailureMode","ExecStop","ExecStopPost"];uid,options...)
    after = native_properties(runner,path,"Unit",fields;uid,options...)
    native_require(all(key->before[key].data == after[key].data,fields[1:4]))
    state = systemd_property(after,"ActiveState","s")
    pid = systemd_property(service,"MainPID","u")
    cgroup = systemd_property(service,"ControlGroup","s")
    native_require(state in ("active","activating","deactivating","inactive","failed") &&
        pid isa Integer && !(pid isa Bool) && 0 <= pid <= typemax(Int) &&
        (cgroup == native_cgroup(receipt,uid) || (isempty(cgroup) && pid == 0 && state in ("inactive","failed"))) &&
        systemd_property(service,"Slice","s") == "app.slice" &&
        systemd_property(service,"Type","s") == "exec" && systemd_property(service,"ExitType","s") == "main" &&
        systemd_property(service,"RemainAfterExit","b") === false && systemd_property(service,"Restart","s") == "no" &&
        systemd_property(service,"KillMode","s") == "control-group" && systemd_property(service,"KillSignal","i") == 2 &&
        systemd_property(service,"SendSIGKILL","b") === true && systemd_property(service,"FinalKillSignal","i") == 9 &&
        systemd_property(service,"TimeoutStopUSec","t") == 2_000_000 &&
        systemd_property(service,"TimeoutStopFailureMode","s") == "terminate")
    for key in ("ExecStop","ExecStopPost")
        commands = systemd_property(service,key,"a(sasbttttuii)")
        native_require(commands isa AbstractVector && isempty(commands))
    end
    # A transient unit can be queued before service_start assigns its invocation
    # ID. Only an unbound receipt, inactive state and zero process/group admit
    # this pre-activation identity. A bound receipt can never revert to it.
    invocation === nothing && native_require(receipt.physical_id === nothing &&
        state == "inactive" && pid == 0 && isempty(cgroup))
    return NativeUnitIdentity(resource_name(receipt),invocation,native_cgroup(receipt,uid),Int(pid),Symbol(state))
end

function native_jobs_clear(runner,receipt;options...)
    values = native_bus_result(runner,["call","org.freedesktop.systemd1","/org/freedesktop/systemd1",
        "org.freedesktop.systemd1.Manager","ListJobs"],"a(usssoo)";options...)
    native_require(values isa AbstractVector && length(values) <= 4096)
    for item in values
        native_require(item isa AbstractVector && length(item) == 6 && item[1] isa Integer &&
            !(item[1] isa Bool) && item[1] > 0 && all(v->v isa String && ncodeunits(v)<=1024,item[2:6]))
        item[2] == resource_name(receipt) && return false
    end
    return true
end

function native_group_empty(path)
    native_require(Base.Filesystem.diskstat("/sys/fs/cgroup").ftype == 0x63677270)
    directory = "/sys/fs/cgroup" * path
    metadata = lstat(directory)
    metadata.ioerrno == Base.UV_ENOENT && return true
    native_require(metadata.ioerrno == 0 && realpath(directory) == directory && isdir(metadata))
    content = open(joinpath(directory,"cgroup.events"),"r") do io; read(io,4097); end
    native_require(length(content) <= 4096)
    entries = [split(line) for line in split(String(content),'\n';keepempty=false)]
    native_require(all(pair->length(pair)==2,entries))
    populated = filter(pair->pair[1]=="populated",entries)
    native_require(length(populated)==1 && only(populated)[2] in ("0","1"))
    return only(populated)[2] == "0"
end

function native_service_command(runner,arguments;uid=native_uid(),which=Sys.which,
        invoke=command->run_owned_command!(runner,command))
    executable = which("systemctl")
    native_require(executable !== nothing)
    environment = container_command_environment()
    environment["DBUS_SESSION_BUS_ADDRESS"] = native_bus_address(uid)
    result = invoke(setenv(Cmd([executable,"--user","--no-pager",arguments...]),environment))
    native_require(result isa CommandResult)
    return result.exitcode == 0
end

"""
    remove_owned_native!(journal, runner, receipt; timeout_seconds=10) -> Bool

Retire only a receipt-matched transient native unit. Bind a started crash-gap intent
to its full invocation before requesting stop. A verified inactive, never-started
unit can have its queued start canceled while its receipt remains unbound.
Recheck scope and invocation before
every mutation; never signal saved PIDs or clear global jobs/failed services.

Require no pending unit job and an empty or removed exact kernel group before
forgetting the receipt. The unit must be absent or receipt-matched and stopped;
inactive metadata can remain cached while another unit references it.
Inspection errors, replacement invocations
and unfinished cleanup retain ownership. This operation neither launches work nor
requires CPU/memory controllers that may have disappeared since launch.
"""
function remove_owned_native!(journal::ResourceJournal,runner::CommandRunner,receipt::ResourceReceipt;
        timeout_seconds=10,uid=native_uid(),machine_id=strip(read("/etc/machine-id",String)),
        socket_check=native_local_bus,group_empty=native_group_empty,clock=()->time_ns()/1e9,
        pause=sleep,control=arguments->native_service_command(runner,arguments;uid),options...)
    native_require(receipt.backend == :native && timeout_seconds isa Real && !(timeout_seconds isa Bool) &&
        isfinite(timeout_seconds) && 0 < timeout_seconds <= 30)
    current = lock(journal.lock) do
        receipt.journal_id == journal.id || throw(ArgumentError("resource belongs to another journal"))
        any(r->r.id==receipt.id,resource_receipts(journal)) ? current_receipt(journal,receipt) : nothing
    end
    target = something(current,receipt)
    check_scope() = native_require(native_scope(runner;uid,machine_id,socket_check,options...) == target.scope)
    check_scope()
    deadline = clock() + timeout_seconds
    stop_requested = false
    reset_requested = false
    while clock() < deadline
        identity = inspect_native_unit(runner,target;uid,options...)
        if identity === nothing
            if native_jobs_clear(runner,target;uid,options...) && group_empty(native_cgroup(target,uid))
                check_scope()
                # Confirm the unit is still absent after the independent job,
                # kernel and scope observations; a late start retains ownership.
                if native_unit_listing(runner,target;uid,options...) === nothing
                    current === nothing || forget_resource!(journal,target)
                    return true
                end
            end
        elseif current !== nothing
            identity.invocation === nothing || (target = bind_resource!(journal,target,identity.invocation))
            check_scope()
            confirm = inspect_native_unit(runner,target;uid,options...)
            confirm === nothing && continue
            # If it started between observations, first persist that invocation
            # and repeat the identity check before any control command.
            identity.invocation === nothing && confirm.invocation !== nothing && continue
            if confirm.state in (:inactive,:failed) && group_empty(confirm.cgroup) &&
                    native_jobs_clear(runner,target;uid,options...)
                # Clear only this completed owned failed record; GC can then
                # remove it. A successful command is not itself absence proof.
                if confirm.state == :failed && !reset_requested
                    control(["reset-failed",resource_name(target)]) || return false
                    reset_requested = true
                elseif confirm.state == :inactive && !stop_requested && !reset_requested
                    control(["stop","--no-block",resource_name(target)]) || return false
                    stop_requested = true
                elseif confirm.state == :inactive && confirm.pid == 0
                    # A stopped transient definition may stay cached because
                    # another unit refers to it. Do not stop that other unit to
                    # force GC. Exact identity, no queued job and an empty kernel
                    # group prove physical retirement without erasing references.
                    check_scope()
                    final = inspect_native_unit(runner,target;uid,options...)
                    if final !== nothing && final.state == :inactive && final.pid == 0 &&
                            native_jobs_clear(runner,target;uid,options...) && group_empty(final.cgroup)
                        forget_resource!(journal,target)
                        return true
                    end
                end
            elseif !stop_requested
                control(["stop","--no-block",resource_name(target)]) || return false
                stop_requested = true
            end
        else
            # A removed receipt cannot authorize another mutation, but repeated
            # cleanup can confirm the same stopped cached definition read-only.
            if identity.state == :inactive && identity.pid == 0 && native_jobs_clear(runner,target;uid,options...) &&
                    group_empty(identity.cgroup)
                check_scope()
                final = inspect_native_unit(runner,target;uid,options...)
                return (final === nothing || (final.state == :inactive && final.pid == 0)) &&
                    native_jobs_clear(runner,target;uid,options...) && group_empty(identity.cgroup)
            end
            return false
        end
        clock() < deadline || return false
        pause(min(0.1,max(0.0,deadline-clock())))
    end
    return false
end

"""
    recover_native!(journal, runner)

Attempt every native receipt and retain unresolved ownership. Other backends are
untouched. The root agent still reconciles all backends before advertising any
profile; successful cleanup grants no live lease or preparation.
"""
function recover_native!(journal::ResourceJournal,runner::CommandRunner;options...)
    completed = true
    for receipt in filter(r->r.backend == :native,resource_receipts(journal))
        completed &= try remove_owned_native!(journal,runner,receipt;options...) catch; false end
    end
    completed || throw(ArgumentError("owned native recovery remains unresolved"))
    return nothing
end
