function local_container_storage_identity(path)
    root = realpath(path)
    metadata = stat(root)
    isdir(metadata) || throw(ArgumentError("container storage root is not a directory"))
    return (root, string(metadata.device), string(metadata.inode))
end

"""
    container_scope(runner, host) -> String

Inspect the exact local engine used by host.command and hash its identity with
this machine ID and effective UID. Docker uses its daemon ID; Podman uses its
persistent graph-root identity. A changed context/store cannot establish absence
of a resource recorded in another scope. This grants cleanup targeting only, not
permission to launch on a host with missing resource limits.
"""
function container_scope(runner::CommandRunner, host::ContainerHostCheck;
        probe=arguments -> container_probe(runner, arguments),
        machine_id=strip(read("/etc/machine-id", String)),
        storage_identity=local_container_storage_identity)
    occursin(r"^[a-f0-9]{32}$", machine_id) || throw(ArgumentError("local machine identity is unavailable"))
    isempty(host.command) && throw(ArgumentError("local container command is unavailable"))
    :local_engine_required in host.failures && throw(ArgumentError("container cleanup requires a local engine"))
    :linux_required in host.failures && throw(ArgumentError("container cleanup requires a Linux engine"))
    format = host.engine.name == :podman ? "json" : "{{json .}}"
    info = inspected_container_json(runner, [host.command; "info"; "--format"; format]; probe)
    return container_scope_identity(host,info;machine_id,storage_identity)
end

function container_scope_identity(host,info;
        machine_id=strip(read("/etc/machine-id",String)),storage_identity=local_container_storage_identity)
    occursin(r"^[a-f0-9]{32}$",machine_id) || throw(ArgumentError("local machine identity is unavailable"))
    reference = if host.engine.name == :docker
        id = get(info, :ID, nothing)
        id isa AbstractString && occursin(r"^[A-Za-z0-9:._-]{1,256}$", id) ||
            throw(ArgumentError("Docker daemon identity is unavailable"))
        get(info, :OSType, nothing) == "linux" || throw(ArgumentError("container engine changed operating system"))
        ("docker", String(id))
    elseif host.engine.name == :podman
        details = get(info, :host, nothing)
        store = get(info, :store, nothing)
        details isa JSON3.Object && get(details, :serviceIsRemote, nothing) === false &&
            get(details, :os, nothing) == "linux" || throw(ArgumentError("Podman engine is not local Linux"))
        store isa JSON3.Object || throw(ArgumentError("Podman storage identity is unavailable"))
        path = get(store, :graphRoot, nothing)
        path isa String && isabspath(path) && ncodeunits(path) <= 4096 && !occursin('\0', path) ||
            throw(ArgumentError("Podman storage identity is unavailable"))
        ("podman", storage_identity(path)...)
    else
        throw(ArgumentError("unsupported container engine"))
    end
    return bytes2hex(sha256(JSON3.write((machine_id, string(ccall(:geteuid, Cuint, ())), reference))))
end

function scoped_container_command(runner, host, arguments; invoke=arguments ->
        run_owned_command!(runner, setenv(Cmd(arguments), container_command_environment())))
    result = invoke([host.command; arguments])
    result isa CommandResult || throw(ArgumentError("container command returned an invalid result"))
    return result
end

function inspect_owned_container(runner, host, receipt, scope; invoke)
    target = something(receipt.physical_id, resource_name(receipt))
    result = scoped_container_command(runner, host,
        ["container", "inspect", "--format", "{{json .}}", target]; invoke)
    result.exitcode == 0 || return nothing
    ncodeunits(result.output) <= 1024^2 || throw(ArgumentError("container inspection exceeds its byte limit"))
    object = try JSON3.read(result.output) catch; nothing end
    object isa JSON3.Object || throw(ArgumentError("container inspection has an unsupported schema"))
    id = get(object, :Id, nothing)
    name = get(object, :Name, nothing)
    config = get(object, :Config, nothing)
    config isa JSON3.Object || throw(ArgumentError("container identity is unavailable"))
    labels = get(config, :Labels, nothing)
    name isa AbstractString && startswith(name, "/") && (name = name[2:end])
    matches_resource(receipt, scope, id, name, labels) ||
        throw(ArgumentError("container identity does not match its ownership receipt"))
    return String(id)
end

function container_absent(runner, host, receipt; invoke)
    # An inspect error is not absence: a working engine must return an empty
    # bounded inventory for the exact recorded ID, or exact generated name.
    filter = receipt.physical_id === nothing ?
        "name=" * resource_name(receipt) : "id=" * receipt.physical_id
    result = scoped_container_command(runner, host,
        ["container", "ls", "--all", "--no-trunc", "--filter", filter,
            "--format", "{{.ID}} {{.Names}}"]; invoke)
    result.exitcode == 0 || return false
    ncodeunits(result.output) <= 1024^2 || return false
    lines = split(strip(result.output), '\n'; keepempty=false)
    length(lines) <= 256 || return false
    for line in lines
        parts = split(strip(line); limit=2)
        length(parts) == 2 && occursin(r"^[a-f0-9]{64}$", parts[1]) || return false
        names = split(parts[2], ',')
        all(name -> occursin(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$", name), names) || return false
        if receipt.physical_id === nothing
            resource_name(receipt) in names && return false
        else
            parts[1] == receipt.physical_id && return false
        end
    end
    return true
end

"""
    remove_owned_container!(journal, runner, host, receipt) -> Bool

Retire only a receipt-matched container in the freshly inspected engine scope.
Unbound acquisition intents resolve through their generated name and every owned
label before being bound to the full ID. Stop and removal always target that ID,
never a name, wildcard, label filter or a general prune command.

A failed command is not proof of absence. Forget the receipt only after a working
engine reports the exact resource absent in the same scope. Return false on an
unresolved removal; identity mismatches throw and retain the receipt. This cleanup
path does not require CPU/memory launch prerequisites to remain available.
"""
function remove_owned_container!(journal::ResourceJournal, runner::CommandRunner,
        host::ContainerHostCheck, receipt::ResourceReceipt;
        probe=arguments -> container_probe(runner, arguments),
        invoke=arguments -> run_owned_command!(runner,
            setenv(Cmd(arguments), container_command_environment())),
        machine_id=strip(read("/etc/machine-id", String)),
        storage_identity=local_container_storage_identity)
    receipt.backend == host.engine.name || throw(ArgumentError("container receipt belongs to another engine"))
    current = lock(journal.lock) do
        receipts = resource_receipts(journal)
        receipt.journal_id == journal.id || throw(ArgumentError("resource belongs to another journal"))
        any(value -> value.id == receipt.id, receipts) ? current_receipt(journal, receipt) : nothing
    end
    scope = container_scope(runner, host; probe, machine_id, storage_identity)
    if current === nothing
        scope == receipt.scope || throw(ArgumentError("container engine scope changed; original resource is unresolved"))
        absent = container_absent(runner, host, receipt; invoke)
        return absent && container_scope(runner, host; probe, machine_id, storage_identity) == receipt.scope
    end
    scope == current.scope || throw(ArgumentError("container engine scope changed; original resource is unresolved"))
    id = inspect_owned_container(runner, host, current, scope; invoke)
    if id === nothing
        container_absent(runner, host, current; invoke) || return false
    else
        current = bind_resource!(journal, current, id)
        # Scientific/terminal owners normally request graceful protocol shutdown
        # first. Physical teardown still provides a finite stop then hard removal.
        for arguments in (["container", "stop", "--time", "2", id], ["container", "rm", "--force", id])
            container_scope(runner, host; probe, machine_id, storage_identity) == current.scope ||
                throw(ArgumentError("container engine changed during cleanup"))
            scoped_container_command(runner, host, arguments; invoke)
        end
        container_absent(runner, host, current; invoke) || return false
    end
    # Do not let a replaced engine's empty inventory release the original receipt.
    container_scope(runner, host; probe, machine_id, storage_identity) == current.scope ||
        throw(ArgumentError("container engine changed during cleanup"))
    forget_resource!(journal, current)
    return true
end

"""
    recover_containers!(journal, runner, host)

Attempt every recorded container belonging to this engine, joining failures before
reporting unresolved recovery. Other backends' receipts remain untouched for their
own adapters. The root physical driver must reconcile all backends before its
agent can announce availability. No preparation or live lease is restored.
"""
function recover_containers!(journal::ResourceJournal, runner::CommandRunner, host::ContainerHostCheck;
        options...)
    receipts = filter(receipt -> receipt.backend == host.engine.name, resource_receipts(journal))
    complete = true
    for receipt in receipts
        success = try
            remove_owned_container!(journal, runner, host, receipt; options...)
        catch
            false
        end
        complete &= success
    end
    complete || throw(ArgumentError("owned container recovery remains unresolved"))
    return nothing
end
