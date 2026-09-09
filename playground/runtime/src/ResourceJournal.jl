"""
    ResourceReceipt

Persist an acquisition intent before launching a physical resource, then bind its
immutable runtime ID after creation. A receipt identifies cleanup targets only:
it is neither a live lease nor preparation evidence.
"""
struct ResourceReceipt
    "Fresh acquisition UUID, also used in the generated resource name."
    id::UUID
    "Persistent identity of the private journal."
    journal_id::UUID
    "Exact assignment which requested the resource."
    fence::AssignmentFence
    "Podman, Docker, or a native systemd unit."
    backend::Symbol
    "SHA-256 identity of the inspected local engine or native supervisor scope."
    scope::String
    "Full container ID or systemd invocation ID; nothing before binding."
    physical_id::Union{Nothing,String}
end
Base.show(io::IO, receipt::ResourceReceipt) = print(io, "ResourceReceipt(", receipt.id, ", ", receipt.backend, ", <private>)")
Base.show(io::IO, ::MIME"text/plain", receipt::ResourceReceipt) = show(io, receipt)

"""
    ResourceJournal(root, worker_id; capacity=256)

Open a private Linux receipt directory under an exclusive kernel-held lock. The
directory must belong to this effective UID, have mode 0700, and contain no linked
path components. Receipt files are private, bounded JSON; no process is launched.
A fresh persistent UUID distinguishes this journal from another agent directory.

Recovery reads stale intents and bound identities; it never reconstructs process
authority from a saved PID. Unknown, malformed or foreign files fail closed and
are not deleted. Closing the journal releases its lock, not its recorded resources.
Physical drivers must finish their own teardown before closing it.
"""
mutable struct ResourceJournal
    "Dedicated operator-owned directory, never a user upload destination."
    root::String
    "Worker subject identity for every receipt."
    worker_id::String
    "Persistent UUID matching the private owner marker."
    id::UUID
    "Maximum recorded resources, including unresolved acquisition/cleanup."
    capacity::Int
    "Kernel-held single-writer lock."
    ownership::Base.Filesystem.File
    "Opened directory descriptor used to verify identity and sync mutations."
    directory::Base.Filesystem.File
    "Serialize this journal's bounded local operations."
    lock::ReentrantLock
    "Whether this handle can no longer authorize journal mutations."
    closed::Bool
end
Base.show(io::IO, journal::ResourceJournal) = print(io, "ResourceJournal(", journal.worker_id, ", <private>)")
Base.show(io::IO, ::MIME"text/plain", journal::ResourceJournal) = show(io, journal)

const RESOURCE_RECEIPT_BYTES = 8192
const RESOURCE_OWNER_FILE = "owner.json"
const RESOURCE_LOCK_FILE = "owner.lock"

function receipt_private_stat(metadata; directory=false)
    metadata.uid == ccall(:geteuid, Cuint, ()) ||
        throw(ArgumentError("resource journal ownership does not match this user"))
    if directory
        isdir(metadata) && metadata.mode & 0o7777 == 0o700 ||
            throw(ArgumentError("resource journal directory must have mode 0700"))
    else
        isfile(metadata) && metadata.mode & 0o7777 == 0o600 && metadata.nlink == 1 ||
            throw(ArgumentError("resource journal files must be private, regular and unlinked"))
    end
    return nothing
end

function receipt_open(path::String; create=false, exclusive=false, writable=false, directory=false)
    flags = Base.Filesystem.JL_O_CLOEXEC | Base.Filesystem.JL_O_NOFOLLOW
    directory || (flags |= Base.Filesystem.JL_O_NONBLOCK)
    flags |= writable ? Base.Filesystem.JL_O_RDWR : Base.Filesystem.JL_O_RDONLY
    create && (flags |= Base.Filesystem.JL_O_CREAT)
    exclusive && (flags |= Base.Filesystem.JL_O_EXCL)
    directory && (flags |= Base.Filesystem.JL_O_DIRECTORY)
    return Base.Filesystem.open(path, flags, 0o600)
end

function receipt_sync(file::Base.Filesystem.File)
    # A slow storage flush must not monopolize the agent's Julia scheduler.
    result = Base.@threadcall(:fsync, Cint, (Cint,), fd(file))
    result == 0 || throw(ArgumentError("resource journal synchronization failed"))
    return nothing
end

function receipt_read(path::String)
    receipt_private_stat(lstat(path))
    file = receipt_open(path)
    try
        before = stat(file)
        receipt_private_stat(before)
        0 < before.size <= RESOURCE_RECEIPT_BYTES ||
            throw(ArgumentError("resource receipt exceeds its size contract"))
        payload = read(file, RESOURCE_RECEIPT_BYTES + 1)
        length(payload) == before.size && eof(file) ||
            throw(ArgumentError("resource receipt changed during inspection"))
        after = stat(file)
        (before.size, before.mtime, before.ctime) == (after.size, after.mtime, after.ctime) ||
            throw(ArgumentError("resource receipt changed during inspection"))
        parsed = try JSON3.read(payload) catch; nothing end
        parsed isa JSON3.Object || throw(ArgumentError("resource receipt must be valid JSON"))
        return parsed
    finally
        close(file)
    end
end

function receipt_shape(object, expected)
    length(object) == length(expected) && Set(keys(object)) == Set(expected) ||
        throw(ArgumentError("resource receipt has missing, duplicate or unsupported fields"))
    return nothing
end

function journal_owner(object, worker_id)
    receipt_shape(object, (:schema_version, :kind, :worker_id, :journal_id))
    object.schema_version === 1 && object.kind == "lcm-resource-journal" &&
        object.worker_id == worker_id || throw(ArgumentError("resource journal identity does not match"))
    object.journal_id isa String || throw(ArgumentError("resource journal UUID must be text"))
    return UUID(Protocol.runtime_uuid(object.journal_id))
end

function receipt_payload(receipt::ResourceReceipt)
    return Dict("schema_version"=>1, "kind"=>"lcm-resource-receipt",
        "resource_id"=>string(receipt.id), "journal_id"=>string(receipt.journal_id),
        "fence"=>receipt.fence, "backend"=>string(receipt.backend), "scope"=>receipt.scope,
        "physical_id"=>receipt.physical_id)
end

function decode_resource_receipt(journal::ResourceJournal, object)
    receipt_shape(object, (:schema_version, :kind, :resource_id, :journal_id, :fence, :backend, :scope, :physical_id))
    object.schema_version === 1 && object.kind == "lcm-resource-receipt" ||
        throw(ArgumentError("unsupported resource receipt version"))
    object.resource_id isa String && object.journal_id isa String ||
        throw(ArgumentError("resource receipt UUIDs must be text"))
    id = UUID(Protocol.runtime_uuid(object.resource_id))
    UUID(Protocol.runtime_uuid(object.journal_id)) == journal.id ||
        throw(ArgumentError("resource receipt belongs to another journal"))
    object.backend in ("podman", "docker", "native") ||
        throw(ArgumentError("unsupported resource receipt backend"))
    fence = Protocol.decode_runtime_message(AssignmentFence, JSON3.write(object.fence))
    fence.worker_id == journal.worker_id || throw(ArgumentError("resource receipt belongs to another worker"))
    physical = object.physical_id
    backend = Symbol(object.backend)
    object.scope isa String && occursin(r"^[a-f0-9]{64}$", object.scope) ||
        throw(ArgumentError("resource receipt supervisor scope is invalid"))
    physical === nothing || (physical isa String && valid_physical_id(backend, physical)) ||
        throw(ArgumentError("resource receipt physical identity is invalid"))
    return ResourceReceipt(id, journal.id, fence, backend, object.scope, physical)
end

valid_physical_id(backend::Symbol, value::String) = occursin(
    backend == :native ? r"^[a-f0-9]{32}$" : r"^[a-f0-9]{64}$", value)

function validate_journal!(journal::ResourceJournal)
    journal.closed && throw(ArgumentError("resource journal is closed"))
    realpath(journal.root) == journal.root || throw(ArgumentError("resource journal path changed"))
    actual = lstat(journal.root)
    receipt_private_stat(actual; directory=true)
    original = stat(journal.directory)
    (actual.device, actual.inode) == (original.device, original.inode) ||
        throw(ArgumentError("resource journal directory was replaced"))
    current_lock = lstat(joinpath(journal.root, RESOURCE_LOCK_FILE))
    receipt_private_stat(current_lock)
    held_lock = stat(journal.ownership)
    (current_lock.device, current_lock.inode) == (held_lock.device, held_lock.inode) ||
        throw(ArgumentError("resource journal ownership lock was replaced"))
    journal_owner(receipt_read(joinpath(journal.root, RESOURCE_OWNER_FILE)), journal.worker_id) == journal.id ||
        throw(ArgumentError("resource journal marker changed"))
    return nothing
end

function receipt_atomic_write(root::String, directory::Base.Filesystem.File, name::String, object;
        replacing=false)
    payload = JSON3.write(object)
    ncodeunits(payload) <= RESOURCE_RECEIPT_BYTES || throw(ArgumentError("resource receipt exceeds its size contract"))
    path = joinpath(root, name)
    if ispath(path) || islink(path)
        replacing || throw(ArgumentError("resource receipt already exists"))
        receipt_private_stat(lstat(path))
    elseif replacing
        throw(ArgumentError("resource receipt disappeared before update"))
    end
    temporary = joinpath(root, ".pending-" * string(uuid4()))
    file = receipt_open(temporary; create=true, exclusive=true, writable=true)
    try
        write(file, payload) == ncodeunits(payload) || throw(ArgumentError("resource receipt write was incomplete"))
        receipt_sync(file)
        close(file)
        # Same directory, same filesystem, one kernel-locked writer. Do not use
        # mv(force=true), which can unlink the destination before renaming.
        Base.Filesystem.rename(temporary, path)
        receipt_sync(directory)
    finally
        isopen(file) && close(file)
        if isfile(temporary) && !islink(temporary)
            receipt_private_stat(lstat(temporary))
            rm(temporary)
            receipt_sync(directory)
        end
    end
    return nothing
end

function ResourceJournal(root::AbstractString, worker_id::AbstractString; capacity=256)
    Sys.islinux() || throw(ArgumentError("resource journals require Linux"))
    worker = Protocol.runtime_token(worker_id)
    capacity isa Integer && !(capacity isa Bool) && 1 <= capacity <= 256 ||
        throw(ArgumentError("resource journal capacity must be in 1:256"))
    path = abspath(root)
    path in ("/", homedir(), pwd(), dirname(pwd()), tempdir()) &&
        throw(ArgumentError("resource journal requires a dedicated directory"))
    ispath(path) || mkpath(path; mode=0o700)
    realpath(path) == path || throw(ArgumentError("resource journal cannot contain linked path components"))
    receipt_private_stat(lstat(path); directory=true)
    initial_names = readdir(path)
    if !(RESOURCE_OWNER_FILE in initial_names)
        all(==(RESOURCE_LOCK_FILE), initial_names) ||
            throw(ArgumentError("nonempty resource directory lacks its ownership marker"))
    end
    directory = receipt_open(path; directory=true)
    ownership = nothing
    journal = nothing
    try
        ownership = receipt_open(joinpath(path, RESOURCE_LOCK_FILE); create=true, writable=true)
        receipt_private_stat(stat(ownership))
        ccall(:flock, Cint, (Cint, Cint), fd(ownership), 6) == 0 ||
            throw(ArgumentError("another agent owns this resource journal"))
        names = readdir(path)
        owner_path = joinpath(path, RESOURCE_OWNER_FILE)
        if !ispath(owner_path) && !islink(owner_path)
            names == [RESOURCE_LOCK_FILE] ||
                throw(ArgumentError("nonempty resource directory lacks its ownership marker"))
            id = uuid4()
            receipt_atomic_write(path, directory, RESOURCE_OWNER_FILE,
                Dict("schema_version"=>1, "kind"=>"lcm-resource-journal",
                    "worker_id"=>worker, "journal_id"=>string(id)))
        else
            id = journal_owner(receipt_read(owner_path), worker)
        end
        journal = ResourceJournal(path, worker, id, capacity, ownership, directory, ReentrantLock(), false)
        # Validate every committed receipt before permitting any recovery action.
        resource_receipts(journal; recovering=true)
        return journal
    catch
        journal === nothing || (journal.closed = true)
        ownership === nothing || close(ownership)
        close(directory)
        rethrow()
    end
end

"""
    resource_receipts(journal) -> Vector{ResourceReceipt}

Read and validate every recorded acquisition. Pending atomic-write fragments are
removed only during exclusive startup recovery, after all committed records and
filenames validate. Unknown entries remain untouched and prevent startup.
"""
function resource_receipts(journal::ResourceJournal; recovering=false)
    lock(journal.lock) do
        validate_journal!(journal)
        names = readdir(journal.root)
        length(names) <= 2 * journal.capacity + 3 ||
            throw(ArgumentError("resource journal contains too many entries"))
        receipts = ResourceReceipt[]
        pending = String[]
        for name in names
            name in (RESOURCE_OWNER_FILE, RESOURCE_LOCK_FILE) && continue
            if occursin(r"^\.pending-[a-f0-9-]{36}$", name)
                recovering || throw(ArgumentError("resource journal contains an unresolved write"))
                Protocol.runtime_uuid(name[10:end])
                receipt_private_stat(lstat(joinpath(journal.root, name)))
                filesize(joinpath(journal.root, name)) <= RESOURCE_RECEIPT_BYTES ||
                    throw(ArgumentError("pending resource receipt exceeds its size contract"))
                push!(pending, name)
                continue
            end
            occursin(r"^[a-f0-9-]{36}\.json$", name) ||
                throw(ArgumentError("resource journal contains an unknown entry"))
            id = UUID(Protocol.runtime_uuid(name[1:end-5]))
            receipt = decode_resource_receipt(journal, receipt_read(joinpath(journal.root, name)))
            receipt.id == id || throw(ArgumentError("resource receipt filename does not match"))
            push!(receipts, receipt)
        end
        length(receipts) <= journal.capacity || throw(ArgumentError("resource journal capacity exceeded"))
        allunique(receipt.fence.lease_id for receipt in receipts) ||
            throw(ArgumentError("resource journal contains duplicate assignment ownership"))
        for name in pending
            rm(joinpath(journal.root, name))
        end
        isempty(pending) || receipt_sync(journal.directory)
        return receipts
    end
end

"""
    resource_name(receipt) -> String

Derive the only permitted physical name from its acquisition UUID. Native drivers
use a transient systemd service name; container drivers use the same fixed prefix.
"""
resource_name(receipt::ResourceReceipt) =
    "lcm-exec-" * string(receipt.id) * (receipt.backend == :native ? ".service" : "")

"""
    reserve_resource!(journal, fence, backend, scope) -> ResourceReceipt

Commit a new acquisition intent before the physical driver's launch command.
The exact fence may reuse an existing intent; a different fence with the same
lease ID is rejected. Backend is one of :podman, :docker or :native. Scope is
the driver's SHA-256 identity of its inspected local supervisor/engine, not a
browser value or an assumption derived only from the executable name.
"""
function reserve_resource!(journal::ResourceJournal, fence::AssignmentFence, backend::Symbol,
        scope::AbstractString)
    Protocol.validate(fence)
    fence.worker_id == journal.worker_id || throw(ArgumentError("resource fence belongs to another worker"))
    backend in (:podman, :docker, :native) || throw(ArgumentError("unsupported resource backend"))
    occursin(r"^[a-f0-9]{64}$", scope) || throw(ArgumentError("invalid resource supervisor scope"))
    return lock(journal.lock) do
        receipts = resource_receipts(journal)
        existing = findfirst(receipt -> receipt.fence.lease_id == fence.lease_id, receipts)
        if existing !== nothing
            receipt = receipts[existing]
            receipt.fence == fence && receipt.backend == backend && receipt.scope == scope ||
                throw(ArgumentError("resource lease already has a different acquisition"))
            return receipt
        end
        length(receipts) < journal.capacity || throw(ArgumentError("resource journal is full"))
        receipt = ResourceReceipt(uuid4(), journal.id, fence, backend, String(scope), nothing)
        receipt_atomic_write(journal.root, journal.directory, string(receipt.id) * ".json", receipt_payload(receipt))
        return receipt
    end
end

function current_receipt(journal::ResourceJournal, receipt::ResourceReceipt)
    validate_journal!(journal)
    receipt.journal_id == journal.id || throw(ArgumentError("resource belongs to another journal"))
    path = joinpath(journal.root, string(receipt.id) * ".json")
    current = decode_resource_receipt(journal, receipt_read(path))
    current.id == receipt.id && current.fence == receipt.fence &&
        current.backend == receipt.backend && current.scope == receipt.scope ||
        throw(ArgumentError("resource receipt identity changed"))
    return current
end

"""
    bind_resource!(journal, receipt, physical_id) -> ResourceReceipt

Atomically bind an intent to a full container ID or native unit invocation ID.
An existing binding may be confirmed but never replaced by a different process.
A failed write leaves the intent available for name-and-label recovery.
"""
function bind_resource!(journal::ResourceJournal, receipt::ResourceReceipt, physical_id::AbstractString)
    value = String(physical_id)
    valid_physical_id(receipt.backend, value) || throw(ArgumentError("invalid physical resource identity"))
    return lock(journal.lock) do
        current = current_receipt(journal, receipt)
        if current.physical_id !== nothing
            current.physical_id == value || throw(ArgumentError("resource is already bound to another physical identity"))
            return current
        end
        receipt.physical_id === nothing || throw(ArgumentError("resource binding was lost"))
        bound = ResourceReceipt(current.id, journal.id, current.fence, current.backend, current.scope, value)
        receipt_atomic_write(journal.root, journal.directory, string(bound.id) * ".json",
            receipt_payload(bound); replacing=true)
        return bound
    end
end

"""
    forget_resource!(journal, receipt)

Remove exactly the current receipt after the physical driver has independently
proved the resource absent. This method performs no signaling or container removal
and cannot release a lease. A stale pre-binding receipt cannot discard a binding.
"""
function forget_resource!(journal::ResourceJournal, receipt::ResourceReceipt)
    lock(journal.lock) do
        validate_journal!(journal)
        receipt.journal_id == journal.id || throw(ArgumentError("resource belongs to another journal"))
        path = joinpath(journal.root, string(receipt.id) * ".json")
        (ispath(path) || islink(path)) || return nothing
        current = current_receipt(journal, receipt)
        current.physical_id == receipt.physical_id || throw(ArgumentError("resource binding changed before receipt removal"))
        rm(path)
        receipt_sync(journal.directory)
    end
    return nothing
end

"""
    resource_labels(receipt) -> Dict{String,String}

Produce fixed ownership labels for physical launch and recovery checks. These
correlation identities are not credentials and do not establish live authority.
"""
function resource_labels(receipt::ResourceReceipt)
    return Dict(
        "org.linecablemodels.kind"=>"executor-v1",
        "org.linecablemodels.journal"=>string(receipt.journal_id),
        "org.linecablemodels.resource"=>string(receipt.id),
        "org.linecablemodels.scope"=>receipt.scope,
        "org.linecablemodels.worker"=>receipt.fence.worker_id,
        "org.linecablemodels.boot"=>receipt.fence.worker_boot,
        "org.linecablemodels.lease"=>receipt.fence.lease_id,
        "org.linecablemodels.fingerprint"=>receipt.fence.fingerprint)
end

"""
    matches_resource(receipt, scope, physical_id, name, labels) -> Bool

Require the observed supervisor/engine scope, generated name, every ownership
label, and any already-bound full physical identity to match. An intent with no
binding still requires exact labels.
A matching receipt identifies cleanup ownership, not current execution authority.
"""
function matches_resource(receipt::ResourceReceipt, scope, physical_id, name, labels)
    scope == receipt.scope || return false
    physical_id isa AbstractString && valid_physical_id(receipt.backend, String(physical_id)) || return false
    name == resource_name(receipt) || return false
    receipt.physical_id === nothing || receipt.physical_id == physical_id || return false
    labels isa AbstractDict || labels isa JSON3.Object || return false
    return all(get(labels, key, nothing) == value for (key, value) in resource_labels(receipt))
end

"""Release the journal lock; leave all acquisition receipts intact for recovery."""
function Base.close(journal::ResourceJournal)
    lock(journal.lock) do
        journal.closed && return nothing
        close(journal.ownership)
        close(journal.directory)
        journal.closed = true
    end
    return nothing
end
