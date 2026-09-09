"""
    EnvironmentFingerprint

Describe a bounded, read-only digest of a native profile's declared source
closure. This identifies inputs to a trusted Julia process, not installed
artifact integrity, successful preparation or operating-system isolation.
"""
struct EnvironmentFingerprint
    "Lowercase SHA-256 digest, including Julia version and target platform."
    digest::String
    "Root Julia package name; never evaluated by the fingerprint reader."
    package::String
    "Root Julia package UUID."
    uuid::UUID
    "Number of source/configuration files read."
    files::Int
    "Total source/configuration bytes read."
    bytes::Int
end

"""Maintain finite source-file accounting during one read-only inspection."""
mutable struct SourceInventory
    "Stable logical file names and content digests; no source contents retained."
    hashes::Dict{String,String}
    "Source identity/metadata retained until the complete inspection is rechecked."
    stamps::Dict{String,Tuple}
    "Total file and directory visits, including empty directories."
    visits::Int
    "Number of bytes already read."
    bytes::Int
    "Maximum number of files and directory entries per directory."
    max_files::Int
    "Maximum bytes in one file."
    max_file_bytes::Int
    "Maximum aggregate bytes."
    max_bytes::Int
end

function source_path(path::AbstractString)
    absolute = abspath(path)
    current = absolute
    while true
        islink(current) && throw(ArgumentError("native source paths cannot contain symbolic links"))
        parent = dirname(current)
        parent == current && break
        current = parent
    end
    return absolute
end

source_stamp(info) = (info.device, info.inode, info.size, info.mtime, info.ctime, info.mode)

function remember_source!(inventory::SourceInventory, path::String)
    inventory.visits += 1
    inventory.visits <= inventory.max_files || throw(ArgumentError("native source entry limit exceeded"))
    source_path(path)
    stamp = source_stamp(lstat(path))
    previous = get(inventory.stamps, path, stamp)
    previous == stamp || throw(ArgumentError("native source changed during fingerprint inspection"))
    inventory.stamps[path] = stamp
    return stamp
end

function verify_source_stamps(inventory::SourceInventory)
    for (path, stamp) in inventory.stamps
        source_path(path)
        source_stamp(lstat(path)) == stamp ||
            throw(ArgumentError("native source changed during fingerprint inspection"))
    end
    return nothing
end

function fingerprint_file!(inventory::SourceInventory, path::String, label::String)
    haskey(inventory.hashes, label) && throw(ArgumentError("duplicate native source identity"))
    length(inventory.hashes) < inventory.max_files || throw(ArgumentError("native source file limit exceeded"))
    source_path(path)
    remember_source!(inventory, path)
    before = lstat(path)
    isfile(before) || throw(ArgumentError("native source must be a regular file"))
    before.size <= inventory.max_file_bytes || throw(ArgumentError("native source file exceeds its byte limit"))
    inventory.bytes + before.size <= inventory.max_bytes || throw(ArgumentError("native source byte limit exceeded"))
    bytes = open(path, "r") do io
        # Cap the read even if an operator changes a file after the size check.
        read(io, min(inventory.max_file_bytes, inventory.max_bytes - inventory.bytes) + 1)
    end
    source_stamp(before) == source_stamp(lstat(path)) && length(bytes) == before.size ||
        throw(ArgumentError("native source changed during fingerprint inspection"))
    inventory.bytes += length(bytes)
    inventory.hashes[label] = bytes2hex(sha256(bytes))
    return bytes
end

function fingerprint_toml!(inventory::SourceInventory, path::String, label::String)
    bytes = fingerprint_file!(inventory, path, label)
    return try
        TOML.parse(String(bytes))
    catch
        # Parser errors may include private file contents or local paths.
        throw(ArgumentError("invalid native source TOML"))
    end
end

function fingerprint_tree!(inventory::SourceInventory, directory::String, label::String, depth::Int=0)
    depth < 48 || throw(ArgumentError("native source directory depth exceeded"))
    source_path(directory)
    isdir(directory) || throw(ArgumentError("native source directory is missing"))
    remember_source!(inventory, directory)
    before = source_stamp(lstat(directory))
    names = readdir(directory)
    length(names) <= inventory.max_files || throw(ArgumentError("native source directory entry limit exceeded"))
    for name in names
        path = joinpath(directory, name)
        source_path(path)
        child = label * "/" * name
        if isdir(path)
            fingerprint_tree!(inventory, path, child, depth + 1)
        else
            fingerprint_file!(inventory, path, child)
        end
    end
    before == source_stamp(lstat(directory)) ||
        throw(ArgumentError("native source directory changed during fingerprint inspection"))
    return nothing
end

function source_package(data)
    name = get(data, "name", nothing)
    name isa String && occursin(r"^[A-Za-z][A-Za-z0-9_]{0,127}$", name) ||
        throw(ArgumentError("native profile requires a named Julia package"))
    uuid = try
        UUID(get(data, "uuid", ""))
    catch
        throw(ArgumentError("native source package UUID is invalid"))
    end
    return name, uuid
end

function fingerprint_package!(inventory::SourceInventory, directory::String, label::String)
    source_path(directory)
    remember_source!(inventory, directory)
    (ispath(joinpath(directory, "JuliaProject.toml")) || islink(joinpath(directory, "JuliaProject.toml"))) &&
        throw(ArgumentError("native profile source inspection requires Project.toml, not JuliaProject.toml"))
    project = fingerprint_toml!(inventory, joinpath(directory, "Project.toml"), label * "/Project.toml")
    haskey(project, "manifest") && throw(ArgumentError("custom native manifest locations are unsupported"))
    identity = source_package(project)
    entrypoint = joinpath(directory, "src", first(identity) * ".jl")
    isfile(entrypoint) || throw(ArgumentError("native source package entry point is missing"))
    fingerprint_tree!(inventory, joinpath(directory, "src"), label * "/src")
    for name in ("ext", "deps")
        path = joinpath(directory, name)
        (ispath(path) || islink(path)) && fingerprint_tree!(inventory, path, label * "/" * name)
    end
    for name in ("Artifacts.toml", "JuliaArtifacts.toml", "LocalPreferences.toml", "JuliaLocalPreferences.toml")
        path = joinpath(directory, name)
        (ispath(path) || islink(path)) && fingerprint_file!(inventory, path, label * "/" * name)
    end
    extras = joinpath(directory, "RuntimeSources.toml")
    if ispath(extras) || islink(extras)
        data = fingerprint_toml!(inventory, extras, label * "/RuntimeSources.toml")
        Set(keys(data)) == Set(("schema_version", "files")) && get(data, "schema_version", nothing) === 1 ||
            throw(ArgumentError("unsupported native extra-source declaration"))
        files = data["files"]
        files isa Vector && length(files) <= 64 && all(p -> p isa String, files) && allunique(files) ||
            throw(ArgumentError("native extra sources must be a bounded, distinct file list"))
        targets = Set{String}()
        for file in sort(files)
            !isempty(file) && ncodeunits(file) <= 512 && !isabspath(file) &&
                !occursin(r"[\x00\r\n*?\[\]\\]", file) ||
                throw(ArgumentError("native extra source must be an explicit relative file"))
            target = source_path(joinpath(directory, file))
            target in targets && throw(ArgumentError("duplicate native extra source target"))
            push!(targets, target)
            fingerprint_file!(inventory, target, label * "/extra/" * file)
        end
    end
    return project
end

function native_manifest(directory::String)
    version = "$(VERSION.major).$(VERSION.minor)"
    for name in ("JuliaManifest-v$version.toml", "JuliaManifest.toml")
        (ispath(joinpath(directory, name)) || islink(joinpath(directory, name))) &&
            throw(ArgumentError("native profile source inspection requires Manifest.toml naming"))
    end
    versioned = "Manifest-v$version.toml"
    return ispath(joinpath(directory, versioned)) || islink(joinpath(directory, versioned)) ? versioned : "Manifest.toml"
end

function native_stdlib_identities()
    identities = Dict{UUID,String}()
    for name in readdir(Sys.STDLIB)
        path = joinpath(Sys.STDLIB, name, "Project.toml")
        isfile(path) || continue
        package = TOML.parsefile(path)
        identities[UUID(package["uuid"])] = package["name"]
    end
    return identities
end

"""
    native_environment_fingerprint(project; max_files=16384,
        max_file_bytes=33554432, max_bytes=268435456) -> EnvironmentFingerprint

Hash the selected manifest, root package, and every manifest path dependency's
Project.toml, src/, ext/, deps/, artifact declarations and local preferences.
RuntimeSources.toml may list up to 64 additional relative files for deliberate
cross-package includes. Package dependencies outside the local source tree are
identified by the manifest's pinned Git tree hashes, not imported or downloaded.

The digest includes this Julia version and target platform. Files are sorted by
logical package identity, not absolute checkout location. Manifest/Project bytes
remain literal: relocating an environment with absolute source paths changes its
digest. Alternate JuliaProject/JuliaManifest names and symbolic links fail closed.
Version-specific Manifest-vMAJOR.MINOR.toml takes precedence over Manifest.toml.

This is an operator-owned source check, not a sandbox or an attestation of a
mutable package depot. Inspection never connects, instantiates, prepares or runs
an environment. Launch must recheck the digest; changing trusted sources while
an executor is running requires retiring that executor, not retaining readiness.

# Errors

Throw ArgumentError for malformed, unsupported, changed or excessive sources.
Filesystem errors remain filesystem errors; callers must not expose their paths
to browser diagnostics. A source inspection alone cannot advertise readiness.
"""
function native_environment_fingerprint(project::AbstractString;
        max_files=16384, max_file_bytes=32 * 1024^2, max_bytes=256 * 1024^2)
    all(v -> v isa Integer && !(v isa Bool) && v > 0, (max_files, max_file_bytes, max_bytes)) &&
        max_files <= 16384 && max_file_bytes <= 32 * 1024^2 && max_bytes <= 256 * 1024^2 ||
        throw(ArgumentError("native fingerprint limits exceed supported bounds"))
    directory = source_path(project)
    inventory = SourceInventory(Dict{String,String}(), Dict{String,Tuple}(), 0, 0,
        max_files, max_file_bytes, max_bytes)
    root = fingerprint_package!(inventory, directory, "root")
    name, uuid = source_package(root)
    manifest_name = native_manifest(directory)
    manifest = fingerprint_toml!(inventory, joinpath(directory, manifest_name), "root/" * manifest_name)
    get(manifest, "manifest_format", nothing) == "2.0" || throw(ArgumentError("native profile requires manifest format 2.0"))
    get(manifest, "julia_version", nothing) == string(VERSION) ||
        throw(ArgumentError("native profile manifest must match this exact Julia version"))
    records = get(manifest, "deps", nothing)
    records isa Dict && length(records) <= 1024 || throw(ArgumentError("native manifest dependencies are invalid or excessive"))
    identities = Dict{UUID,String}()
    packages = [(name, uuid, root)]
    stdlibs = native_stdlib_identities()
    for dependency in sort!(collect(keys(records)))
        entries = records[dependency]
        entries isa Vector && 1 <= length(entries) <= 8 || throw(ArgumentError("invalid native manifest package entries"))
        for entry in entries
            entry isa Dict || throw(ArgumentError("invalid native manifest package entry"))
            id = try
                UUID(get(entry, "uuid", ""))
            catch
                throw(ArgumentError("invalid native manifest package UUID"))
            end
            haskey(identities, id) && throw(ArgumentError("duplicate native manifest package UUID"))
            identities[id] = dependency
            if haskey(entry, "path")
                path = entry["path"]
                path isa String && !isempty(path) && ncodeunits(path) <= 4096 &&
                    !occursin(r"[\x00\r\n]", path) || throw(ArgumentError("invalid native manifest source path"))
                package = fingerprint_package!(inventory, source_path(joinpath(directory, path)), "package/" * string(id))
                source_package(package) == (dependency, id) || throw(ArgumentError("native manifest source identity mismatch"))
                !haskey(entry, "version") || get(package, "version", nothing) == entry["version"] ||
                    throw(ArgumentError("native manifest source version mismatch"))
                push!(packages, (dependency, id, package))
            else
                tree = get(entry, "git-tree-sha1", nothing)
                (tree isa String && occursin(r"^[a-f0-9]{40}$", tree)) ||
                    (tree === nothing && get(stdlibs, id, nothing) == dependency) ||
                    throw(ArgumentError("native dependency lacks a pinned source tree"))
            end
        end
    end
    for (_, _, package) in packages
        dependencies = get(package, "deps", Dict{String,Any}())
        dependencies isa Dict || throw(ArgumentError("invalid native package dependencies"))
        for (dependency, id) in dependencies
            resolved = try
                UUID(id)
            catch
                throw(ArgumentError("invalid native package dependency UUID"))
            end
            get(identities, resolved, nothing) == dependency ||
                throw(ArgumentError("native package dependencies do not match the selected manifest"))
        end
    end
    verify_source_stamps(inventory)
    encoded = IOBuffer()
    println(encoded, "lcm-native-source-v1\n", VERSION, '\n', Sys.MACHINE, '\n', Sys.WORD_SIZE)
    for label in sort!(collect(keys(inventory.hashes)))
        # A length prefix avoids ambiguous labels when filenames contain whitespace.
        println(encoded, ncodeunits(label), ':', label, ':', inventory.hashes[label])
    end
    return EnvironmentFingerprint(bytes2hex(sha256(take!(encoded))), name, uuid,
        length(inventory.hashes), inventory.bytes)
end

"""
    verify_native_environment(profile) -> EnvironmentFingerprint

Recompute and compare an approved native scientific profile's source identity.
Reject other isolation kinds and mismatched inputs. This check must precede a
native launch, but does not replace resource preflight or executor preparation.
"""
function verify_native_environment(profile::ProfileDefinition)
    profile.kind == :scientific && profile.isolation == :trusted_process ||
        throw(ArgumentError("native verification requires a trusted scientific profile"))
    found = native_environment_fingerprint(profile.environment)
    found.digest == profile.fingerprint || throw(ArgumentError("native environment fingerprint does not match approval"))
    return found
end
