"""
    PublishedAssets

Share bounded, committed publication snapshots between the Bonito publisher and
the runtime gateway. No network services or application resources are owned here.
"""
module PublishedAssets
using UUIDs

const COMMIT_FILE = ".lcm-publication"
const MAX_FILE_BYTES = 64 * 1024^2
const MAX_SITE_BYTES = 128 * 1024^2

function revision(directory)
    file = joinpath(directory, COMMIT_FILE)
    isfile(file) || return ""
    islink(file) && throw(ArgumentError("publication marker must not be a symbolic link"))
    value = strip(open(io -> String(read(io, 64)), file))
    tryparse(UUID, value) === nothing && throw(ArgumentError("invalid publication marker"))
    return value
end

function snapshot(directory)
    isfile(joinpath(directory, "index.html")) || throw(ArgumentError("published site has no index.html"))
    files = Dict{String,String}()
    assets = Dict{String,Vector{UInt8}}()
    total = 0
    for (parent, dirs, names) in walkdir(directory; follow_symlinks=false)
        any(name -> islink(joinpath(parent, name)), [dirs; names]) &&
            throw(ArgumentError("published site must not contain symbolic links"))
        for name in names
            name == COMMIT_FILE && continue
            file = joinpath(parent, name)
            route = "/" * replace(relpath(file, directory), '\\'=>'/')
            any(prefix -> route == prefix || startswith(route, prefix * "/"), ("/runtime", "/applications")) &&
                throw(ArgumentError("publication contains a reserved runtime route"))
            realpath(file) == file || throw(ArgumentError("publication escaped its directory"))
            bytes = open(io -> read(io, min(MAX_FILE_BYTES, MAX_SITE_BYTES - total) + 1), file)
            length(bytes) <= MAX_FILE_BYTES || throw(ArgumentError("published file exceeds 64 MiB"))
            total += length(bytes)
            total <= MAX_SITE_BYTES || throw(ArgumentError("publication exceeds 128 MiB"))
            assets[route] = bytes
            files[route] = file
            if name == "index.html"
                prefix = dirname(route)
                for alias in (prefix, prefix == "/" ? "/" : prefix * "/")
                    files[alias] = file
                    assets[alias] = bytes
                end
            end
        end
    end
    return files, assets
end

"""
    PublishedSite(directory)

Capture a validated publication of at most 128 MiB, with a 64 MiB per-file limit.
Requests use immutable byte snapshots, never partially rebuilt files. A new
`publish!` marker replaces the complete snapshot on the next request. Keep the
previous revision's non-HTML assets for in-flight pages; removed page routes do
not survive. An invalid replacement retains the last valid publication.
"""
mutable struct PublishedSite
    "Validated, non-symlink publication directory."
    directory::String
    "Current exact public routes and their original file paths."
    files::Dict{String,String}
    "Current immutable response bodies, including index aliases."
    assets::Dict{String,Vector{UInt8}}
    "Previous revision's assets, without document routes."
    previous::Dict{String,Vector{UInt8}}
    "Last successfully loaded commit marker."
    revision::String
    "Rejected marker, avoiding repeated validation of the same failed build."
    rejected::String
    "Serialize publication replacement and request lookup."
    lock::ReentrantLock
end

function PublishedSite(directory::AbstractString)
    root = realpath(directory)
    before = revision(root)
    files, assets = snapshot(root)
    revision(root) == before || throw(ArgumentError("publication changed while indexing"))
    return PublishedSite(root, files, assets, Dict{String,Vector{UInt8}}(), before, "", ReentrantLock())
end

"""Return committed bytes for an exact route, or `nothing`; never resolve request paths on disk."""
function published_asset(site::PublishedSite, path::String)
    lock(site.lock) do
        current = try revision(site.directory) catch; return get(site.assets, path, nothing) end
        if !isempty(current) && current != site.revision && current != site.rejected
            try
                files, assets = snapshot(site.directory)
                revision(site.directory) == current || return get(site.assets, path, nothing)
                previous = Dict(route=>bytes for (route, bytes) in site.assets
                    if haskey(site.files, route) && lowercase(splitext(site.files[route])[2]) != ".html")
                site.files, site.assets, site.previous = files, assets, previous
                site.revision = current
            catch
                site.rejected = current
                @warn "Rejected incomplete or invalid publication; retaining the last valid snapshot"
            end
        end
        return get(site.assets, path) do
            get(site.previous, path, nothing)
        end
    end
end

"""
    publish!(directory)

Validate a completed build and atomically replace its commit marker. Existing
publishers adopt it without restarting applications. Call only after a successful
full-site or presentation render; failed renders must not publish a marker.
Return the new revision string. Validation errors leave the previous marker intact.
"""
function publish!(directory::AbstractString)
    root = realpath(directory)
    snapshot(root)
    value = string(uuid4())
    path, io = mktemp(root)
    try
        write(io, value)
        close(io)
        # rename replaces one file atomically; no delete-then-move window.
        Base.Filesystem.rename(path, joinpath(root, COMMIT_FILE))
    finally
        isopen(io) && close(io)
        isfile(path) && rm(path)
    end
    return value
end

export PublishedSite, published_asset, publish!
end
