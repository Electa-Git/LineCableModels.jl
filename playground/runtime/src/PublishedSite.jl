"""
    PublishedSite(directory)

Index an already built, read-only publication. Exact routes and index aliases
are captured once; symlinks and runtime-reserved routes cannot escape or bypass
the owned gateway. Constructing this index starts no services.
"""
struct PublishedSite
    "Exact public route to validated file mapping."
    files::Dict{String,String}
    function PublishedSite(directory::AbstractString)
        root = realpath(directory)
        isfile(joinpath(root, "index.html")) || throw(ArgumentError("published site has no index.html"))
        files = Dict{String,String}()
        for (parent, dirs, names) in walkdir(root; follow_symlinks=false)
            any(name -> islink(joinpath(parent, name)), [dirs; names]) &&
                throw(ArgumentError("published site must not contain symbolic links"))
            for name in names
                file = joinpath(parent, name)
                route = "/" * replace(relpath(file, root), '\\'=>'/')
                (route in ("/runtime", "/applications") ||
                    startswith(route, "/runtime/") || startswith(route, "/applications/")) &&
                    throw(ArgumentError("publication contains a reserved runtime route"))
                files[route] = file
                if name == "index.html"
                    prefix = dirname(route)
                    files[prefix] = file
                    files[prefix == "/" ? "/" : prefix * "/"] = file
                end
            end
        end
        return new(files)
    end
end

const STATIC_MIME = Dict(".html"=>"text/html; charset=utf-8", ".css"=>"text/css",
    ".js"=>"text/javascript", ".json"=>"application/json", ".svg"=>"image/svg+xml",
    ".png"=>"image/png", ".jpg"=>"image/jpeg", ".jpeg"=>"image/jpeg", ".gif"=>"image/gif",
    ".woff"=>"font/woff", ".woff2"=>"font/woff2", ".ttf"=>"font/ttf", ".pdf"=>"application/pdf")

function serve_published(stream, site::PublishedSite, path::String)
    file = get(site.files, path, nothing)
    isnothing(file) && return false
    # A file changed to a symlink after indexing must still fail closed.
    isfile(file) && !islink(file) && realpath(file) == file ||
        throw(AccessDenied(503, "Published asset is unavailable"))
    filesize(file) <= 64 * 1024^2 || throw(AccessDenied(413, "Published asset exceeds its size bound"))
    mime = get(STATIC_MIME, lowercase(splitext(file)[2]), "application/octet-stream")
    gateway_response(stream, 200, read(file); content_type=mime)
    return true
end
