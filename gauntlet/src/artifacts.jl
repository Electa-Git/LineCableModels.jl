function _tree_files(root::AbstractString)
    result = String[]
    for (directory, _, names) in walkdir(root), name in names

        push!(result, relpath(joinpath(directory, name), root))
    end
    return sort!(result)
end

"""Return a deterministic SHA-256 identity for file names and contents."""
function tree_hash(root::AbstractString)
    isdir(root) || throw(ArgumentError("artifact tree is absent: $root"))
    context = SHA.SHA2_256_CTX()
    for relative in _tree_files(root)
        SHA.update!(context, codeunits(replace(relative, '\\' => '/')))
        SHA.update!(context, UInt8[0])
        open(joinpath(root, relative), "r") do io
            buffer = Vector{UInt8}(undef, 1024 * 1024)
            while !eof(io)
                count = readbytes!(io, buffer)
                SHA.update!(context, @view buffer[1:count])
            end
        end
        SHA.update!(context, UInt8[0])
    end
    return bytes2hex(SHA.digest!(context))
end

function _archive(source::AbstractString, target::AbstractString)
    mkpath(dirname(target))
    open(target, "w") do output
        stream = CodecZlib.GzipCompressorStream(output; level = 9)
        try
            Tar.create(source, stream; portable = true)
        finally
            close(stream)
        end
    end
    return target
end

function _verify_archive(archive, expected_hash; dataset::Bool)
    return mktempdir() do destination
        open(archive, "r") do input
            stream = CodecZlib.GzipDecompressorStream(input)
            try
                Tar.extract(stream, destination)
            finally
                close(stream)
            end
        end
        observed = tree_hash(destination)
        observed == expected_hash || throw(ArgumentError(
            "archive extraction changed its tree hash: expected $expected_hash, got $observed",
        ))
        dataset && Dataset(destination)
        return observed
    end
end

function stage_artifacts(name::Symbol, args...; kwargs...)
    stage_artifacts(datasource(name), args...; kwargs...)
end

stage_artifacts(::FEM, args...; kwargs...) = _fem_not_implemented()

function stage_smoke(name::Symbol, args...; kwargs...)
    stage_smoke(datasource(name), args...; kwargs...)
end

stage_smoke(::FEM, args...; kwargs...) = _fem_not_implemented()
