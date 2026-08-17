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

function _verify_archive(archive, expected_hash; corpus::Bool)
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
        corpus && Corpus(destination)
        return observed
    end
end

"""Stage verified raw and normalized PSCAD corpora as local immutable archives."""
function stage_artifacts(
        staging::AbstractString;
        output::AbstractString = joinpath(_GAUNTLET_ROOT, ".artifacts"),
        repeat_normalized::Union{Nothing, AbstractString} = nothing
)
    raw = joinpath(abspath(staging), "raw")
    normalized = joinpath(abspath(staging), "normalized")
    isfile(joinpath(raw, "index.toml")) || throw(ArgumentError(
        "staging directory has no complete raw corpus",
    ))
    Corpus(normalized)
    raw_hash = tree_hash(raw)
    normalized_hash = tree_hash(normalized)
    index = TOML.parsefile(joinpath(normalized, "index.toml"))
    aliases = TOML.parsefile(joinpath(normalized, "aliases.toml"))["aliases"]
    if repeat_normalized !== nothing
        repeated_hash = tree_hash(repeat_normalized)
        repeated_hash == normalized_hash || throw(ArgumentError(
            "independent normalizations differ: $normalized_hash and $repeated_hash",
        ))
    end

    raw_archive = _archive(raw, joinpath(output, "pscad-raw-v1.tar.gz"))
    normalized_archive = _archive(
        normalized, joinpath(output, "pscad-normalized-v1.tar.gz")
    )
    _verify_archive(raw_archive, raw_hash; corpus = false)
    _verify_archive(normalized_archive, normalized_hash; corpus = true)
    result = Dict{String, Any}(
        "schema_version" => 1,
        "raw" => Dict(
            "archive" => basename(raw_archive),
            "archive_sha256" => _sha256(raw_archive),
            "tree_sha256" => raw_hash,
            "files" => length(_tree_files(raw)),
            "bytes" => sum(filesize(joinpath(raw, path)) for path in _tree_files(raw))
        ),
        "normalized" => Dict(
            "archive" => basename(normalized_archive),
            "archive_sha256" => _sha256(normalized_archive),
            "tree_sha256" => normalized_hash,
            "files" => length(_tree_files(normalized)),
            "bytes" => sum(
                filesize(joinpath(normalized, path)) for path in _tree_files(normalized)
            ),
            "successes" => length(index["cases"]),
            "rejections" => length(index["rejections"]),
            "aliases" => length(aliases)
        ),
        "reproduce" => "julia --project=gauntlet gauntlet/bin/stage_artifacts.jl $(abspath(staging))"
    )
    _case_toml(joinpath(output, "manifest.toml"), result)
    return result
end

"""Copy the fixed smoke selection from a verified normalized corpus."""
function stage_smoke(
        source::AbstractString,
        destination::AbstractString = _SMOKE_ROOT;
        configuration::AbstractString = joinpath(
            _GAUNTLET_ROOT, "config", "smoke.toml"
        )
)
    corpus = Corpus(source)
    selected = String.(TOML.parsefile(configuration)["cases"])
    all(id -> id in keys(corpus), selected) || throw(KeyError(
        "the smoke selection contains a case absent from the normalized corpus",
    ))
    source_index = TOML.parsefile(joinpath(corpus.root, "index.toml"))
    success = Dict{String, String}()
    rejections = Dict{String, String}()
    for id in selected
        if haskey(corpus.cases, id)
            relative = corpus.cases[id]
            target = joinpath(destination, relative)
            mkpath(dirname(target))
            cp(joinpath(corpus.root, relative), target; force = true)
            success[id] = relative
        else
            relative = corpus.rejections[id]
            target = joinpath(destination, relative)
            mkpath(dirname(target))
            cp(joinpath(corpus.root, relative), target; force = true)
            rejections[id] = relative
        end
    end
    selected_set = Set(selected)
    aliases = Dict(
        alias => target for (alias, target) in corpus.aliases if target in selected_set
    )
    family_counts = Dict{String, Int}()
    detailed = 0
    sparse = 0
    for id in keys(success)
        case = corpus[id]
        family = lowercase(String(nameof(typeof(case.family))))
        family_counts[family] = get(family_counts, family, 0) + 1
        length(frequencies(case.reference.phase)) == 1 ? (sparse += 1) : (detailed += 1)
    end
    index = Dict{String, Any}(
        "schema_version" => 1,
        "source_manifest_sha256" => source_index["source_manifest_sha256"],
        "cases" => success,
        "rejections" => rejections,
        "excluded" => source_index["excluded"],
        "family_counts" => family_counts,
        "detailed_cases" => detailed,
        "sparse_cases" => sparse,
        "recovered_fits" => count(
            id -> any(
                assumption -> assumption.subject === :artifact_recovery,
                corpus[id].assumptions
            ),
            keys(success))
    )
    _case_toml(joinpath(destination, "index.toml"), index)
    _case_toml(joinpath(destination, "aliases.toml"), Dict("aliases" => aliases))
    cp(joinpath(corpus.root, "schema.toml"), joinpath(destination, "schema.toml");
        force = true)
    staged = Corpus(destination)
    length(staged) == length(selected) || throw(DimensionMismatch(
        "staged smoke corpus does not contain all selected cases",
    ))
    return (;
        destination = abspath(destination),
        cases = length(success),
        rejections = length(rejections),
        aliases = length(aliases),
        tree_sha256 = tree_hash(destination)
    )
end
