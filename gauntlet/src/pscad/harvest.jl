"""Describe the detailed matrix evidence required from one PSCAD run."""
struct PSCADRequirement{T <: Real}
    suffixes::NTuple{4, String}
    rows::Int
    frequency_range::Tuple{T, T}
end

function required(::PSCAD)
    PSCADRequirement(
        ("_zm.out", "_zp.out", "_ym.out", "_yp.out"),
        101,
        (1.0e-3, 1.0e7)
    )
end

"""Record one identity-checked amendment to a PSCAD campaign."""
struct PSCADAmendment
    case_id::String
    specification_sha256::String
    manifest::String
    matrix_dimension::Int
end

"""Summarize completeness verification for a PSCAD amendment set."""
struct HarvestReport
    expected::Int
    verified::Int
    amendments::Vector{PSCADAmendment}
    errors::Vector{String}
end

function Base.show(io::IO, report::HarvestReport)
    print(io, "HarvestReport(", report.verified, '/', report.expected,
        " verified, errors=", length(report.errors), ')')
end

function _artifact_names(manifest)
    return lowercase.(basename.(_native_path.(String.(getindex.(manifest["artifacts"], "path")))))
end

"""Return canonical successful records missing declared PSCAD matrix evidence."""
function pending(::PSCAD, dataset_path::AbstractString; requirement = required(PSCAD()))
    dataset = _json_dict(abspath(dataset_path))
    records, _ = _canonical_records(dataset)
    root = dirname(abspath(dataset_path))
    result = NamedTuple[]
    for (id, record) in sort!(collect(records); by = first)
        String(record["status"]) == "success" || continue
        case_root = joinpath(root, _native_path(String(record["path"])))
        manifest_path = joinpath(case_root, "manifest.json")
        manifest = _json_dict(manifest_path)
        names = _artifact_names(manifest)
        missing = [suffix
                   for suffix in requirement.suffixes
                   if !any(endswith(suffix), names)]
        isempty(missing) || push!(result,
            (;
                case_id = id,
                specification_sha256 = String(record["specification_sha256"]),
                source_sha256 = String(record["source_sha256"]),
                manifest = manifest_path,
                missing
            ))
    end
    return result
end

pending(name::Symbol, args...; kwargs...) = pending(datasource(name), args...; kwargs...)

function _amendment_manifests(root::AbstractString)
    paths = String[]
    for (directory, _, names) in walkdir(abspath(root))
        "manifest.json" in names && push!(paths, joinpath(directory, "manifest.json"))
    end
    return sort!(paths)
end

function _matrix_tables(files, requirement::PSCADRequirement)
    tables = PSCADIO.Detailed[]
    dimensions = Int[]
    for suffix in requirement.suffixes
        path = _suffix_file(files, suffix)
        path === nothing && throw(ArgumentError("missing required PSCAD artifact $suffix"))
        table = PSCADIO.read_detailed(path)
        size(table.values, 1) == requirement.rows || throw(DimensionMismatch(
            "$(basename(path)) has $(size(table.values, 1)) rows; expected $(requirement.rows)",
        ))
        dimension = isqrt(size(table.values, 2))
        dimension^2 == size(table.values, 2) || throw(DimensionMismatch(
            "$(basename(path)) does not contain a square row-major matrix",
        ))
        push!(dimensions, dimension)
        push!(tables, table)
    end
    all(==(first(dimensions)), dimensions) || throw(DimensionMismatch(
        "PSCAD matrix dimensions differ within the required quartet",
    ))
    frequency = first(tables).frequency
    all(table -> table.frequency == frequency, tables) || throw(DimensionMismatch(
        "PSCAD matrix frequencies differ within the required quartet",
    ))
    first(frequency) ≈ first(requirement.frequency_range) || throw(DomainError(
        first(frequency), "PSCAD sweep starts outside the declared range"
    ))
    last(frequency) ≈ last(requirement.frequency_range) || throw(DomainError(
        last(frequency), "PSCAD sweep ends outside the declared range"
    ))
    return first(dimensions)
end

"""Verify identity, hashes, shapes, and frequency coverage of PSCAD amendments."""
function verify(
        ::PSCAD,
        dataset_path::AbstractString,
        amendment_root::AbstractString;
        requirement = required(PSCAD())
)
    expected_rows = pending(PSCAD(), dataset_path; requirement)
    expected = Dict(row.specification_sha256 => row for row in expected_rows)
    seen = Set{String}()
    amendments = PSCADAmendment[]
    errors = String[]
    for manifest_path in _amendment_manifests(amendment_root)
        try
            manifest = _json_dict(manifest_path)
            specification = String(manifest["specification_sha256"])
            specification in seen && throw(ArgumentError(
                "duplicate amendment specification $specification",
            ))
            push!(seen, specification)
            haskey(expected, specification) || throw(ArgumentError(
                "amendment specification is not missing from the canonical campaign",
            ))
            original = expected[specification]
            String(manifest["case_id"]) == original.case_id || throw(ArgumentError(
                "amendment case identifier differs from the canonical record",
            ))
            String(manifest["source_sha256"]) == original.source_sha256 ||
                throw(ArgumentError(
                    "amendment donor hash differs from the canonical record",
                ))
            String(manifest["status"]) == "success" || throw(ArgumentError(
                "amendment status is not success",
            ))
            files, _ = _artifact_paths(dirname(manifest_path), manifest)
            dimension = _matrix_tables(files, requirement)
            push!(amendments, PSCADAmendment(
                original.case_id, specification, manifest_path, dimension
            ))
        catch error
            error isa InterruptException && rethrow()
            push!(errors, "$(manifest_path): $(sprint(showerror, error))")
        end
    end
    missing = setdiff(keys(expected), seen)
    isempty(missing) || push!(errors,
        "missing $(length(missing)) expected amendment specifications")
    unexpected = setdiff(seen, keys(expected))
    isempty(unexpected) || push!(errors,
        "found $(length(unexpected)) unexpected amendment specifications")
    return HarvestReport(
        length(expected), length(amendments), sort!(amendments; by = x -> x.case_id), errors
    )
end

verify(name::Symbol, args...; kwargs...) = verify(datasource(name), args...; kwargs...)

function _amendment_map(dataset_path::AbstractString, amendment_root::AbstractString)
    report = verify(PSCAD(), dataset_path, amendment_root)
    isempty(report.errors) || throw(ErrorException(
        "PSCAD amendment verification failed:\n" * join(report.errors, '\n'),
    ))
    report.verified == report.expected || throw(DimensionMismatch(
        "verified $(report.verified) of $(report.expected) PSCAD amendments",
    ))
    return Dict(row.specification_sha256 => row.manifest for row in report.amendments)
end

struct LineKey
    iid::Int
    definition::String
    component_type::String
    name::String
end

struct Variant
    resistivity_ohm_m::Union{Nothing, Float64}
    formulations::Dict{String, Any}
    reduction::String
end

function _variant_dict(variant::Variant)
    Dict{String, Any}(
        "resistivity_ohm_m" => variant.resistivity_ohm_m,
        "formulations" => variant.formulations,
        "reduction" => variant.reduction
    )
end

function _canonical_json(value)
    value isa AbstractDict && return '{' *
           join(
               (JSON3.write(String(key)) * ':' * _canonical_json(value[key])
               for key in sort!(collect(keys(value)); by = string)), ','
           ) * '}'
    value isa Union{AbstractVector, Tuple} && return '[' *
           join((_canonical_json(item) for item in value), ',') * ']'
    return JSON3.write(value)
end

_sha256_json(value) = bytes2hex(SHA.sha256(_canonical_json(value)))

function _slug(value::AbstractString, limit::Int = 52)
    result = strip(replace(value, r"[^A-Za-z0-9._-]+" => "_"), ['.', '_', '-'])
    isempty(result) && (result = "unnamed")
    return first(result, min(length(result), limit))
end

function _write_pscad_json(path::AbstractString, value)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.pretty(io, value)
        write(io, '\n')
    end
    return path
end

function _project_identity(path::AbstractString)
    document = EzXML.readxml(path)
    root = EzXML.root(document)
    return (
        name = haskey(root, "name") ? root["name"] : splitext(basename(path))[1],
        version = haskey(root, "version") ? root["version"] : nothing
    )
end

function _normalized_lcp(path::AbstractString)
    bytes = read(path)
    text = try
        String(bytes)
    catch
        transcode(String, bytes)
    end
    normalized = replace(text, "\r\n" => "\n", '\r' => '\n')
    return join(rstrip.(split(normalized, '\n')), '\n') |> rstrip |> value -> value * "\n"
end

function _case_id(version::AbstractString, inputs)
    context = SHA.SHA2_256_CTX()
    SHA.update!(context, codeunits(version))
    SHA.update!(context, UInt8[0])
    for path in sort!(collect(inputs); by = path -> lowercase(basename(path)))
        SHA.update!(context, codeunits(lowercase(basename(path))))
        SHA.update!(context, UInt8[0])
        SHA.update!(context, codeunits(_normalized_lcp(path)))
        SHA.update!(context, UInt8[0])
    end
    return bytes2hex(SHA.digest!(context))
end

function _snapshot(root::AbstractString; recursive::Bool = true)
    isdir(root) || return Dict{String, Tuple{Float64, Int64}}()
    result = Dict{String, Tuple{Float64, Int64}}()
    iterator = recursive ? walkdir(root) : ((root, String[], readdir(root)),)
    for (directory, _, names) in iterator, name in names

        path = joinpath(directory, name)
        isfile(path) || continue
        result[relpath(path, root)] = (mtime(path), filesize(path))
    end
    return result
end

function _changed(before, after)
    return sort!([name for (name, stamp) in after if get(before, name, nothing) != stamp])
end

function _wait_for_artifacts(
        temp::AbstractString,
        before;
        timeout::Real,
        quiet::Real,
        process_root::AbstractString = pwd(),
        process_before = _snapshot(process_root; recursive = false),
        project_root::AbstractString,
        project_before
)
    deadline = time() + timeout
    last = (
        temp = _snapshot(temp),
        process = _snapshot(process_root; recursive = false),
        project = _snapshot(project_root; recursive = false)
    )
    changed_at = time()
    while time() < deadline
        sleep(0.25)
        current = (
            temp = _snapshot(temp),
            process = _snapshot(process_root; recursive = false),
            project = _snapshot(project_root; recursive = false)
        )
        if current != last
            last = current
            changed_at = time()
        end
        changed = vcat(
            _changed(before, current.temp),
            _changed(process_before, current.process),
            _changed(project_before, current.project)
        )
        input = any(name -> lowercase(splitext(name)[2]) in (".cli", ".tli"), changed)
        output = any(
            name -> lowercase(splitext(name)[2]) in (".clo", ".tlo", ".out"),
            changed
        )
        input && output && time() - changed_at >= quiet && return current
        for (root, baseline, snapshot) in (
            (temp, before, current.temp),
            (process_root, process_before, current.process),
            (project_root, project_before, current.project)
        )
            for name in _changed(baseline, snapshot)
                lowercase(splitext(name)[2]) == ".log" || continue
                text = lowercase(read(joinpath(root, name), String))
                any(marker -> occursin(marker, text),
                    ("tline.exe: error", "emttl ending!", "forrtl: severe")) &&
                    return current
            end
        end
    end
    return (
        temp = _snapshot(temp),
        process = _snapshot(process_root; recursive = false),
        project = _snapshot(project_root; recursive = false)
    )
end

function _solver_rejection(record::AbstractString, artifacts)
    markers = ("tline.exe: error", "emttl ending!", "forrtl: severe", "singularity occurs")
    for artifact in artifacts
        artifact["kind"] == "lcp_diagnostics" || continue
        path = joinpath(record, artifact["path"])
        text = lowercase(read(path, String))
        matched = filter(marker -> occursin(marker, text), markers)
        isempty(matched) || return Dict(
            "markers" => collect(matched),
            "log" => artifact["path"],
            "tail" =>
                join(last(split(text, '\n'), min(30, length(split(text, '\n')))), '\n')
        )
    end
    return nothing
end

function _classify_artifact(path::AbstractString, location::Symbol)
    extension = lowercase(splitext(path)[2])
    extension == ".out" && location in (:process_cwd, :project_cwd) &&
        return "detailed_frequency_output"
    return get(
        Dict(
            ".tli" => "normalized_lcp_input",
            ".cli" => "normalized_lcp_input",
            ".tlo" => "emtdc_constants",
            ".clo" => "emtdc_constants",
            ".out" => "human_readable_matrices",
            ".log" => "lcp_diagnostics"
        ),
        extension,
        "detailed_or_auxiliary_output")
end

function _copy_changed(source, names, destination, location::Symbol)
    rows = Dict{String, Any}[]
    subdirectory = location === :temp ? "" : String(location)
    for name in names
        input = joinpath(source, name)
        isfile(input) || continue
        relative = isempty(subdirectory) ? name : joinpath(subdirectory, name)
        output = joinpath(destination, relative)
        mkpath(dirname(output))
        cp(input, output; force = true)
        push!(rows,
            Dict(
                "path" => _native_path(joinpath("artifacts", relative)),
                "size" => filesize(output),
                "sha256" => _sha256(output),
                "kind" => _classify_artifact(output, location)
            ))
    end
    return rows
end

function _source_catalog(path::AbstractString)
    values = TOML.parsefile(path)
    return get(values, "sources", Any[])
end

function _extract_zip(path::AbstractString, destination::AbstractString)
    root = abspath(destination)
    mkpath(root)
    reader = ZipFile.Reader(path)
    try
        for file in reader.files
            target = normpath(joinpath(root, file.name))
            startswith(target, root * Base.Filesystem.path_separator) || target == root ||
                throw(ArgumentError("archive member escapes its destination: $(file.name)"))
            endswith(file.name, '/') && (mkpath(target); continue)
            mkpath(dirname(target))
            open(target, "w") do io
                write(io, read(file))
            end
        end
    finally
        close(reader)
    end
    return root
end

function harvest_sources(
        root::AbstractString,
        catalog::AbstractString;
        download::Bool = true,
        extract::Bool = true
)
    destination_root = joinpath(abspath(root), "donors", "_downloads")
    mkpath(destination_root)
    rows = Dict{String, Any}[]
    for source in _source_catalog(catalog)
        row = Dict{String, Any}(String(key) => value for (key, value) in source)
        destination = joinpath(destination_root, String(source["filename"]))
        if !isfile(destination)
            if !download
                row["status"] = "missing"
                push!(rows, row)
                continue
            end
            temporary = destination * ".part"
            Downloads.download(String(source["url"]), temporary)
            mv(temporary, destination; force = true)
        end
        observed = _sha256(destination)
        row["path"] = destination
        row["actual_sha256"] = observed
        if lowercase(observed) != lowercase(String(source["sha256"]))
            row["status"] = "hash_mismatch"
        else
            row["status"] = "verified"
            if extract && get(source, "kind", "") == "archive" &&
               haskey(source, "extract_to")
                row["extracted_to"] = _extract_zip(
                    destination, joinpath(abspath(root), String(source["extract_to"]))
                )
            end
        end
        push!(rows, row)
    end
    report = Dict(
        "created_at" => string(now(UTC)),
        "root" => abspath(root),
        "sources" => rows
    )
    _write_pscad_json(joinpath(abspath(root), "donors", "source_manifest.json"), report)
    return report
end

function catalog(root::AbstractString)
    projects = Dict{String, Any}[]
    for (directory, _, names) in walkdir(abspath(root)), name in names

        lowercase(splitext(name)[2]) == ".pscx" || continue
        path = joinpath(directory, name)
        push!(projects,
            Dict(
                "path" => path,
                "relative_path" => relpath(path, abspath(root)),
                "size" => filesize(path),
                "sha256" => _sha256(path)
            ))
    end
    sort!(projects; by = row -> lowercase(String(row["relative_path"])))
    return Dict("root" => abspath(root), "project_count" => length(projects),
        "projects" => projects)
end

function automate end

function automate(::PSCAD, ::Val{action}, options) where {action}
    throw(ErrorException(
        "PSCAD live automation requires PythonCall and the external mhi.pscad package; " *
        "load PythonCall before running `linecablebenchmark $action`"
    ))
end
