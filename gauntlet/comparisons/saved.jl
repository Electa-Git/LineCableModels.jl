"Read a checksummed completed calculation without reconstructing a design or loading a solver."
function read_calculation(path::AbstractString; sha256_expected = nothing)
    path = abspath(path)
    isfile(path) && isfile(path * ".sha256") || throw(ArgumentError(
        "completed calculation and checksum required: $path"))
    digest = bytes2hex(open(sha256, path))
    checksum = split(read(path * ".sha256", String))
    !isempty(checksum) && first(checksum) == digest &&
    (sha256_expected === nothing || sha256_expected == digest) || throw(ArgumentError(
        "calculation checksum mismatch: $path"))
    # Read plain numerical and selection records only. In particular, do not
    # deserialize backend-specific computation_details from older implementations.
    document = jldopen(path, "r") do file
        names = ("schema_version", "kind", "status", "case_id", "backend", "problem",
            "formulation", "calculation", "repository", "active_project", "selection", "frequencies", "basis", "domain", "port_order",
            "Z", "Y", "moments", "comparison_unsupported", "data_sha256",
            "implementation", "computation_signature", "elapsed_at_completion_seconds",
            "batch_selection_count", "propagation", "sampling", "parameter_manifest",
            "applied_variation", "correlation", "computation_details", "retained_files", "points", "axes")
        Dict(name => file[name] for name in names if haskey(file, name))
    end
    document["schema_version"] in (1, 2) && document["status"] === :complete ||
        throw(ArgumentError("unsupported or incomplete calculation: $path"))
    for file in get(document, "retained_files", ())
        (isabspath(file.path) || first(splitpath(normpath(file.path))) == "..") &&
            throw(ArgumentError("retained evidence must be inside its calculation directory"))
        evidence=joinpath(dirname(path), file.path)
        isfile(evidence) && bytes2hex(open(sha256, evidence)) == file.sha256 ||
            throw(ArgumentError("retained backend evidence is missing or changed: $evidence"))
    end
    kind = document["kind"]
    result = if kind === :gauntlet_calculation
        document["domain"] === :PhaseDomain ||
            throw(ArgumentError("only phase-domain calculations are supported"))
        LineParameters(LineCableModels.PhaseDomain, document["Z"], document["Y"],
            document["frequencies"]; basis = document["basis"],
            details = merge(get(document, "computation_details", (;)),
                (comparison_unsupported = get(document, "comparison_unsupported", (;)),
                    files = [(path = file.path,
                                 source = joinpath(dirname(path), file.path),
                                 sha256 = file.sha256)
                             for file in get(document, "retained_files", ())])))
    elseif kind === :gauntlet_result_space
        values=[LineParameters(
                    PhaseDomain, point.Z, point.Y, point.frequencies; basis = point.basis,
                    details = merge(point.computation_details,
                        (files = [(path = file.path,
                                      source = joinpath(dirname(path), file.path), sha256 = file.sha256)
                                  for file in get(document, "retained_files", ())
                                  if file.point == index],)))
                for (index, point) in enumerate(document["points"])]
        for (point, value) in zip(document["points"], values)
            semantic_sha256(value, (; port_order=document["port_order"])) == point.data_sha256 ||
                throw(ArgumentError("result-space point data changed"))
        end
        ParametricResult(nothing, values, document["axes"], (;))
    elseif kind === :gauntlet_moments
        MomentResult(document["moments"])
    else
        throw(ArgumentError("not a saved Gauntlet calculation: $path"))
    end
    ports = String.(document["port_order"])
    allunique(ports) &&
    length(ports) == (kind === :gauntlet_calculation ?
     size(result.Z, 1) :
     kind === :gauntlet_result_space ? size(first(result).Z, 1) :
     length(result.port_order)) ||
        throw(DimensionMismatch("invalid saved terminal order: $path"))
    coordinates=(port_order=ports, axes=get(document, "axes", nothing))
    data_digest = semantic_sha256(result, coordinates)
    get(document, "data_sha256", data_digest) == data_digest ||
        throw(ArgumentError("calculation data digest mismatch: $path"))
    selection=document["selection"]
    if selection isa AbstractDict
        selection=(; (Symbol(k)=>v for (k,v) in selection)...)
    elseif selection isa Union{AbstractString, Symbol}
        selection=(id=Symbol(selection),)
    end
    selection isa NamedTuple || throw(ArgumentError("saved selection must be a named record"))
    formulation=document["formulation"]
    if formulation isa NamedTuple && all(key -> haskey(formulation,key),(:id,:input_sha256,:formulation,:options))
        formulation=formulation.formulation
    end
    metadata = (path, sha256 = digest, data_sha256 = data_digest,
        case_id = string(document["case_id"]),
        backend = string(document["backend"]), selection = selection,
        axes=coordinates.axes,
        formulation, calculation=get(document,"calculation",nothing),
        repository=get(document,"repository",nothing),active_project=get(document,"active_project",nothing), implementation = get(document, "implementation", nothing),
        computation_signature = get(document, "computation_signature", nothing),
        input_sha256 = semantic_sha256(document["problem"]), port_order = ports,
        frequencies = document["frequencies"], basis = document["basis"], domain = document["domain"],
        propagation = get(document, "propagation", nothing), sampling = get(document, "sampling", nothing),
        uncertainty = (parameters = get(document, "parameter_manifest", nothing),
            variation = get(document, "applied_variation", nothing), correlation = get(
                document, "correlation", nothing)),
        timing = (
            scope = document["schema_version"] == 1 ? :batch_elapsed_at_completion :
                    :compute_wall,
            seconds = get(document, "elapsed_at_completion_seconds", missing),
            batch_selections = get(document, "batch_selection_count", missing)))
    return (; result, metadata)
end

function validate(::typeof(read_calculation), result::Union{AbstractCoreResult, MomentResult}, metadata::NamedTuple)
    semantic_sha256(result, metadata) == metadata.data_sha256 || throw(ArgumentError(
        "result was modified after loading; saved RMS values no longer describe it"))
    return nothing
end

function validate(::typeof(read_calculation), result::AbstractParametricResult, metadata::NamedTuple)
    NamedTuple(result).axes == metadata.axes &&
        semantic_sha256(result, metadata) == metadata.data_sha256 || throw(ArgumentError(
            "result or axes were modified after loading; saved RMS values no longer describe them"))
    return nothing
end

"Persist comparisons of saved operands; no numerical execution or CI-reference approval is performed."
function compare_saved(benchmark::BenchmarkDefinition; directory::AbstractString)
    validate(Base.write,directory)
    validate(benchmark)
    publication=report(BenchmarkTableDefinition(;benchmark.comparison_settings...),
        (reference=benchmark.reference.problem, candidate=benchmark.candidate.problem,
            context=(id=benchmark.id,case_id=benchmark.case_id,collection=benchmark.collection)))
    return record_benchmark(benchmark, publication; directory)
end

"""
    record_benchmark(benchmark, comparisons; directory)

Write completed comparisons with their exact settings and checksummed operands.
No comparison or solver is executed. Existing identical records are verified and reused.
"""
function record_benchmark(benchmark::BenchmarkDefinition, publication::ReportArtifact; directory::AbstractString)
    validate(Base.write,directory)
    errors=publication.published.comparisons
    bytes2hex(open(sha256, benchmark.source_file)) == benchmark.source_sha256 ||
        throw(ArgumentError("benchmark definition changed after loading"))
    reference = benchmark.reference.problem
    candidate = benchmark.candidate.problem
    a, b = reference.metadata, candidate.metadata
    a.port_order == b.port_order ||
        throw(ArgumentError("benchmark $(benchmark.id): saved inputs or coordinates differ; no implicit conversion is permitted"))
    # Persist the existing RMSError payload as plain arrays/details. Readers do
    # not depend on the concrete parametric type used when it was calculated.
    comparisons = [(quantity = row.quantity, statistic = row.statistic,
                       reference_index = row.reference_index, candidate_index = row.candidate_index,
                       absolute = row.error.absolute, relative = row.error.relative,
                       details = row.error.details) for row in errors]
    for error in errors
        all(name -> haskey(error.error.details, name),
            (:band, :normalization, :actual_bounds, :sample_count, :indices, :status)) ||
            throw(ArgumentError("completed comparisons must retain their actual settings and sample coordinates"))
    end
    for operand in (a, b)
        bytes2hex(open(sha256, operand.path)) == operand.sha256 ||
            throw(ArgumentError("saved operand changed during comparison: $(operand.path)"))
    end
    analysis_id=semantic_sha256((reference = a.sha256, candidate = b.sha256,
        comparison = benchmark.comparison_settings))
    path = joinpath(abspath(directory), string(benchmark.id), analysis_id, "snapshot.jld2")
    if isfile(path)
        read_benchmark(path)
        return path
    end
    ispath(path) &&
        throw(ArgumentError("benchmark record already exists: $path; select a new output directory"))
    mkpath(dirname(path))
    temporary = tempname(dirname(path))
    try
        JLD2.jldsave(temporary; schema_version = 2, kind = :gauntlet_benchmark,
            analysis_id,
            summary=[merge(NamedTuple(row),(snapshot=analysis_id,)) for row in eachrow(publication.table.maxima)],
            formulations=[_selection_value(NamedTuple(row)) for row in eachrow(publication.table.formulations)],
            benchmark_id = string(benchmark.id),
            case_id = string(benchmark.case_id), description = benchmark.model isa
                                                               LoadedCase ?
                                                               benchmark.model.definition.description :
                                                               benchmark.model.description,
            collection = string(benchmark.collection), benchmark_source_sha256 = benchmark.source_sha256,
            calculations = (
                reference = merge(a, (path = relpath(a.path, dirname(path)),
                    id = benchmark.reference.id)),
                candidate = merge(b, (path = relpath(b.path, dirname(path)),
                    id = benchmark.candidate.id))),
            comparison_settings = benchmark.comparison_settings,
            tolerances = benchmark.tolerances, port_order = a.port_order, frequencies = a.frequencies,
            basis = a.basis, domain = a.domain, reference_comparison = comparisons,
            timings = (reference = a.timing, candidate = b.timing), recorded_at_utc = string(now(UTC)))
        mv(temporary, path; force = false)
        write(path * ".sha256", bytes2hex(open(sha256, path)) * "  snapshot.jld2\n")
    finally
        isfile(temporary) && rm(temporary)
    end
    return path
end

"""
    read_benchmark(path; load_results=false)

Read an explicitly ordered benchmark and verify the snapshot and operand hashes.
A benchmark directory returns `(id, reference, candidate, analyses)`, with loaded
results and metadata. For a snapshot file, the default returns the recorded
dictionary; `load_results=true` returns the same loaded bundle as a directory,
resolving operand paths relative to the snapshot. Loading performs no solve.

The loaded bundle supports `report(BenchmarkTableDefinition(false), benchmark)`
and `plot(benchmark, ydata; ...)` for REPL tables and matrix-cell overlays.
"""
function read_benchmark(path::AbstractString; load_results::Bool = false, previous::Bool=false)
    if isdir(path)
        state_path=joinpath(path,"state.toml")
        if isfile(state_path)
            state=TOML.parsefile(state_path)
            if haskey(state,"current") || haskey(state,"attempt")
                (state["state"] == "complete" || previous) || throw(ArgumentError(
                    "latest draft is $(state["state"]); use previous=true to inspect the previous complete result"))
                haskey(state,"current") || throw(ArgumentError("no previous completed result"))
                current=state["current"]
                isabspath(current) || first(splitpath(normpath(current))) == ".." ?
                    throw(ArgumentError("invalid draft path")) : nothing
                return read_benchmark(joinpath(path,current);load_results)
            end
        end
        reference=read_calculation(joinpath(path, "reference", "calculation.jld2"))
        candidate=read_calculation(joinpath(path, "candidate", "calculation.jld2"))
        analyses=[read_benchmark(joinpath(folder, name))
                  for (folder, _, names) in
                      walkdir(joinpath(path, "analyses"))
                  for name in sort(names) if name == "snapshot.jld2"]
        isempty(analyses) &&
            throw(ArgumentError("benchmark has no retained analysis: $path"))
        return (id = Symbol(first(analyses)["benchmark_id"]), reference, candidate, analyses)
    end
    isfile(path) && isfile(path * ".sha256") ||
        throw(ArgumentError("benchmark record and checksum required: $path"))
    checksum = split(read(path * ".sha256", String))
    !isempty(checksum) && first(checksum) == bytes2hex(open(sha256, path)) ||
        throw(ArgumentError("benchmark checksum mismatch: $path"))
    record = JLD2.load(path)
    if haskey(record, "comparison_policy")
        record["comparison_settings"] = pop!(record, "comparison_policy")
    end
    settings=record["comparison_settings"]
    if haskey(settings, :kind)
        settings.kind in (:line_parameters, :uq_moments) || throw(ArgumentError("unsupported saved comparison kind"))
        statistics=settings.kind === :uq_moments ? (:mean, :std) : (:value,)
        record["comparison_settings"]=merge(Base.structdiff(settings, (;kind=settings.kind)), (;statistics))
    end
    record["schema_version"] in (1,2) && record["kind"] === :gauntlet_benchmark ||
        throw(ArgumentError("explicit benchmark record required; calculation artifacts do not declare references: $path"))
    record["reference_comparison"] = [merge((reference_index=1, candidate_index=1), row)
        for row in record["reference_comparison"]]
    operands = record["calculations"]
    keys(operands) == (:reference, :candidate) ||
        throw(ArgumentError("benchmark requires one explicit reference and candidate"))
    for operand in operands
        source = isabspath(operand.path) ? operand.path :
                 normpath(joinpath(dirname(path), operand.path))
        isfile(source) && bytes2hex(open(sha256, source)) == operand.sha256 ||
            throw(ArgumentError("benchmark operand missing or checksum mismatch: $source"))
    end
    if load_results
        reference, candidate = map(operands) do operand
            source = isabspath(operand.path) ? operand.path :
                     normpath(joinpath(dirname(path), operand.path))
            read_calculation(source; sha256_expected = operand.sha256)
        end
        return (id = Symbol(record["benchmark_id"]), reference, candidate, analyses = [record])
    end
    return record
end

"Compare explicit file bindings in a TOML benchmark definition; missing operands are errors, never inferred."
function compare_saved(definition::AbstractString; directory::AbstractString)
    source = abspath(definition)
    plan = TOML.parsefile(source)
    get(plan, "schema_version", nothing) == 1 ||
        throw(ArgumentError("unsupported benchmark definition schema"))
    entries = get(plan, "benchmarks", nothing)
    entries isa Vector && !isempty(entries) ||
        throw(ArgumentError("define at least one explicit benchmark"))
    allunique(entry["id"] for entry in entries) ||
        throw(ArgumentError("duplicate benchmark IDs"))
    # Resolve every declared operand before publishing comparisons.
    benchmarks = map(entries) do entry
        calculations = map((:reference, :candidate)) do role
            haskey(entry, string(role)) ||
                throw(ArgumentError("benchmark $(entry["id"]) needs an explicit $role"))
            binding = entry[string(role)]
            value = read_calculation(normpath(joinpath(dirname(source), binding["path"]));
                sha256_expected = binding["sha256"])
            BenchmarkCalculation(role, value, value.metadata.formulation)
        end
        settings = get(entry, "comparison", get(plan, "comparison", Dict{String, Any}()))
        kind=get(settings, "kind", "line_parameters")
        kind in ("line_parameters", "uq_moments") || throw(ArgumentError("unknown comparison kind"))
        isempty(setdiff(keys(settings), ("kind","quantities","statistics","bands","normalizations",
            "atol","fundamental","harmonics","unsupported","pairing"))) ||
            throw(ArgumentError("unknown comparison settings"))
        # Convert only declared file values. ReportBuilder owns all omitted defaults.
        comparison=(; (Symbol(key)=>
            (key in ("quantities","statistics","normalizations") ? Tuple(Symbol.(value)) :
             key == "bands" ? Tuple(band isa String ? Symbol(band) : Tuple(band) for band in value) :
             key == "pairing" ? Tuple(Tuple(pair) for pair in value) :
             key in ("atol","unsupported") && value isa AbstractDict ? (;(Symbol(k)=>v for (k,v) in value)...) : value)
            for (key,value) in settings if key != "kind")...)
        if kind == "uq_moments"
            get(comparison,:statistics,(:mean,:std)) == (:mean,:std) ||
                throw(ArgumentError("moment kind requires mean/std statistics"))
            comparison=merge(comparison,(statistics=(:mean,:std),))
        end
        model = (id = Symbol(entry["case"]),
            description = get(entry, "description", entry["case"]))
        benchmark_definition(
            Symbol(entry["id"]), model.id, Symbol(get(plan, "collection", "manual")),
            source, model, calculations..., comparison, get(entry, "tolerances", (;)))
    end
    foreach(validate, benchmarks)
    return [compare_saved(benchmark; directory) for benchmark in benchmarks]
end

"""Hash the complete retained operands, comparisons and declaration of one benchmark."""
function semantic_sha256(::typeof(read_benchmark),directory::AbstractString)
    root=abspath(directory)
    read_benchmark(root)
    state_path=joinpath(root,"state.toml")
    if isfile(state_path)
        state=TOML.parsefile(state_path)
        haskey(state,"current") && (root=joinpath(root,state["current"]))
    end
    files=[joinpath(root,role,"calculation.jld2") for role in ("reference","candidate")]
    declaration=joinpath(root,"declarations.jld2")
    isfile(declaration) && push!(files,declaration)
    for (folder,_,names) in walkdir(joinpath(root,"analyses"))
        "snapshot.jld2" in names && push!(files,joinpath(folder,"snapshot.jld2"))
    end
    return semantic_sha256(Dict(relpath(path,root)=>bytes2hex(open(sha256,path)) for path in files))
end
