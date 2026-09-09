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
            "formulation", "selection", "frequencies", "basis", "domain", "port_order",
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
    metadata = (path, sha256 = digest, data_sha256 = data_digest,
        case_id = string(document["case_id"]),
        backend = string(document["backend"]), selection = selection,
        axes=coordinates.axes,
        formulation = document["formulation"], implementation = get(document, "implementation", nothing),
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

function benchmark_comparisons(
        settings::NamedTuple, reference::AbstractCoreResult, candidate::AbstractCoreResult)
    settings.statistics == (:value,) || throw(ArgumentError("line parameters require value comparisons"))
    return [(quantity,
                statistic = :value, reference_index = 1, candidate_index = 1,
                error = compare(reference, candidate,
                    getproperty(LineCableModels, quantity); band, normalization,
                    atol = settings.atol, fundamental = settings.fundamental, harmonics = settings.harmonics,
                    unsupported = settings.unsupported)) for band in settings.bands
            for quantity in settings.quantities for normalization in settings.normalizations]
end

function benchmark_comparisons(settings::NamedTuple,
        reference::AbstractCoreResult, candidate::AbstractParametricResult)
    return [merge(row, (reference_index = 1, candidate_index = index))
            for (index, value) in enumerate(candidate)
            for row in benchmark_comparisons(settings, reference, value)]
end

function benchmark_comparisons(settings::NamedTuple,
        reference::AbstractParametricResult, candidate::AbstractCoreResult)
    return [merge(row, (reference_index = index, candidate_index = 1))
            for (index, value) in enumerate(reference)
            for row in benchmark_comparisons(settings, value, candidate)]
end

function benchmark_comparisons(settings::NamedTuple,
        reference::AbstractParametricResult, candidate::AbstractParametricResult)
    length(reference) == length(candidate) ||
        throw(DimensionMismatch("paired result spaces must have equal cardinality"))
    length(NamedTuple(reference).axes.problems) == length(NamedTuple(candidate).axes.problems) ||
        throw(DimensionMismatch("paired result spaces must have equal problem-axis lengths"))
    return [merge(row, (reference_index = index, candidate_index = index))
            for (index, (left, right)) in enumerate(zip(reference, candidate))
            for row in benchmark_comparisons(settings, left, right)]
end

function benchmark_comparisons(settings::NamedTuple, reference::MomentResult, candidate::MomentResult)
    settings.statistics == (:mean, :std) || throw(ArgumentError("moments require mean/std comparisons"))
    comparison = compare(reference, candidate)
    return [(quantity, statistic, reference_index = 1, candidate_index = 1,
                error = getproperty(getproperty(comparison.errors, quantity), statistic))
            for quantity in keys(comparison.errors) for statistic in (:mean, :std)]
end

"Persist comparisons of saved operands; no numerical execution or CI-reference approval is performed."
function compare_saved(benchmark::BenchmarkDefinition; directory::AbstractString)
    validate(benchmark)
    comparison=benchmark_comparisons(benchmark.comparison_settings,
        benchmark.reference.problem.result, benchmark.candidate.problem.result)
    return record_benchmark(benchmark, comparison; directory)
end

"""
    record_benchmark(benchmark, comparisons; directory)

Write completed comparisons with their exact settings and checksummed operands.
No comparison or solver is executed. Existing identical records are verified and reused.
"""
function record_benchmark(benchmark::BenchmarkDefinition, errors; directory::AbstractString)
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
        JLD2.jldsave(temporary; schema_version = 1, kind = :gauntlet_benchmark,
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
and `plot(benchmark, requests; ...)` for REPL tables and matrix-cell overlays.
"""
function read_benchmark(path::AbstractString; load_results::Bool = false)
    if isdir(path)
        reference=read_calculation(joinpath(path, "reference", "calculation.jld2"))
        candidate=read_calculation(joinpath(path, "candidate", "calculation.jld2"))
        analyses=[read_benchmark(joinpath(folder, name))
                  for (folder, _, names) in
                      walkdir(joinpath(path, "analyses"))
                  for name in sort(names) if name == "snapshot.jld2"]
        isempty(analyses) &&
            throw(ArgumentError("benchmark has no retained analysis: $path"))
        return (id = Symbol(basename(path)), reference, candidate, analyses)
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
    record["schema_version"] == 1 && record["kind"] === :gauntlet_benchmark ||
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
        isempty(setdiff(keys(settings), ("kind", "quantities", "bands", "normalizations",
            "atol", "fundamental", "harmonics", "unsupported"))) ||
            throw(ArgumentError("unknown comparison settings"))
        atol=get(settings, "atol", nothing)
        comparison=(
            quantities=Tuple(Symbol.(get(settings, "quantities", kind == "uq_moments" ? ["R", "L", "C", "G"] : ["Z", "Y"]))),
            statistics=kind == "uq_moments" ? (:mean, :std) : (:value,),
            bands=Tuple(b isa String ? Symbol(b) : Tuple(b) for b in get(settings, "bands", ["all"])),
            normalizations=Tuple(Symbol.(get(settings, "normalizations", ["reference_rms"]))),
            atol=atol isa AbstractDict ? (; (Symbol(k)=>v for (k,v) in atol)...) : atol,
            fundamental=get(settings, "fundamental", 50.0), harmonics=get(settings, "harmonics", 50),
            unsupported=(; (Symbol(k)=>v for (k,v) in get(settings, "unsupported", Dict()))...))
        model = (id = Symbol(entry["case"]),
            description = get(entry, "description", entry["case"]))
        benchmark_definition(
            Symbol(entry["id"]), model.id, Symbol(get(plan, "collection", "manual")),
            source, model, calculations..., comparison, get(entry, "tolerances", (;)))
    end
    foreach(validate, benchmarks)
    return [compare_saved(benchmark; directory) for benchmark in benchmarks]
end
