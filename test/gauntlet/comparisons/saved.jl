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
            "applied_variation", "correlation")
        Dict(name => file[name] for name in names if haskey(file, name))
    end
    document["schema_version"] == 1 && document["status"] === :complete ||
        throw(ArgumentError("unsupported or incomplete calculation: $path"))
    kind = document["kind"]
    result = if kind === :gauntlet_calculation
        document["domain"] === :PhaseDomain ||
            throw(ArgumentError("only phase-domain calculations are supported"))
        LineParameters(LineCableModels.PhaseDomain, document["Z"], document["Y"],
            document["frequencies"]; basis = document["basis"],
            details = (comparison_unsupported = get(document, "comparison_unsupported", (;)),))
    elseif kind === :gauntlet_moments
        MomentResult(document["moments"])
    else
        throw(ArgumentError("not a saved Gauntlet calculation: $path"))
    end
    ports = String.(document["port_order"])
    allunique(ports) &&
    length(ports) == (kind === :gauntlet_calculation ?
     size(result.Z, 1) : length(result.port_order)) ||
        throw(DimensionMismatch("invalid saved terminal order: $path"))
    data = kind === :gauntlet_calculation ?
           (Z = vec(document["Z"]), Y = vec(document["Y"]),
        frequencies = document["frequencies"],
        port_order = ports, basis = document["basis"]) : document["moments"]
    data_digest = semantic_sha256(data)
    get(document, "data_sha256", data_digest) == data_digest ||
        throw(ArgumentError("calculation data digest mismatch: $path"))
    metadata = (path, sha256 = digest, data_sha256 = data_digest,
        case_id = string(document["case_id"]),
        backend = string(document["backend"]), selection = document["selection"],
        formulation = document["formulation"], implementation = get(document, "implementation", nothing),
        computation_signature = get(document, "computation_signature", nothing),
        input_sha256 = semantic_sha256(document["problem"]), port_order = ports,
        frequencies = document["frequencies"], basis = document["basis"], domain = document["domain"],
        propagation = get(document, "propagation", nothing), sampling = get(document, "sampling", nothing),
        uncertainty = (parameters = get(document, "parameter_manifest", nothing),
            variation = get(document, "applied_variation", nothing), correlation = get(
                document, "correlation", nothing)),
        timing = (scope = :batch_elapsed_at_completion,
            seconds = get(document, "elapsed_at_completion_seconds", missing),
            batch_selections = get(document, "batch_selection_count", missing)))
    return (; result, metadata)
end

function benchmark_comparisons(
        policy::LineParametersPolicy, reference::LineParameters, candidate::LineParameters)
    return [(quantity,
                statistic = :value,
                error = compare(reference, candidate,
                    getproperty(LineCableModels, quantity); band, normalization,
                    atol = policy.atol, fundamental = policy.fundamental, harmonics = policy.harmonics,
                    unsupported = policy.unsupported)) for band in policy.bands
            for quantity in policy.quantities for normalization in policy.normalizations]
end

function benchmark_comparisons(::UQMomentPolicy, reference::MomentResult, candidate::MomentResult)
    comparison = compare(reference, candidate)
    return [(quantity, statistic,
                error = getproperty(getproperty(comparison.errors, quantity), statistic))
            for quantity in keys(comparison.errors) for statistic in (:mean, :std)]
end

function comparison_policy_record(policy::LineParametersPolicy)
    return (kind = :line_parameters, quantities = policy.quantities, bands = policy.bands,
        normalizations = policy.normalizations, atol = policy.atol,
        fundamental = policy.fundamental, harmonics = policy.harmonics, unsupported = policy.unsupported)
end
comparison_policy_record(::UQMomentPolicy) = (kind = :uq_moments,)

"Persist comparisons of saved operands; no numerical execution or CI-reference approval is performed."
function compare_saved(benchmark::OwnedBenchmark; directory::AbstractString)
    bytes2hex(open(sha256, benchmark.source_file)) == benchmark.source_sha256 ||
        throw(ArgumentError("benchmark definition changed after loading"))
    reference = benchmark.reference.problem
    candidate = benchmark.candidate.problem
    a, b = reference.metadata, candidate.metadata
    a.case_id == b.case_id == string(benchmark.case_id) &&
    a.input_sha256 == b.input_sha256 && a.port_order == b.port_order &&
    a.frequencies == b.frequencies && a.basis == b.basis && a.domain == b.domain ||
        throw(ArgumentError("benchmark $(benchmark.id): saved inputs or coordinates differ; no implicit conversion is permitted"))
    errors = benchmark_comparisons(benchmark.comparison_policy, reference.result, candidate.result)
    # Persist the existing RMSError payload as plain arrays/details. Readers do
    # not depend on the concrete parametric type used when it was calculated.
    comparisons = [(quantity = row.quantity, statistic = row.statistic,
                       absolute = row.error.absolute, relative = row.error.relative,
                       details = isempty(row.error.details) ?
                                 (normalization = :reference_rms, band = :all,
                           actual_bounds = extrema(a.frequencies),
                           sample_count = length(a.frequencies), indices = 1:length(a.frequencies),
                           status = fill(:compared, size(row.error.relative)), reason = nothing, atol = nothing) :
                                 row.error.details) for row in errors]
    for operand in (a, b)
        bytes2hex(open(sha256, operand.path)) == operand.sha256 ||
            throw(ArgumentError("saved operand changed during comparison: $(operand.path)"))
    end
    path = joinpath(abspath(directory), string(benchmark.id), "snapshot.jld2")
    ispath(path) &&
        throw(ArgumentError("benchmark record already exists: $path; select a new output directory"))
    mkpath(dirname(path))
    temporary = tempname(dirname(path))
    try
        JLD2.jldsave(temporary; schema_version = 1, kind = :gauntlet_benchmark,
            numerical_reference_approval = :unreviewed, benchmark_id = string(benchmark.id),
            case_id = string(benchmark.case_id), description = benchmark.model.description,
            collection = string(benchmark.collection), benchmark_source_sha256 = benchmark.source_sha256,
            calculations = (
                reference = merge(a, (
                    id = benchmark.reference.id, owner = benchmark.reference.owner)),
                candidate = merge(b, (
                    id = benchmark.candidate.id, owner = benchmark.candidate.owner))),
            comparison_policy = comparison_policy_record(benchmark.comparison_policy),
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

"Read an explicitly ordered benchmark, checking its record and both numerical operand files."
function read_benchmark(path::AbstractString)
    isfile(path) && isfile(path * ".sha256") ||
        throw(ArgumentError("benchmark record and checksum required: $path"))
    checksum = split(read(path * ".sha256", String))
    !isempty(checksum) && first(checksum) == bytes2hex(open(sha256, path)) ||
        throw(ArgumentError("benchmark checksum mismatch: $path"))
    record = JLD2.load(path)
    record["schema_version"] == 1 && record["kind"] === :gauntlet_benchmark ||
        throw(ArgumentError("explicit benchmark record required; calculation artifacts do not declare references: $path"))
    operands = record["calculations"]
    keys(operands) == (:reference, :candidate) ||
        throw(ArgumentError("benchmark requires one explicit reference and candidate"))
    for operand in operands
        source = isabspath(operand.path) ? operand.path :
                 normpath(joinpath(dirname(path), operand.path))
        isfile(source) && bytes2hex(open(sha256, source)) == operand.sha256 ||
            throw(ArgumentError("benchmark operand missing or checksum mismatch: $source"))
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
            benchmark_calculation(role, Symbol(binding["owner"]), value, value.metadata.formulation)
        end
        settings = get(entry, "comparison", get(plan, "comparison", Dict{String, Any}()))
        policy = if get(settings, "kind", "line_parameters") == "uq_moments"
            isempty(setdiff(keys(settings), ("kind",))) || throw(ArgumentError(
                "UQMomentPolicy retains full-band R/L/C/G mean/std reference-RMS comparisons; custom policy settings are not supported"))
            UQMomentPolicy()
        else
            get(settings, "kind", "line_parameters") == "line_parameters" ||
                throw(ArgumentError("unknown comparison policy"))
            isempty(setdiff(keys(settings),
                ("kind", "quantities", "bands", "normalizations",
                    "atol", "fundamental", "harmonics", "unsupported"))) ||
                throw(ArgumentError("unknown line-parameter comparison settings"))
            tolerances = get(settings, "atol", nothing)
            unsupported = get(settings, "unsupported", Dict{String, Any}())
            LineParametersPolicy(
                quantities = Symbol.(get(settings, "quantities", ["Z", "Y"])),
                bands = map(b -> b isa String ? Symbol(b) : Tuple(b), get(settings, "bands", ["all"])),
                normalizations = Symbol.(get(settings, "normalizations", ["reference_rms"])),
                atol = tolerances isa Dict ?
                       (; (Symbol(k)=>v for (k, v) in tolerances)...) : tolerances,
                fundamental = get(settings, "fundamental", 50.0), harmonics = get(settings, "harmonics", 50),
                unsupported = (; (Symbol(k)=>v for (k, v) in unsupported)...))
        end
        model = (id = Symbol(entry["case"]),
            description = get(entry, "description", entry["case"]))
        benchmark_definition(
            Symbol(entry["id"]), model.id, Symbol(get(plan, "collection", "manual")),
            source, model, calculations..., policy, get(entry, "tolerances", (;)))
    end
    return [compare_saved(benchmark; directory) for benchmark in benchmarks]
end
