const DEFAULT_REPORT_FLOOR = (Z = 1.0e-6, Y = 1.0e-9)

_snapshot_digest(path::AbstractString) = bytes2hex(sha256(read(path)))

function _snapshot_document(path::AbstractString, collection::Symbol)
    isfile(path) || throw(ArgumentError("Gauntlet snapshot is missing: $path"))
    digest_path = joinpath(dirname(path), "snapshot.sha256")
    isfile(digest_path) || throw(ArgumentError(
        "Gauntlet snapshot digest is missing: $digest_path",
    ))
    expected_digest = first(split(strip(read(digest_path, String))))
    observed_digest = _snapshot_digest(path)
    expected_digest == observed_digest || throw(ArgumentError(
        "Gauntlet snapshot SHA-256 does not match $digest_path",
    ))
    snapshot = JLD2.load(path)
    required = (
        "schema_version", "gauntlet_version", "benchmark_id", "case_id",
        "collection", "case_source_sha256", "benchmark_source_sha256",
        "parameter_manifest", "applied_variation", "correlation",
        "calculations", "comparison_policy", "tolerances", "port_order", "frequencies",
        "reference_comparison", "recorded_at_utc"
    )
    missing = filter(key -> !haskey(snapshot, key), required)
    isempty(missing) || throw(ArgumentError(
        "Gauntlet snapshot $path is missing fields: $(join(missing, ", "))",
    ))
    snapshot["schema_version"] == 2 || throw(ArgumentError(
        "Gauntlet snapshot $path does not use schema 2",
    ))
    snapshot["gauntlet_version"] == string(GAUNTLET_VERSION) || throw(ArgumentError(
        "Gauntlet snapshot $path uses version $(snapshot["gauntlet_version"]), " *
        "expected $(GAUNTLET_VERSION)",
    ))
    snapshot["collection"] == string(collection) || throw(ArgumentError(
        "Gauntlet snapshot $path belongs to collection $(snapshot["collection"]), " *
        "not $collection",
    ))
    snapshot["benchmark_id"] == basename(dirname(path)) || throw(ArgumentError(
        "Gauntlet snapshot $path does not match its benchmark directory",
    ))
    return (; snapshot, digest = observed_digest)
end

function _primitive_values(parameters::LineParameters)
    return (
        frequency = observe(parameters, frequencies),
        series_impedance = observe(parameters, Z),
        shunt_admittance = observe(parameters, Y)
    )
end

function _comparison_values(comparison::LineParametersBenchmark)
    return (
        series_impedance_absolute_error = observe(comparison, Z, absolute_error),
        series_impedance_relative_error = observe(comparison, Z, relative_error),
        shunt_admittance_absolute_error = observe(comparison, Y, absolute_error),
        shunt_admittance_relative_error = observe(comparison, Y, relative_error)
    )
end

function _maximum(values::AbstractMatrix)
    any(isnan, values) && throw(ArgumentError("RMS statistics cannot contain NaN"))
    index = argmax(values)
    return (value = values[index], row = index[1], column = index[2])
end

function _visible_maximum(absolute_values, relative_values, floor::Real)
    indices = findall(absolute_values .> floor)
    isempty(indices) && return (value = missing, row = missing, column = missing)
    index = indices[argmax(relative_values[indices])]
    return (value = relative_values[index], row = index[1], column = index[2])
end

function _error_summary(absolute_values, relative_values, floor::Real)
    return (
        absolute = _maximum(absolute_values),
        relative = _visible_maximum(absolute_values, relative_values, floor),
        relative_raw = _maximum(relative_values)
    )
end

function _method_description(formulation, owner::Symbol, field::Symbol)
    hasproperty(formulation, owner) || return missing
    methods = getproperty(formulation, owner)
    hasproperty(methods, field) || return missing
    method = getproperty(methods, field)
    hasproperty(method, :description) || return missing
    return String(method.description)
end

function _line_snapshot(path::AbstractString, collection::Symbol)
    loaded = _snapshot_document(path, collection)
    snapshot = loaded.snapshot
    snapshot["comparison_policy"] === :line_parameters || throw(ArgumentError(
        "Gauntlet snapshot $path is not a line-parameter comparison",
    ))
    required = (
        "reference_execution", "reference", "accepted", "julia_benchmark",
        "formulation"
    )
    missing = filter(key -> !haskey(snapshot, key), required)
    isempty(missing) || throw(ArgumentError(
        "Gauntlet line snapshot $path is missing fields: $(join(missing, ", "))",
    ))
    reference = snapshot["reference"]
    accepted = snapshot["accepted"]
    comparison = snapshot["reference_comparison"]
    reference isa LineParameters || throw(ArgumentError(
        "Gauntlet snapshot $path has no reference LineParameters",
    ))
    accepted isa LineParameters || throw(ArgumentError(
        "Gauntlet snapshot $path has no accepted LineParameters",
    ))
    comparison isa LineParametersBenchmark || throw(ArgumentError(
        "Gauntlet snapshot $path has no LineParametersBenchmark",
    ))
    observed = compare(reference, accepted)
    isequal(_comparison_values(comparison), _comparison_values(observed)) ||
        throw(ArgumentError(
            "stored comparison does not match the recorded results in $path",
        ))
    reference_values = _primitive_values(reference)
    accepted_values = _primitive_values(accepted)
    snapshot["frequencies"] == reference_values.frequency &&
        snapshot["frequencies"] == accepted_values.frequency || throw(ArgumentError(
            "Gauntlet snapshot frequencies do not match its results: $path",
        ))
    length(snapshot["port_order"]) == size(reference_values.series_impedance, 1) ||
        throw(ArgumentError(
            "Gauntlet snapshot port order does not match matrix dimensions: $path",
        ))
    execution = snapshot["reference_execution"]
    execution.backend === collection || throw(ArgumentError(
        "Gauntlet snapshot execution backend does not match $collection: $path",
    ))
    iszero(execution.exit_code) || throw(ArgumentError(
        "Gauntlet snapshot records reference exit code $(execution.exit_code): $path",
    ))
    path_fields = filter(
        field -> endswith(string(field), "_path") || field === :output_dir,
        propertynames(execution)
    )
    isempty(path_fields) || throw(ArgumentError(
        "Gauntlet snapshot retains local diagnostic paths: $(join(path_fields, ", "))",
    ))
    return (; loaded..., comparison, reference, reference_values)
end

function _line_report_row(path::AbstractString, collection::Symbol, zero_atol)
    loaded = _line_snapshot(path, collection)
    snapshot = loaded.snapshot
    comparison = _comparison_values(loaded.comparison)
    impedance = _error_summary(
        comparison.series_impedance_absolute_error,
        comparison.series_impedance_relative_error,
        zero_atol.Z
    )
    admittance = _error_summary(
        comparison.shunt_admittance_absolute_error,
        comparison.shunt_admittance_relative_error,
        zero_atol.Y
    )
    execution = snapshot["reference_execution"]
    benchmark = snapshot["julia_benchmark"]
    formulation = snapshot["formulation"]
    environment = benchmark.environment
    reference_version = hasproperty(execution, :version) ? execution.version : missing
    return (
        benchmark = snapshot["benchmark_id"],
        case = snapshot["case_id"],
        conductors = size(loaded.reference_values.series_impedance, 1),
        frequencies = length(loaded.reference_values.frequency),
        frequency_start_hz = first(loaded.reference_values.frequency),
        frequency_end_hz = last(loaded.reference_values.frequency),
        basis = string(basis(loaded.reference)),
        domain = string(nameof(domain(loaded.reference))),
        reference_earth = _method_description(formulation, :reference, :earth_impedance),
        linecablemodels_earth = _method_description(
            formulation, :candidate, :earth_impedance
        ),
        linecablemodels_insulation = _method_description(
            formulation, :candidate, :insulation_admittance
        ),
        Z_zero_atol = zero_atol.Z,
        Z_rms_absolute = impedance.absolute.value,
        Z_absolute_row = impedance.absolute.row,
        Z_absolute_column = impedance.absolute.column,
        Z_rms_relative = impedance.relative.value,
        Z_relative_row = impedance.relative.row,
        Z_relative_column = impedance.relative.column,
        Z_rms_relative_raw = impedance.relative_raw.value,
        Y_zero_atol = zero_atol.Y,
        Y_rms_absolute = admittance.absolute.value,
        Y_absolute_row = admittance.absolute.row,
        Y_absolute_column = admittance.absolute.column,
        Y_rms_relative = admittance.relative.value,
        Y_relative_row = admittance.relative.row,
        Y_relative_column = admittance.relative.column,
        Y_rms_relative_raw = admittance.relative_raw.value,
        julia_minimum_ms = 1.0e3 * benchmark.minimum_seconds,
        julia_median_ms = 1.0e3 * benchmark.median_seconds,
        julia_bytes = benchmark.bytes,
        julia_allocations = benchmark.allocations,
        julia_samples = benchmark.samples,
        reference_seconds = execution.elapsed_seconds,
        reference_version,
        julia_version = environment.julia_version,
        kernel = environment.kernel,
        architecture = environment.architecture,
        threads = environment.threads,
        snapshot_sha256 = loaded.digest,
        recorded_at_utc = snapshot["recorded_at_utc"]
    )
end

function _moment_rms(reference::AbstractArray{<:Real, 3}, candidate::AbstractArray{<:Real, 3})
    absolute = Matrix{Float64}(undef, size(reference, 1), size(reference, 2))
    relative = similar(absolute)
    for row in axes(reference, 1), column in axes(reference, 2)
        left = @view reference[row, column, :]
        right = @view candidate[row, column, :]
        difference_norm = sum(abs2, left .- right)
        reference_norm = sum(abs2, left)
        absolute[row, column] = sqrt(difference_norm / length(left))
        relative[row, column] = iszero(reference_norm) ?
                                (iszero(difference_norm) ? 0.0 : Inf) :
                                sqrt(difference_norm / reference_norm)
    end
    return (; absolute, relative)
end

function _validate_moment_records(snapshot, path)
    reference = snapshot["accepted_reference"]
    candidate = snapshot["accepted_candidate"]
    for record in (reference, candidate)
        keys(record.values) == (:R, :L, :C, :G) || throw(ArgumentError(
            "UQ moment snapshot $path must contain exactly R, L, C, and G",
        ))
        record.frequencies == snapshot["frequencies"] || throw(ArgumentError(
            "UQ moment snapshot $path has inconsistent frequencies",
        ))
        record.port_order == snapshot["port_order"] || throw(ArgumentError(
            "UQ moment snapshot $path has inconsistent terminal order",
        ))
        record.domain === :PhaseDomain || throw(ArgumentError(
            "UQ moment snapshot $path must use PhaseDomain",
        ))
        expected = (
            length(record.port_order),
            length(record.port_order),
            length(record.frequencies)
        )
        for product in values(record.values)
            size(product.mean) == expected || throw(DimensionMismatch(
                "UQ moment mean in $path does not match $expected",
            ))
            size(product.std) == expected || throw(DimensionMismatch(
                "UQ moment std in $path does not match $expected",
            ))
            all(isfinite, product.mean) || throw(ArgumentError(
                "UQ moment means in $path must be finite",
            ))
            all(value -> isfinite(value) && value >= 0, product.std) ||
                throw(ArgumentError("UQ moment standard deviations in $path are invalid"))
        end
    end
    reference.basis === candidate.basis || throw(ArgumentError(
        "UQ moment bases in $path do not match",
    ))
    errors = map(reference.values, candidate.values) do left, right
        (
            mean = _moment_rms(left.mean, right.mean),
            std = _moment_rms(left.std, right.std)
        )
    end
    stored = snapshot["reference_comparison"]
    for quantity in keys(errors), statistic in (:mean, :std)
        observed = getproperty(getproperty(errors, quantity), statistic)
        accepted = getproperty(getproperty(stored, quantity), statistic)
        isequal(observed.absolute, accepted.absolute) &&
            isequal(observed.relative, accepted.relative) || throw(ArgumentError(
                "stored UQ comparison does not match moment products in $path",
            ))
    end
    return (; reference, candidate, errors)
end

function _uq_report_row(path::AbstractString, collection::Symbol)
    loaded = _snapshot_document(path, collection)
    snapshot = loaded.snapshot
    snapshot["comparison_policy"] === :uq_moments || throw(ArgumentError(
        "Gauntlet snapshot $path is not a UQ moment comparison",
    ))
    required = (
        "accepted_reference", "accepted_candidate", "monte_carlo", "timings",
        "environment", "tolerances"
    )
    missing = filter(key -> !haskey(snapshot, key), required)
    isempty(missing) || throw(ArgumentError(
        "Gauntlet UQ snapshot $path is missing fields: $(join(missing, ", "))",
    ))
    moments = _validate_moment_records(snapshot, path)
    pairs = Pair{Symbol, Any}[]
    failing = Tuple{Symbol, Symbol, Int, Int}[]
    for quantity in (:R, :L, :C, :G), statistic in (:mean, :std)
        error = getproperty(getproperty(moments.errors, quantity), statistic)
        absolute = _maximum(error.absolute)
        relative_raw = _maximum(error.relative)
        limit = getproperty(
            getproperty(snapshot["tolerances"].reference, statistic), quantity
        )
        relative = _visible_maximum(
            error.absolute,
            error.relative,
            limit.absolute
        )
        prefix = "$(quantity)_$(statistic)_rms"
        append!(pairs, [
            Symbol(prefix, "_absolute") => absolute.value,
            Symbol(prefix, "_relative") => relative.value,
            Symbol(prefix, "_row") => relative.row,
            Symbol(prefix, "_column") => relative.column,
            Symbol(prefix, "_relative_raw") => relative_raw.value
        ])
        failed = findall(.!((error.absolute .<= limit.absolute) .|
                            (error.relative .<= limit.relative)))
        append!(failing, ((quantity, statistic, index[1], index[2]) for index in failed))
    end
    monte_carlo = snapshot["monte_carlo"]
    timings = snapshot["timings"]
    hasproperty(timings, :execution) && hasproperty(timings, :benchmark) ||
        throw(ArgumentError(
            "Gauntlet UQ snapshot $path has no execution/timing benchmark record",
        ))
    performance = timings.benchmark
    performance === nothing && throw(ArgumentError(
        "Gauntlet UQ snapshot $path has no LEP/Monte Carlo timing benchmark",
    ))
    fixed = (
        benchmark = snapshot["benchmark_id"],
        case = snapshot["case_id"],
        conductors = length(snapshot["port_order"]),
        frequencies = length(snapshot["frequencies"]),
        frequency_start_hz = first(snapshot["frequencies"]),
        frequency_end_hz = last(snapshot["frequencies"]),
        basis = string(moments.reference.basis),
        domain = string(moments.reference.domain),
        monte_carlo_trials = monte_carlo.trials,
        monte_carlo_seed = string(monte_carlo.seed),
        failing_matrix_locations = isempty(failing) ? Base.missing :
                                   join(string.(failing), ";"),
        reference_seconds = timings.execution.reference,
        candidate_seconds = timings.execution.candidate,
        lep_minimum_seconds = performance.reference.minimum_seconds,
        lep_median_seconds = performance.reference.median_seconds,
        monte_carlo_minimum_seconds = performance.candidate.minimum_seconds,
        monte_carlo_median_seconds = performance.candidate.median_seconds,
        monte_carlo_over_lep_speedup = performance.speedup,
        minimum_required_speedup = performance.settings.minimum_speedup,
        timing_comparable = performance.comparable,
        lep_timing_samples = performance.reference.samples,
        monte_carlo_timing_samples = performance.candidate.samples,
        julia_version = snapshot["environment"].julia_version,
        kernel = snapshot["environment"].kernel,
        architecture = snapshot["environment"].architecture,
        threads = snapshot["environment"].threads,
        snapshot_sha256 = loaded.digest,
        recorded_at_utc = snapshot["recorded_at_utc"]
    )
    return merge(fixed, (; pairs...))
end

function _report_floor(zero_atol::NamedTuple)
    propertynames(zero_atol) == (:Z, :Y) || throw(ArgumentError(
        "zero_atol must contain exactly Z and Y",
    ))
    for quantity in (:Z, :Y)
        floor = getproperty(zero_atol, quantity)
        floor isa Real && isfinite(floor) && floor >= zero(floor) ||
            throw(ArgumentError("zero_atol.$quantity must be finite and nonnegative"))
    end
    return zero_atol
end

"""Read and validate a recorded Gauntlet v2 collection."""
function report(
        path::AbstractString;
        backend::Symbol,
        zero_atol::NamedTuple = DEFAULT_REPORT_FLOOR
)
    _report_floor(zero_atol)
    benchmarks = joinpath(path, "benchmarks")
    isdir(benchmarks) || throw(ArgumentError(
        "Gauntlet collection has no benchmarks directory: $path",
    ))
    directories = sort!(filter(isdir, readdir(benchmarks; join = true)))
    isempty(directories) && throw(ArgumentError(
        "Gauntlet collection has no recorded benchmarks: $path",
    ))
    documents = [JLD2.load(joinpath(directory, "snapshot.jld2"))
                 for directory in directories]
    policies = unique(document["comparison_policy"] for document in documents)
    length(policies) == 1 || throw(ArgumentError(
        "a Gauntlet collection cannot mix comparison policies",
    ))
    rows = if only(policies) === :line_parameters
        [_line_report_row(joinpath(directory, "snapshot.jld2"), backend, zero_atol)
         for directory in directories]
    elseif only(policies) === :uq_moments
        [_uq_report_row(joinpath(directory, "snapshot.jld2"), backend)
         for directory in directories]
    else
        throw(ArgumentError("unsupported Gauntlet comparison policy $(only(policies))"))
    end
    frame = DataFrame(rows)
    metadata!(frame, "collection", backend, style = :note)
    metadata!(frame, "gauntlet_version", GAUNTLET_VERSION, style = :note)
    metadata!(frame, "zero_atol", zero_atol, style = :note)
    return frame
end

function report(
        backend::Symbol;
        artifact_root::AbstractString = ARTIFACT_ROOT,
        zero_atol::NamedTuple = DEFAULT_REPORT_FLOOR
)
    return report(backend_stage(backend; artifact_root); backend, zero_atol)
end

function _report_value(value)
    value === missing && return "missing"
    text = value isa AbstractFloat ? repr(value) : string(value)
    return replace(text, '\t' => ' ', '\n' => ' ', '\r' => ' ')
end

function _write_report(
        backend::Symbol,
        stage::AbstractString;
        zero_atol::NamedTuple = DEFAULT_REPORT_FLOOR
)
    frame = report(stage; backend, zero_atol)
    jld2_path = joinpath(stage, "report.jld2")
    temporary = jld2_path * ".new"
    isfile(temporary) && rm(temporary; force = true)
    JLD2.jldsave(
        temporary;
        format_version = 2,
        collection = string(backend),
        gauntlet_version = string(GAUNTLET_VERSION),
        zero_atol,
        report = frame,
        validation = (
            snapshot_digests = true,
            calculation_exit_codes = true,
            comparisons_recomputed = true
        )
    )
    mv(temporary, jld2_path; force = true)
    tsv_path = joinpath(stage, "report.tsv")
    open(tsv_path, "w") do io
        println(io, join(string.(names(frame)), '\t'))
        for row in eachrow(frame)
            println(io, join(_report_value.(Tuple(row)), '\t'))
        end
    end
    digest_path = joinpath(stage, "report.sha256")
    open(digest_path, "w") do io
        println(io, "$(_snapshot_digest(jld2_path))  report.jld2")
        println(io, "$(_snapshot_digest(tsv_path))  report.tsv")
    end
    return (; frame, jld2_path, tsv_path, digest_path)
end
