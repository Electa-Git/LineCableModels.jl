const DEFAULT_REPORT_FLOOR = (Z = 1.0e-6, Y = 1.0e-9)

function _snapshot_digest(path::AbstractString)
    return bytes2hex(sha256(read(path)))
end

function _snapshot_case(path::AbstractString, backend::Symbol)
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
        "gauntlet_version", "case_name", "backend", "port_order", "frequencies",
        "formulation", "reference_execution", "reference", "accepted",
        "reference_comparison", "julia_benchmark", "recorded_at_utc"
    )
    missing = filter(key -> !haskey(snapshot, key), required)
    isempty(missing) || throw(ArgumentError(
        "Gauntlet snapshot $path is missing fields: $(join(missing, ", "))",
    ))
    snapshot["gauntlet_version"] == string(GAUNTLET_VERSION) || throw(ArgumentError(
        "Gauntlet snapshot $path uses version $(snapshot["gauntlet_version"]), " *
        "expected $(GAUNTLET_VERSION)",
    ))
    snapshot["backend"] == string(backend) || throw(ArgumentError(
        "Gauntlet snapshot $path belongs to backend $(snapshot["backend"]), " *
        "not $backend",
    ))
    snapshot["case_name"] == basename(dirname(path)) || throw(ArgumentError(
        "Gauntlet snapshot $path does not match its case directory",
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
    for quantity in (:Z, :Y)
        stored_error = getproperty(comparison, quantity)
        observed_error = getproperty(observed, quantity)
        isequal(stored_error.absolute, observed_error.absolute) &&
        isequal(stored_error.relative, observed_error.relative) ||
            throw(ArgumentError(
                "stored $quantity comparison does not match the recorded results in $path",
            ))
    end

    frequencies_value = snapshot["frequencies"]
    frequencies_value == reference.f || throw(ArgumentError(
        "Gauntlet snapshot frequencies do not match its reference result: $path",
    ))
    frequencies_value == accepted.f || throw(ArgumentError(
        "Gauntlet snapshot frequencies do not match its accepted result: $path",
    ))
    port_order = snapshot["port_order"]
    length(port_order) == size(reference.Z.values, 1) || throw(ArgumentError(
        "Gauntlet snapshot port order does not match its matrix dimensions: $path",
    ))

    execution = snapshot["reference_execution"]
    execution.backend === backend || throw(ArgumentError(
        "Gauntlet snapshot execution backend does not match $backend: $path",
    ))
    hasproperty(execution, :exit_code) || throw(ArgumentError(
        "Gauntlet snapshot has no reference exit code: $path",
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
    return (; snapshot, comparison, digest = observed_digest)
end

function _maximum(values::AbstractMatrix)
    any(isnan, values) && throw(ArgumentError("RMS statistics cannot contain NaN"))
    index = argmax(values)
    return (value = values[index], row = index[1], column = index[2])
end

function _visible_maximum(error::RMSError, floor::Real)
    indices = findall(error.absolute .> floor)
    isempty(indices) && return (value = missing, row = missing, column = missing)
    index = indices[argmax(error.relative[indices])]
    return (value = error.relative[index], row = index[1], column = index[2])
end

function _error_summary(error::RMSError, floor::Real)
    absolute = _maximum(error.absolute)
    relative_raw = _maximum(error.relative)
    relative = _visible_maximum(error, floor)
    return (; absolute, relative, relative_raw)
end

function _method_description(formulation, owner::Symbol, field::Symbol)
    hasproperty(formulation, owner) || return missing
    methods = getproperty(formulation, owner)
    hasproperty(methods, field) || return missing
    method = getproperty(methods, field)
    hasproperty(method, :description) || return missing
    return String(method.description)
end

function _report_row(path::AbstractString, backend::Symbol, zero_atol::NamedTuple)
    loaded = _snapshot_case(path, backend)
    snapshot = loaded.snapshot
    comparison = loaded.comparison
    impedance = _error_summary(comparison.Z, zero_atol.Z)
    admittance = _error_summary(comparison.Y, zero_atol.Y)
    reference = snapshot["reference"]
    execution = snapshot["reference_execution"]
    benchmark = snapshot["julia_benchmark"]
    formulation = snapshot["formulation"]
    environment = benchmark.environment
    reference_version = hasproperty(execution, :pscad_version) ?
                        execution.pscad_version :
                        hasproperty(execution, :version) ? execution.version : missing
    return (
        case = snapshot["case_name"],
        conductors = size(reference.Z.values, 1),
        frequencies = length(reference.f),
        frequency_start_hz = first(reference.f),
        frequency_end_hz = last(reference.f),
        basis = string(basis(reference)),
        domain = string(nameof(domain(reference))),
        reference_earth = _method_description(formulation, :reference, :earth_impedance),
        linecablemodels_earth = _method_description(
            formulation,
            :candidate,
            :earth_impedance
        ),
        linecablemodels_insulation = _method_description(
            formulation,
            :candidate,
            :insulation_admittance
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

"""
    report(path::AbstractString; backend, zero_atol=DEFAULT_REPORT_FLOOR)
    report(backend::Symbol; artifact_root=ARTIFACT_ROOT, zero_atol=DEFAULT_REPORT_FLOOR)

Read a recorded Gauntlet backend collection and return one aggregate row per case.

The report verifies every snapshot digest, successful reference exit status, frequency and port dimensions, and stored `LineParametersBenchmark`. Relative RMS display values exclude matrix terms whose absolute RMS is at or below the corresponding `zero_atol` floor. Raw relative maxima remain in the `Z_rms_relative_raw` and `Y_rms_relative_raw` columns.

# Arguments

- `path`: Collection directory containing `cases/`.
- `backend`: Recorded reference backend.

# Keywords

- `artifact_root`: Root of ignored staged collections.
- `zero_atol`: Nonnegative absolute RMS display floors for Z \\[Ω/m or Ω\\] and Y \\[S/m or S\\], according to the stored result basis.

# Returns

- A `DataFrame` containing numerical, structural, execution, performance, and integrity statistics for every case.

# Errors

- `ArgumentError` when the collection is empty or any snapshot, digest, result, comparison, dimension, backend, execution status, or display floor is invalid.
"""
function report(
        path::AbstractString;
        backend::Symbol,
        zero_atol::NamedTuple = DEFAULT_REPORT_FLOOR
)
    _report_floor(zero_atol)
    cases = joinpath(path, "cases")
    isdir(cases) ||
        throw(ArgumentError("Gauntlet collection has no cases directory: $path"))
    directories = sort!(filter(isdir, readdir(cases; join = true)))
    isempty(directories) && throw(ArgumentError(
        "Gauntlet collection has no recorded cases: $path",
    ))
    rows = [_report_row(joinpath(directory, "snapshot.jld2"), backend, zero_atol)
            for directory in directories]
    frame = DataFrame(rows)
    metadata!(frame, "backend", backend, style = :note)
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
        format_version = 1,
        backend = string(backend),
        gauntlet_version = string(GAUNTLET_VERSION),
        zero_atol,
        report = frame,
        validation = (
            snapshot_digests = true,
            reference_exit_codes = true,
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
