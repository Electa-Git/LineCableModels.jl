function _snapshot_document(path::AbstractString, collection::Symbol)
    isfile(path) || throw(ArgumentError("Gauntlet snapshot is missing: $path"))
    digest_path = joinpath(dirname(path), "snapshot.sha256")
    isfile(digest_path) || throw(ArgumentError(
        "Gauntlet snapshot digest is missing: $digest_path",
    ))
    expected_digest = split(strip(read(digest_path, String)))
    observed_digest = bytes2hex(open(sha256, path))
    !isempty(expected_digest) && first(expected_digest) == observed_digest ||
        throw(ArgumentError(
            "Gauntlet snapshot SHA-256 does not match $digest_path",
        ))
    snapshot = JLD2.load(path)
    if haskey(snapshot, "comparison_policy")
        snapshot["comparison_settings"]=pop!(snapshot, "comparison_policy")
    end
    required = (
        "schema_version", "benchmark_id", "case_id",
        "collection", "case_source_sha256", "benchmark_source_sha256",
        "parameter_manifest", "applied_variation", "correlation",
        "calculations", "comparison_settings", "tolerances", "port_order", "frequencies",
        "reference_comparison", "recorded_at_utc"
    )
    missing = filter(key -> !haskey(snapshot, key), required)
    isempty(missing) || throw(ArgumentError(
        "Gauntlet snapshot $path is missing fields: $(join(missing, ", "))",
    ))
    snapshot["schema_version"] == SNAPSHOT_SCHEMA_VERSION || throw(ArgumentError(
        "Gauntlet snapshot $path does not use schema $SNAPSHOT_SCHEMA_VERSION",
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

function _line_snapshot(path::AbstractString, collection::Symbol)
    loaded = _snapshot_document(path, collection)
    snapshot = loaded.snapshot
    snapshot["comparison_settings"] === :line_parameters || throw(ArgumentError(
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

function _validate_moment_records(snapshot, path)
    reference = MomentResult(snapshot["accepted_reference"])
    candidate = MomentResult(snapshot["accepted_candidate"])
    for record in (reference, candidate)
        record.frequencies == snapshot["frequencies"] || throw(ArgumentError(
            "UQ moment snapshot $path has inconsistent frequencies",
        ))
        record.port_order == snapshot["port_order"] || throw(ArgumentError(
            "UQ moment snapshot $path has inconsistent terminal order",
        ))
    end
    errors = compare(reference, candidate).errors
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

"""
    read_moments(path::AbstractString; collection::Symbol=:uq)

Read and verify a stored UQ comparison without loading cases or computing results.
The snapshot checksum, coordinate agreement, and recorded element-wise RMS errors
are checked using the same `MomentResult` comparison as numerical execution.

# Arguments

- `path`: Completed `snapshot.jld2` with its sibling `snapshot.sha256`.

# Keywords

- `collection`: Recorded artifact collection, default `:uq`.

# Returns

- A named tuple containing `reference` and `candidate` moments, their `errors`,
  the stored `snapshot`, and its `digest`. Reading does not approve a reference.
"""
function read_moments(path::AbstractString; collection::Symbol = :uq)
    loaded = _snapshot_document(path, collection)
    get(loaded.snapshot, "comparison_settings", nothing) === :uq_moments ||
        throw(ArgumentError("Gauntlet snapshot $path is not a UQ moment comparison"))
    for field in ("accepted_reference", "accepted_candidate")
        haskey(loaded.snapshot, field) || throw(ArgumentError(
            "Gauntlet UQ snapshot $path is missing $field"))
    end
    return merge(loaded, _validate_moment_records(loaded.snapshot, path))
end

"""
    read_collection(path::AbstractString; collection::Symbol)

Read and validate stored benchmark snapshots without creating reports or running
calculations. Checksums, result coordinates and stored comparisons are checked.

# Arguments

- `path`: Collection directory containing `benchmarks/<id>/snapshot.jld2`.

# Keywords

- `collection`: Expected collection identifier.

# Returns

- Validated snapshot records, in benchmark-directory order.
"""
function read_collection(path::AbstractString; collection::Symbol)
    benchmarks = joinpath(path, "benchmarks")
    isdir(benchmarks) ||
        throw(ArgumentError("collection has no benchmarks directory: $path"))
    directories = sort!(filter(isdir, readdir(benchmarks; join = true)))
    isempty(directories) &&
        throw(ArgumentError("collection has no recorded benchmarks: $path"))
    return map(directories) do directory
        file = joinpath(directory, "snapshot.jld2")
        if !isfile(file)
            record=read_benchmark(directory)
            all(analysis -> analysis["collection"] == string(collection), record.analyses) ||
                throw(ArgumentError("benchmark collection differs: $directory"))
            return record
        end
        settings = jldopen(file, "r") do input
            input[haskey(input, "comparison_settings") ? "comparison_settings" : "comparison_policy"]
        end
        settings === :line_parameters && return _line_snapshot(file, collection)
        settings === :uq_moments && return read_moments(file; collection)
        throw(ArgumentError("unsupported stored comparison settings $settings: $file"))
    end
end
