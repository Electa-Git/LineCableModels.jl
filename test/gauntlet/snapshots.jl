function write_snapshot(
        path::AbstractString,
        case::GauntletCase,
        reference::LineParameters;
        pscad_version::AbstractString,
        pscad_elapsed_seconds::Real,
        julia_benchmark = nothing
)
    validate_structure(case, reference)
    mkpath(dirname(path))
    temporary = path * ".new"
    isfile(temporary) && rm(temporary; force = true)
    try
        JLD2.jldsave(
            temporary;
            format_version = SNAPSHOT_FORMAT_VERSION,
            case_name = string(case.name),
            case_sha256 = case_digest(case),
            line_parameters = reference,
            port_order = case.port_order,
            kron_reduced = case.kron_reduced,
            pscad_version = String(pscad_version),
            pscad_formulation = string(case.pscad_formulation),
            pscad_elapsed_seconds = Float64(pscad_elapsed_seconds),
            julia_benchmark,
            recorded_at_utc = string(now(UTC))
        )
        mv(temporary, path; force = true)
    catch
        isfile(temporary) && rm(temporary; force = true)
        rethrow()
    end
    return path
end

function load_snapshot(
        case::GauntletCase;
        path::AbstractString = snapshot_path(case)
)
    isfile(path) || throw(ArgumentError(
        "Gauntlet snapshot is missing: $path. Record it explicitly with " *
        "LINECABLEMODELS_GAUNTLET_MODE=record.",
    ))
    snapshot = JLD2.load(path)
    required = (
        "format_version", "case_name", "case_sha256", "line_parameters",
        "port_order", "kron_reduced", "pscad_version", "pscad_formulation",
        "pscad_elapsed_seconds", "julia_benchmark", "recorded_at_utc"
    )
    missing = filter(key -> !haskey(snapshot, key), required)
    isempty(missing) || throw(ArgumentError(
        "Gauntlet snapshot $path is missing fields: $(join(missing, ", "))",
    ))
    snapshot["format_version"] == SNAPSHOT_FORMAT_VERSION || throw(ArgumentError(
        "Gauntlet snapshot $path has unsupported format version $(snapshot["format_version"])",
    ))
    snapshot["case_name"] == string(case.name) || throw(ArgumentError(
        "Gauntlet snapshot $path belongs to case $(snapshot["case_name"]), not $(case.name)",
    ))
    snapshot["case_sha256"] == case_digest(case) || throw(ArgumentError(
        "Gauntlet snapshot $path does not match $(case.source_file). " *
        "Review the case and record the snapshot explicitly.",
    ))
    snapshot["pscad_formulation"] == string(case.pscad_formulation) ||
        throw(ArgumentError(
            "Gauntlet snapshot PSCAD formulation does not match the case definition",
        ))
    reference = snapshot["line_parameters"]
    reference isa LineParameters || throw(ArgumentError(
        "Gauntlet snapshot $path does not contain Engine.LineParameters",
    ))
    validate_structure(
        case,
        reference;
        port_order = Vector{String}(snapshot["port_order"]),
        kron_reduced = snapshot["kron_reduced"]
    )
    return (; reference, metadata = snapshot)
end

function run_snapshot(
        case::GauntletCase;
        path::AbstractString = snapshot_path(case)
)
    loaded = load_snapshot(case; path)
    candidate = compute!(case.problem, case.formulation)
    validate_structure(case, candidate)
    frequencies(loaded.reference) == frequencies(candidate) || throw(ArgumentError(
        "$(case.name) snapshot and candidate frequencies differ",
    ))
    comparison = compare(loaded.reference, candidate)
    return (
        mode = :snapshot,
        reference = loaded.reference,
        candidate,
        comparison,
        metadata = loaded.metadata,
        pscad = nothing,
        legacy = nothing,
        timings = nothing
    )
end
