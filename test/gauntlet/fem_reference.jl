using Dates
using Gmsh
using JLD2
using LinearAlgebra
using LineCableModels
using SHA

module FEMReferenceCaseLoader
using SHA
import TOML
import LineCableModels
using LineCableModels: AbstractGrid, Grid, Gridspace

const GAUNTLET_ROOT = @__DIR__
include(joinpath(GAUNTLET_ROOT, "case_loader.jl"))
end

using .FEMReferenceCaseLoader: case_index, load_case

const GETDP = get(
    ENV,
    "LINECABLEMODELS_GETDP",
    something(Sys.which("getdp"), "")
)
isempty(GETDP) && error(
    "LINECABLEMODELS_GETDP is unset and getdp is not on PATH"
)

const ROOT = joinpath(
    pkgdir(LineCableModels),
    ".linecablemodels",
    "fem",
    "gauntlet",
    "corrected_fullband"
)
const GETDP_ASSET_ROOT = joinpath(
    pkgdir(LineCableModels),
    "ext",
    "LineCableModelsGmshExt",
    "getdp"
)
const GETDP_ASSETS = (
    joinpath(GETDP_ASSET_ROOT, "model.pro"),
    joinpath(GETDP_ASSET_ROOT, "jacobian_integration.pro"),
    joinpath(GETDP_ASSET_ROOT, "quasi_tem.pro")
)
const MODEL_SOURCE_SHA256 = bytes2hex(sha256(join(read.(GETDP_ASSETS, String), '\0')))

relative_error(actual, reference) =
    norm(actual - reference) / max(norm(reference), eps(Float64))
symmetry_error(matrix) =
    norm(matrix - transpose(matrix)) / max(norm(matrix), eps(Float64))

function frequency_spectra(z, y, frequencies)
    count = length(frequencies)
    spectra = (
        resistance = Vector{Float64}(undef, count),
        inductance = Vector{Float64}(undef, count),
        conductance = Vector{Float64}(undef, count),
        capacitance = Vector{Float64}(undef, count)
    )
    tolerances = NamedTuple{keys(spectra)}(
        ntuple(_ -> Vector{Float64}(undef, count), length(spectra))
    )
    for index in eachindex(frequencies)
        omega = 2pi * frequencies[index]
        values = (
            resistance = real(@view z[:, :, index]),
            inductance = imag(@view z[:, :, index]) / omega,
            conductance = real(@view y[:, :, index]),
            capacitance = imag(@view y[:, :, index]) / omega
        )
        for quantity in keys(values)
            matrix = getproperty(values, quantity)
            getproperty(spectra, quantity)[index] = eigmin(Symmetric(matrix))
            getproperty(tolerances, quantity)[index] =
                1.0e-8 * max(opnorm(matrix), eps(Float64))
        end
    end
    return spectra, tolerances
end

function validate_result(case_id, model, result)
    frequencies = model.problem.frequencies
    terminal_count = length(model.port_order)
    expected = (terminal_count, terminal_count, length(frequencies))
    size(result.Z.values) == expected || error(
        "$case_id returned Z size $(size(result.Z.values)); expected $expected"
    )
    size(result.Y.values) == expected || error(
        "$case_id returned Y size $(size(result.Y.values)); expected $expected"
    )
    result.f == frequencies || error(
        "$case_id returned a different frequency vector"
    )

    primitive = result.details.fem.primitive
    size(primitive.Z_primitive) == expected || error(
        "$case_id returned primitive Z size $(size(primitive.Z_primitive)); " *
        "expected $expected"
    )
    size(primitive.P_primitive) == expected || error(
        "$case_id returned primitive P size $(size(primitive.P_primitive)); " *
        "expected $expected"
    )
    for (label, values) in (
            (:Z, result.Z.values),
            (:Y, result.Y.values),
            (:Z_primitive, primitive.Z_primitive),
            (:P_primitive, primitive.P_primitive)
        )
        all(isfinite, values) || error("$case_id returned non-finite $label")
    end

    z_symmetry = [
        symmetry_error(@view primitive.Z_primitive[:, :, index])
        for index in eachindex(frequencies)
    ]
    p_symmetry = [
        symmetry_error(@view primitive.P_primitive[:, :, index])
        for index in eachindex(frequencies)
    ]
    y_symmetry = [
        symmetry_error(@view result.Y.values[:, :, index])
        for index in eachindex(frequencies)
    ]
    observations = NamedTuple[]
    for (quantity, values, matrices, threshold) in (
            (:Z_primitive, z_symmetry, primitive.Z_primitive, 1.0e-10),
            (:P_primitive, p_symmetry, primitive.P_primitive, 1.0e-8),
            (:Y, y_symmetry, result.Y.values, 1.0e-10)
        )
        maximum_value, frequency_index = findmax(values)
        maximum_value <= threshold && continue
        matrix = @view matrices[:, :, frequency_index]
        difference = matrix - transpose(matrix)
        _, cartesian_pair = findmax(abs, difference)
        terminal_pair = Tuple(cartesian_pair)
        push!(observations, (;
            mode = :reciprocity,
            quantity,
            maximum = maximum_value,
            threshold,
            frequency_index,
            frequency_hz = frequencies[frequency_index],
            terminal_pair,
            pair_difference = difference[cartesian_pair]
        ))
    end

    inversion_residuals = result.details.fem.inversion_residuals
    condition_numbers = result.details.fem.condition_numbers
    length(inversion_residuals) == length(frequencies) || error(
        "$case_id returned the wrong inversion-residual count"
    )
    all(isfinite, condition_numbers) || error(
        "$case_id returned a non-finite P condition number"
    )
    maximum_residual, residual_index = findmax(inversion_residuals)
    maximum_residual <= 1.0e-10 || push!(observations, (;
        mode = :inversion_residual,
        quantity = :P_primitive,
        maximum = maximum_residual,
        threshold = 1.0e-10,
        frequency_index = residual_index,
        frequency_hz = frequencies[residual_index]
    ))

    spectra, tolerances = frequency_spectra(
        primitive.Z_primitive,
        result.Y.values,
        frequencies
    )
    for quantity in keys(spectra)
        values = getproperty(spectra, quantity)
        limits = getproperty(tolerances, quantity)
        failed = findall(values .< .-limits)
        isempty(failed) && continue
        severity = map(index -> -values[index] / limits[index], failed)
        _, worst_offset = findmax(severity)
        frequency_index = failed[worst_offset]
        push!(observations, (;
            mode = :non_passive,
            quantity,
            minimum_eigenvalue = values[frequency_index],
            tolerance = limits[frequency_index],
            frequency_index,
            frequency_hz = frequencies[frequency_index],
            failed_frequency_count = length(failed)
        ))
    end

    expected_invocations = length(frequencies) * length(result.details.fem.terminal_ids)
    actual_invocations = result.details.fem.run.getdp_invocations
    actual_invocations == expected_invocations || push!(observations, (;
        mode = :invocation_count,
        quantity = :getdp_jobs,
        actual = actual_invocations,
        expected = expected_invocations
    ))
    return (;
        z_symmetry,
        p_symmetry,
        y_symmetry,
        inversion_residuals = copy(inversion_residuals),
        condition_numbers = copy(condition_numbers),
        spectra,
        tolerances,
        observations
    )
end

function artifact_path(case_id)
    joinpath(ROOT, string(case_id), "reference.jld2")
end

function valid_existing(path, model)
    isfile(path) || return false
    try
        document = JLD2.load(path)
        return document["schema_version"] == 3 &&
               document["case_source_sha256"] == model.source_sha256 &&
               document["getdp_model_sha256"] == MODEL_SOURCE_SHA256 &&
               document["frequencies"] == model.problem.frequencies &&
               document["port_order"] == model.port_order
    catch
        return false
    end
end

function record_result(case_id, model, result, elapsed_seconds, checks)
    path = artifact_path(case_id)
    mkpath(dirname(path))
    temporary = path * ".new"
    primitive = result.details.fem.primitive
    JLD2.jldsave(
        temporary;
        schema_version = 3,
        kind = :fem_corrected_fullband_reference,
        status = isempty(checks.observations) ? :validated : :observed_anomalies,
        case_id = string(case_id),
        case_source_sha256 = model.source_sha256,
        getdp_model_sha256 = MODEL_SOURCE_SHA256,
        frequencies = copy(result.f),
        port_order = copy(model.port_order),
        terminal_ids = copy(result.details.fem.terminal_ids),
        Z = copy(result.Z.values),
        Y = copy(result.Y.values),
        Z_primitive = copy(primitive.Z_primitive),
        P_primitive = copy(primitive.P_primitive),
        phase_map = copy(primitive.phase_map),
        checks,
        observations = checks.observations,
        mesh_source = result.details.fem.run.mesh_source,
        mesh_fingerprint = result.details.fem.run.mesh_fingerprint,
        getdp_invocations = result.details.fem.run.getdp_invocations,
        elapsed_seconds,
        recorded_at_utc = string(now(UTC))
    )
    mv(temporary, path; force = true)
    open(path * ".sha256", "w") do io
        println(io, bytes2hex(sha256(read(path))), "  ", basename(path))
    end
    return path
end

function write_hash(path)
    open(path * ".sha256", "w") do io
        println(io, bytes2hex(sha256(read(path))), "  ", basename(path))
    end
    return path
end

function archive_stale_reference(path)
    digest = bytes2hex(sha256(read(path)))
    archived = joinpath(
        dirname(path),
        "reference_stale_$(first(digest, 12)).jld2"
    )
    isfile(archived) || cp(path, archived)
    write_hash(archived)
    return archived
end

function next_failure_attempt(directory)
    indices = Int[]
    for entry in readdir(directory)
        matched = match(r"^failure_attempt_(\d{3})\.jld2$", entry)
        matched === nothing || push!(indices, parse(Int, matched[1]))
    end
    return joinpath(
        directory,
        "failure_attempt_$(lpad(isempty(indices) ? 1 : maximum(indices) + 1, 3, '0')).jld2"
    )
end

function archive_current_failure(directory)
    current = joinpath(directory, "failure.jld2")
    isfile(current) || return nothing
    digest = sha256(read(current))
    for entry in readdir(directory; join = true)
        occursin(r"^failure_attempt_\d{3}\.jld2$", basename(entry)) || continue
        sha256(read(entry)) == digest && return entry
    end
    archived = next_failure_attempt(directory)
    cp(current, archived)
    write_hash(archived)
    return archived
end

function record_failure(case_id, model, result, exception, elapsed_seconds)
    directory = joinpath(ROOT, string(case_id))
    mkpath(directory)
    archive_current_failure(directory)
    path = next_failure_attempt(directory)
    temporary = path * ".new"
    run_directory = if result === nothing
        hasproperty(exception, :run_directory) ? exception.run_directory : nothing
    else
        result.details.fem.run.run_directory
    end
    JLD2.jldsave(
        temporary;
        schema_version = 3,
        kind = :fem_corrected_fullband_failure,
        status = :execution_or_invariant_failure,
        case_id = string(case_id),
        case_source_sha256 = model.source_sha256,
        getdp_model_sha256 = MODEL_SOURCE_SHA256,
        frequencies = copy(model.problem.frequencies),
        port_order = copy(model.port_order),
        exception_type = string(typeof(exception)),
        message = sprint(showerror, exception),
        run_directory,
        Z = result === nothing ? nothing : copy(result.Z.values),
        Y = result === nothing ? nothing : copy(result.Y.values),
        Z_primitive = result === nothing ? nothing :
                      copy(result.details.fem.primitive.Z_primitive),
        P_primitive = result === nothing ? nothing :
                      copy(result.details.fem.primitive.P_primitive),
        elapsed_seconds,
        recorded_at_utc = string(now(UTC))
    )
    mv(temporary, path; force = true)
    write_hash(path)
    latest = joinpath(directory, "failure.jld2")
    cp(path, latest; force = true)
    write_hash(latest)
    return path
end

function run_case(case_id)
    model = load_case(case_id)
    path = artifact_path(case_id)
    if valid_existing(path, model)
        println("REUSE\t", case_id, "\t", path)
        flush(stdout)
        return true
    end
    if isfile(path)
        archived = archive_stale_reference(path)
        println("STALE_ARCHIVED\t", case_id, "\t", archived)
        flush(stdout)
    end
    archive_current_failure(dirname(path))

    frequencies = model.problem.frequencies
    length(frequencies) == 101 || error(
        "$case_id has $(length(frequencies)) frequencies; expected 101"
    )
    first(frequencies) == 1.0 && last(frequencies) == 1.0e6 || error(
        "$case_id does not retain its original 1 Hz--1 MHz range"
    )
    println(
        "BEGIN\t", case_id,
        "\tfrequencies=", length(frequencies),
        "\tterminals=", length(model.port_order)
    )
    flush(stdout)
    formulation = Formulation(
        :LineCableModelsFEM;
        options = (
            reduce_bundle = false,
            kron_reduction = false,
            ideal_transposition = false,
            temperature_correction = true
        ),
        fem_options = (
            getdp_executable = GETDP,
            getdp_verbosity = 0,
            gmsh_verbosity = 0,
            keep_run_directory = true,
            mesh_policy = :reuse,
            plot_field_maps = false
        )
    )
    started = time()
    result = nothing
    try
        result = compute(model.problem, formulation; options = (trace = true,))
        elapsed_seconds = time() - started
        checks = validate_result(case_id, model, result)
        path = record_result(case_id, model, result, elapsed_seconds, checks)
        state = isempty(checks.observations) ? "PASS" : "OBSERVED"
        println(
            state, "\t", case_id,
            "\telapsed_s=", round(elapsed_seconds; digits = 3),
            "\tanomalies=", length(checks.observations),
            "\tmax_cond=", maximum(checks.condition_numbers),
            "\tmax_inv_res=", maximum(checks.inversion_residuals),
            "\trecord=", path
        )
        flush(stdout)
        return true
    catch exception
        elapsed_seconds = time() - started
        path = record_failure(
            case_id,
            model,
            result,
            exception,
            elapsed_seconds
        )
        println(
            stderr,
            "FAIL\t", case_id,
            "\telapsed_s=", round(elapsed_seconds; digits = 3),
            "\t", sprint(showerror, exception),
            "\trecord=", path
        )
        flush(stderr)
        return false
    end
end

available = sort!(collect(keys(case_index())); by = string)
requested = isempty(ARGS) ? available : Symbol.(ARGS)
all(case_id -> case_id in available, requested) || error(
    "unknown case requested; available cases are $(join(available, ", "))"
)

println(
    "PLAN\tcases=", length(requested),
    "\tfrequencies_per_case=101\trange_hz=1:1e6"
)
flush(stdout)
statuses = Bool[]
for case_id in requested
    push!(statuses, run_case(case_id))
end
passed = count(identity, statuses)
failed = length(statuses) - passed
println("COMPLETE\tpassed=", passed, "\tfailed=", failed)
failed == 0 || exit(1)
