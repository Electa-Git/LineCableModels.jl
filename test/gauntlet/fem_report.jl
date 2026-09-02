using Dates
using JLD2
using LineCableModels
using SHA

const ROOT = joinpath(
    pkgdir(LineCableModels),
    ".linecablemodels",
    "fem",
    "gauntlet",
    "corrected_fullband"
)

function first_crossing(values, predicate)
    index = findfirst(predicate, values)
    return index === nothing ? 0 : index
end

function violation_count(values, tolerances)
    count(index -> values[index] < -tolerances[index], eachindex(values))
end

function summary_row(case_id, path, document)
    checks = document["checks"]
    frequencies = document["frequencies"]
    p_first = first_crossing(checks.p_symmetry, >(1.0e-8))
    conductance = checks.spectra.conductance
    conductance_tolerance = checks.tolerances.conductance
    g_first = first_crossing(
        eachindex(conductance),
        index -> conductance[index] < -conductance_tolerance[index]
    )
    p_peak = argmax(checks.p_symmetry)
    g_worst = argmin(conductance)
    return (;
        case = case_id,
        status = string(document["status"]),
        terminals = length(document["port_order"]),
        observations = length(get(document, "observations", NamedTuple[])),
        P_first_reciprocity_index = p_first,
        P_first_reciprocity_hz = p_first == 0 ? NaN : frequencies[p_first],
        P_first_reciprocity_error = p_first == 0 ? NaN :
                                    checks.p_symmetry[p_first],
        P_peak_reciprocity_index = p_peak,
        P_peak_reciprocity_hz = frequencies[p_peak],
        P_peak_reciprocity_error = checks.p_symmetry[p_peak],
        G_first_nonpassive_index = g_first,
        G_first_nonpassive_hz = g_first == 0 ? NaN : frequencies[g_first],
        G_first_minimum_eigenvalue = g_first == 0 ? NaN : conductance[g_first],
        G_first_tolerance = g_first == 0 ? NaN : conductance_tolerance[g_first],
        G_nonpassive_frequency_count = violation_count(
            conductance, conductance_tolerance
        ),
        G_worst_index = g_worst,
        G_worst_hz = frequencies[g_worst],
        G_worst_minimum_eigenvalue = conductance[g_worst],
        G_worst_tolerance = conductance_tolerance[g_worst],
        Z_peak_reciprocity_error = maximum(checks.z_symmetry),
        Y_peak_reciprocity_error = maximum(checks.y_symmetry),
        maximum_condition_number = maximum(checks.condition_numbers),
        maximum_inversion_residual = maximum(checks.inversion_residuals),
        getdp_invocations = document["getdp_invocations"],
        elapsed_seconds = document["elapsed_seconds"],
        artifact = path,
        artifact_sha256 = bytes2hex(sha256(read(path)))
    )
end

function observation_rows(case_id, path, document)
    rows = NamedTuple[]
    for observation in get(document, "observations", NamedTuple[])
        push!(rows,
            (;
                case = case_id,
                status = string(document["status"]),
                mode = string(observation.mode),
                quantity = string(observation.quantity),
                frequency_index = get(observation, :frequency_index, 0),
                frequency_hz = get(observation, :frequency_hz, NaN),
                maximum = get(observation, :maximum, NaN),
                threshold = get(observation, :threshold, NaN),
                minimum_eigenvalue = get(
                    observation, :minimum_eigenvalue, NaN
                ),
                tolerance = get(observation, :tolerance, NaN),
                failed_frequency_count = get(
                    observation, :failed_frequency_count, 0
                ),
                terminal_pair = string(get(observation, :terminal_pair, "")),
                pair_difference = string(get(observation, :pair_difference, "")),
                artifact = path
            ))
    end
    return rows
end

function write_tsv(path, rows)
    isempty(rows) && error("cannot write empty report $path")
    temporary = path * ".new"
    names = keys(first(rows))
    open(temporary, "w") do io
        println(io, join(string.(names), '\t'))
        for row in rows
            values = replace.(
                string.((getproperty(row, name) for name in names)),
                '\t' => ' ',
                '\n' => ' '
            )
            println(io, join(values, '\t'))
        end
    end
    mv(temporary, path; force = true)
    return path
end

function write_hash(path)
    open(path * ".sha256", "w") do io
        println(io, bytes2hex(sha256(read(path))), "  ", basename(path))
    end
    return path
end

summaries = NamedTuple[]
observations = NamedTuple[]
for directory in sort!(filter(isdir, readdir(ROOT; join = true)))
    path = joinpath(directory, "reference.jld2")
    isfile(path) || continue
    document = JLD2.load(path)
    document["schema_version"] == 3 || error(
        "unsupported reference schema at $path"
    )
    case_id = string(document["case_id"])
    push!(summaries, summary_row(case_id, path, document))
    append!(observations, observation_rows(case_id, path, document))
end

isempty(summaries) && error("no FEM references were found in $ROOT")
summary_path = write_hash(write_tsv(
    joinpath(ROOT, "reference_summary.tsv"), summaries
))
observation_path = write_hash(write_tsv(
    joinpath(ROOT, "reference_observations.tsv"), observations
))
jld_path = joinpath(ROOT, "reference_summary.jld2")
temporary = jld_path * ".new"
JLD2.jldsave(
    temporary;
    schema_version = 1,
    kind = :fem_corrected_fullband_summary,
    shunt_processing = :reciprocity_enforced_before_potential_inversion,
    raw_reciprocity_quantity = :P_primitive,
    summaries,
    observations,
    recorded_at_utc = string(now(UTC))
)
mv(temporary, jld_path; force = true)
write_hash(jld_path)
println(
    "COMPLETE\treferences=", length(summaries),
    "\tobservations=", length(observations),
    "\tsummary=", summary_path,
    "\tobservation_report=", observation_path,
    "\tjld=", jld_path
)
