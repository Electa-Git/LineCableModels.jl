# Read DC behaviour archives and write CSV and LaTeX tables.
# Run: `julia --project=test benchmarks/dc_behaviour_tables.jl <run folders>... <output folder>`
using TOML
using LineCableModels.PSCAD: read_pscad_fit_report

length(ARGS) >= 2 || error("Expected run folders followed by an output folder")
folders = abspath.(ARGS[1:end-1])
output = mkpath(abspath(last(ARGS)))
tex_text(value) = replace(string(value), '\\' => raw"\textbackslash{}", '_' => raw"\_\allowbreak{}",
    '%' => raw"\%", '&' => raw"\&", '#' => raw"\#", '$' => raw"\$", '{' => raw"\{", '}' => raw"\}")
csv_text(value) = "\"" * replace(value === nothing ? "" : string(value), '"' => "\"\"") * "\""
tex_cell(::Nothing) = "{---}"
tex_cell(value::AbstractFloat) = isfinite(value) ? "\\num[round-mode=figures,round-precision=3]{$value}" : tex_text(value)
tex_cell(value) = tex_text(value)

# Convert native report data to this table schema.
function fit_report_columns(reports)
    fit_index = findfirst(report -> report.reports.fitting !== nothing, reports)
    scan_index = findfirst(report -> report.scan !== nothing, reports)
    enforcement_index = findfirst(report -> report.enforcement !== nothing, reports)
    fitting = fit_index === nothing ? nothing : reports[fit_index]
    absent_fit = (before_dc = nothing, after_dc = nothing, poles = nothing, poles_source = nothing)
    yc = something(fitting === nothing ? nothing : fitting.admittance, absent_fit)
    h = something(fitting === nothing ? nothing : fitting.propagation, absent_fit)
    absent_errors = (maximum_error_percent = nothing, rms_error_percent = nothing)
    yc_before, yc_after = something(yc.before_dc, absent_errors), something(yc.after_dc, absent_errors)
    h_before, h_after = something(h.before_dc, absent_errors), something(h.after_dc, absent_errors)
    raw = fitting === nothing ? (fitting = nothing, before_dc = nothing, after_dc = nothing) : fitting.reports
    scan = something(scan_index === nothing ? nothing : reports[scan_index].scan,
        (band = nothing, ranges = nothing, smallest_eigenvalue = nothing, outcome = nothing, report = nothing))
    enforcement = something(enforcement_index === nothing ? nothing : reports[enforcement_index].enforcement,
        (ranges = nothing, smallest_eigenvalue = nothing, outcome = nothing, report = nothing,
            initial_smallest_eigenvalue = nothing, final_smallest_eigenvalue = nothing, maximum_abs_error_percent = nothing))
    data = Dict{String, Any}(key => something(value, "not found") for (key, value) in (
        "Yc_before_dc_max_error_percent" => (raw.before_dc === nothing ? nothing : yc_before.maximum_error_percent),
        "Yc_after_dc_max_error_percent" => yc_after.maximum_error_percent,
        "Yc_reported_poles" => yc.poles,
        "H_before_dc_max_error_percent" => (raw.before_dc === nothing ? nothing : h_before.maximum_error_percent),
        "H_before_dc_rms_error_percent" => (raw.before_dc === nothing ? nothing : h_before.rms_error_percent),
        "H_after_dc_max_error_percent" => h_after.maximum_error_percent,
        "H_after_dc_rms_error_percent" => h_after.rms_error_percent,
        "H_final_poles" => (h.poles === nothing ? nothing : join(h.poles, ", ")),
        "maximum_residue_pole_ratio" => (fitting === nothing ? nothing : fitting.maximum_residue_pole_ratio),
        "before_dc_report" => raw.before_dc,
        "after_dc_report" => raw.after_dc,
        "scan_report" => scan.report,
        "scan_violation_ranges" => nothing,
        "scan_smallest_eigenvalue" => scan.smallest_eigenvalue,
        "scan_outcome" => scan.outcome,
        "enforcement_report" => enforcement.report,
        "enforcement_violation_ranges" => nothing,
        "enforcement_smallest_eigenvalue" => enforcement.smallest_eigenvalue,
        "enforcement_outcome" => enforcement.outcome,
        "enforcement_initial_smallest_eigenvalue" => enforcement.initial_smallest_eigenvalue,
        "enforcement_last_reported_smallest_eigenvalue" => enforcement.final_smallest_eigenvalue,
        "enforcement_max_abs_error_percent" => enforcement.maximum_abs_error_percent))
    for (name, ranges) in (("scan", scan.ranges), ("enforcement", enforcement.ranges))
        if ranges !== nothing && !isempty(ranges)
            data[name * "_violation_ranges"] = join(["$(row.lower_Hz) to $(row.upper_Hz) Hz: $(row.smallest_eigenvalue)" for row in ranges], "; ")
        end
    end
    return data
end

dc_correction_label(value) = value == 1 ? "Functional Form" : value == 2 ? "Add Pole/Residue" : "unknown ($value)"
passivity_label(value) = value == 0 ? "Disable" : value == 1 ? "Scan only" : value == 2 ? "Enforcement" : "unknown ($value)"

function write_tables(name, caption, columns, rows; tex_columns = columns, layout = repeat("l", length(tex_columns) - 1))
    open(joinpath(output, "$name.csv"), "w") do io
        println(io, join(csv_text.(columns), ','))
        for row in rows
            println(io, join([csv_text(get(row, key, "")) for key in columns], ','))
        end
    end
    columns = tex_columns
    open(joinpath(output, "$name.tex"), "w") do io
        println(io, raw"\begingroup\scriptsize")
        println(io, "\\begin{longtable}{@{}", layout, "@{}}")
        println(io, "\\caption{", caption, raw"}\\\\")
        println(io, raw"\toprule")
        headers = [tex_text(replace(key,
            "Yc_before_dc_max_error_percent" => "Yc max before (%)",
            "Yc_after_dc_max_error_percent" => "Yc max after (%)",
            "Yc_reported_poles" => "Yc poles", "H_final_poles" => "H poles",
            "H_before_dc_max_error_percent" => "H max before (%)",
            "H_after_dc_max_error_percent" => "H max after (%)",
            "H_before_dc_rms_error_percent" => "H RMS before (%)",
            "H_after_dc_rms_error_percent" => "H RMS after (%)",
            "maximum_residue_pole_ratio" => "max residue to pole ratio",
            "final_window_receiving_voltage_relative_difference" => "V error",
            "fault_current_mean_consistency_relative_error" => "relative residual",
            "fault_current_consistency_error_A" => "residual (A)",
            "measured_receiving_current_orientation" => "measured orientation",
            "measured_receiving_current_A" => "measured I (A)",
            "fault_current_A" => "fault I (A)", "receiving_voltage_V" => "receiving V",
            "fault_resistance_ohm" => "R (ohm)", "passivity_outcome" => "passivity outcome",
            "scaled_rms" => "RMS deviation", "violation_ranges" => "ranges (Hz)", "smallest_eigenvalue" => "smallest eigenvalue", "dc_solution_fault_current_relative_error" => "DC error",
            "slowest_time_constant_seconds" => "slowest time",
            "yc_added_pole_time_constant_seconds" => "Yc 1/a0",
            "h_added_pole_time_constant_seconds" => "H 1/a0",
            "dc_solution_relative_error" => "DC rel. error",
            "final_window_relative_difference" => "window difference",
            "target_dc_relative_error" => "target rel. error",
            "steady_state_fault_current_relative_error" => "steady error",
            "final_window_fault_current_relative_difference" => "window difference",
            "steady_state_sending_voltage_relative_error" => "send V steady error",
            "steady_state_receiving_voltage_relative_error" => "recv V steady error",
            "fault_current_tail_variation" => "current variation",
            "sending_voltage_tail_variation" => "send V variation",
            "receiving_voltage_tail_variation" => "recv V variation", "relative_error" => "rel. error", "relative_deviation" => "error", "fault_current" => "fault I", "sending_voltage" => "send V", "receiving_voltage" => "recv V", "_" => " ")) for key in columns[2:end]]
        println(io, join(headers, " & "), raw"\\\midrule\endfirsthead")
        println(io, join(headers, " & "), raw"\\\midrule\endhead")
        previous_cable = ""
        for row in rows
            if row["cable"] != previous_cable
                previous_cable = row["cable"]
                println(io, "\\multicolumn{", length(columns) - 1, "}{l}{", tex_text(previous_cable), raw"}\\\\")
            end
            println(io, join([tex_cell(get(row, key, nothing)) for key in columns[2:end]], " & "), raw"\\\\")
        end
        println(io, raw"\end{longtable}", '\n', raw"\endgroup")
    end
end

# Later final attempts replace earlier ones. Export-only records cannot replace results.
final_status(status) = status in ("completed", "nonfinite_response", "line_constants_error", "known_station_failure", "julia_error", "unsettled_reference")

function fault_current!(row, resistance)
    row["fault_resistance_ohm"] = resistance
    row["fault_current_definition"] = "receiving voltage divided by fault resistance; into fault"
    sign = get(row, "measured_current_to_fault_sign", startswith(row["model"], "dc_pscad") ? 1 : -1)
    row["measured_receiving_current_orientation"] = get(row, "measured_receiving_current_orientation",
        sign == 1 ? "out of line into fault" : "into line at receiving end")
    row["measured_current_to_fault_sign"] = sign
    if haskey(row, "receiving_voltage_V")
        archived = get(row, "fault_current_A", nothing)
        archived === nothing || (row["archived_fault_current_A"] = archived)
        # Earlier archives stored minus the receiving-current channel for every backend.
        measured = get(row, "measured_receiving_current_A", archived === nothing ? nothing : -archived)
        row["fault_current_A"] = row["receiving_voltage_V"] / resistance
        if measured !== nothing
            row["measured_receiving_current_A"] = measured
            discrepancy = sign * measured - row["fault_current_A"]
            row["fault_current_consistency_error_A"] = discrepancy
            scale = abs(row["fault_current_A"])
            row["fault_current_mean_consistency_relative_error"] = iszero(scale) ?
                (iszero(discrepancy) ? 0.0 : Inf) : abs(discrepancy) / scale
            row["fault_current_consistency_assessed"] = isfinite(discrepancy)
        end
        if haskey(row, "receiving_voltage_tail_variation")
            row["fault_current_tail_variation"] = row["receiving_voltage_tail_variation"]
            row["tail_variation"] = row["fault_current_tail_variation"]
        end
    end
    for prefix in ("dc_solution_", "target_dc_")
        key = prefix * "receiving_voltage_V"
        haskey(row, key) || continue
        row[prefix * "fault_current_A"] = row[key] / resistance
    end
    return row
end

physical, finals, matrix, exports, timings, logs, sources = [Dict{String, Any}[] for _ in 1:7]
summaries = [TOML.parsefile(joinpath(folder, "dc_behaviour.toml")) for folder in folders]
for cable in unique(summary["cable"] for summary in summaries)
    indices = findall(summary -> summary["cable"] == cable, summaries)
    cable_folders = folders[indices]
    push!(physical, merge(last(summaries[indices]), Dict("source" => last(cable_folders))))
    attempts = Dict{String, Vector{Dict{String, Any}}}()
    for i in indices
        folder, summary = folders[i], summaries[i]
        names = filter(name -> isdir(joinpath(folder, name)) && occursin(r"^dc_(?:nlt|ulm_\d+|pscad_\d+)$", name), readdir(folder))
        station_final_runs = 0
        for name in names
            path = joinpath(folder, name, "status.toml")
            status = isfile(path) ? TOML.parsefile(path) : Dict{String, Any}(
                "status" => "not_recorded", "message" => "Archive has no status for this model")
            row = merge(status, Dict("cable" => cable, "model" => name, "source" => folder,
                "status_path" => path, "fault_resistance_ohm" => get(status, "fault_resistance_ohm", get(summary, "fault_resistance_ohm", 0.01))))
            push!(get!(attempts, name, Dict{String, Any}[]), row)
            station_final_runs += startswith(name, "dc_pscad") && final_status(status["status"])
        end
        push!(timings, merge(summary, Dict("source" => folder, "station_final_runs" => station_final_runs, "run_folder" => basename(folder))))
    end
    selected = Dict{String, Dict{String, Any}}()
    for (name, candidates) in attempts
        finished = filter(row -> final_status(row["status"]), candidates)
        chosen = isempty(finished) ? first(candidates) : last(finished)
        selected[name] = copy(chosen)
        for row in candidates
            push!(sources, Dict("cable" => cable, "model" => name, "status" => row["status"],
                "source" => row["source"], "selected" => row === chosen))
        end
    end
    haskey(selected, "dc_nlt") || error("No NLT reference for $cable")
    reference = selected["dc_nlt"]
    fault_current!(reference, reference["fault_resistance_ohm"])
    reference_settled = get(reference, "settling_assessed", false) && get(reference, "settled", false) && reference["status"] == "completed"
    models = ["dc_nlt"; sort(filter(name -> startswith(name, "dc_ulm_"), collect(keys(selected))); by = name -> parse(Int, last(split(name, '_'))));
        sort(filter(name -> startswith(name, "dc_pscad_"), collect(keys(selected))); by = name -> parse(Int, last(split(name, '_'))))]
    for name in models
        row = selected[name]
        # Preserve the original NLT measurement when normalizing the selected reference.
        name == "dc_nlt" || fault_current!(row, row["fault_resistance_ohm"])
        row["final_status"] = final_status(row["status"])
        if !startswith(name, "dc_pscad") || get(row, "scenario", "") == "long"
            assessed = get(row, "settling_assessed", false)
            variation = get(row, "fault_current_tail_variation", Inf)
            settled = assessed && isfinite(variation) && variation < 0.001 && row["status"] == "completed"
            row["settled"] = settled
            row["settling"] = assessed ? (settled ? "settled" : "unsettled") : "not_assessed"
            row["nonfinite_response"] = row["status"] == "nonfinite_response"
            row["divergence"] = row["nonfinite_response"] ? "nonfinite" :
                !assessed ? "not assessed" : settled ? "not detected" : "unsettled (finite)"
            for (field, unit) in (("fault_current", "A"), ("sending_voltage", "V"), ("receiving_voltage", "V"))
                ref_key = "$(field)_$unit"
                haskey(reference, ref_key) || continue
                target = reference[ref_key]
                if haskey(row, ref_key)
                    value = row[ref_key]
                    error = value - target
                    relative = error / target
                    row["$(field)_error_$unit"] = error
                    row["$(field)_relative_error"] = relative
                    row["final_window_$(field)_relative_difference"] = relative
                    if settled && reference_settled
                        row["steady_state_$(field)_relative_error"] = relative
                    else
                        pop!(row, "steady_state_$(field)_relative_error", nothing)
                    end
                end
                for prefix in ("dc_solution_", "target_dc_")
                    haskey(row, prefix * ref_key) || continue
                    row[prefix * "$(field)_error_$unit"] = row[prefix * ref_key] - target
                    row[prefix * "$(field)_relative_error"] = (row[prefix * ref_key] - target) / target
                end
            end
        end
        if startswith(name, "dc_pscad")
            fit = get(row, "exported_fit", get(row, "fit", Dict{String, Any}()))
            merge!(row, fit)
            for i in 1:4
                row["LT$i"] = get(row, "LT", fill("not found", 4))[i]
            end
            row["DCCOR_meaning"] = haskey(fit, "DCCOR") ? dc_correction_label(fit["DCCOR"]) : "not found"
            row["CPASS_meaning"] = haskey(fit, "CPASS") ? passivity_label(fit["CPASS"]) : "not found"
            diagnostics = fit_report_columns([])
            paths = get(row, "line_constants_log_paths", String[])
            if isempty(paths)
                directory = joinpath(row["source"], name)
                paths = [joinpath(root, file) for (root, _, files) in walkdir(directory)
                    for file in files if endswith(lowercase(file), ".log")]
            end
            for log_path in paths
                data = fit_report_columns(isfile(log_path) ? [read_pscad_fit_report(log_path)] : [])
                for (key, value) in data
                    value == "not found" || (diagnostics[key] = value)
                end
                push!(logs, merge(row, data, Dict("log_path" => log_path)))
            end
            isempty(paths) && push!(logs, merge(row, diagnostics, Dict("log_path" => "not found")))
            merge!(row, diagnostics)
            kind = get(row, "CPASS", 0) == 2 ? "enforcement" : "scan"
            row["passivity_report"] = get(row, "CPASS_meaning", "not found")
            row["passivity_outcome"] = row["$(kind)_outcome"]
            row["passivity_violation_ranges"] = row["$(kind)_violation_ranges"]
            row["passivity_smallest_eigenvalue"] = row["$(kind)_smallest_eigenvalue"]
            push!(exports, row)
            get(row, "scenario", "") == "short" ? push!(matrix, row) : push!(finals, row)
        else
            push!(finals, row)
        end
    end
end

write_tables("dc_physical", "Physical conductance in S/m; durations and travel times in seconds. The L/R estimate uses the first three baseline frequencies.",
    ["cable", "G_dc_S_per_m", "G_analytic_S_per_m", "tau_seconds", "duration_seconds", "step_seconds"], physical)
write_tables("dc_final", "Current means in A over the last ten percent. Window differences compare those means with NLT. Steady errors require both the model and NLT to pass the current settling criterion. Unsettled models have no steady error. Relative values are fractions.",
    ["cable", "model", "fault_current_A", "final_window_fault_current_relative_difference", "steady_state_fault_current_relative_error", "settling"], finals)
write_tables("dc_voltages", "Core voltage means in V over the last ten percent. Steady errors are relative fractions and require both the model and NLT to pass the current settling criterion.",
    ["cable", "model", "sending_voltage_V", "steady_state_sending_voltage_relative_error", "receiving_voltage_V", "steady_state_receiving_voltage_relative_error"], finals)
write_tables("dc_settling", "Tail variation is the range over the last ten percent divided by the magnitude of its mean. Values are fractions. A model is settled when its current variation is below 0.001. Missing evidence remains unassessed.",
    ["cable", "model", "settling", "fault_current_tail_variation", "sending_voltage_tail_variation", "receiving_voltage_tail_variation"], finals)
pole_rows = Dict{String, Any}[]
for row in filter(row -> startswith(row["model"], "dc_ulm"), finals)
    base = Dict("cable" => row["cable"], "model" => row["model"])
    push!(pole_rows, merge(base, Dict("function" => "slowest", "time_constant_seconds" => get(row, "slowest_time_constant_seconds", nothing))))
    if haskey(row, "yc_time_constant_seconds")
        push!(pole_rows, merge(base, Dict("function" => "Yc", "time_constant_seconds" => row["yc_time_constant_seconds"])))
        for (index, value) in enumerate(row["h_time_constants_seconds"])
            push!(pole_rows, merge(base, Dict("function" => "H$index", "time_constant_seconds" => value)))
        end
    end
    for quantity in ("yc", "h")
        rate = "$(quantity)_added_pole_a0_per_second"
        if haskey(row, rate)
            push!(pole_rows, merge(base, Dict("function" => "added $quantity", "a0_per_second" => row[rate],
                "time_constant_seconds" => row["$(quantity)_added_pole_time_constant_seconds"])))
        end
    end
end
write_tables("dc_poles", "Pole time constants in seconds: reciprocal of the smallest absolute real part, for Yc and each H group. The slowest row is the maximum across functions. Added poles show their rate in inverse seconds and its reciprocal.",
    ["cable", "model", "function", "time_constant_seconds", "a0_per_second"], pole_rows)
ulm_rows = filter(row -> startswith(row["model"], "dc_ulm"), finals)
write_tables("dc_fits", "ULM fault-current errors relative to settled NLT, as fractions. DC errors use the fitted circuit solution. Window differences use the last ten percent. Time constants and reciprocals of added-pole rates are in seconds. DC errors apply to every fit regardless of its transient settling flag.",
    ["cable", "model", "method", "passivity", "dc_solution_fault_current_relative_error", "final_window_fault_current_relative_difference", "settling", "slowest_time_constant_seconds", "yc_added_pole_time_constant_seconds", "h_added_pole_time_constant_seconds"], ulm_rows; tex_columns = ["cable", "model", "method", "passivity", "final_window_fault_current_relative_difference", "dc_solution_fault_current_relative_error", "settled", "slowest_time_constant_seconds", "yc_added_pole_time_constant_seconds", "h_added_pole_time_constant_seconds"])
solution_rows, check_rows = Dict{String, Any}[], Dict{String, Any}[]
for row in ulm_rows
    for (field, unit) in (("fault_current", "A"), ("sending_voltage", "V"), ("receiving_voltage", "V"))
        base = Dict{String, Any}("cable" => row["cable"], "model" => row["model"], "quantity" => field, "unit" => unit)
        push!(solution_rows, merge(base, Dict(
            "dc_solution" => get(row, "dc_solution_$(field)_$unit", nothing),
            "dc_solution_error" => get(row, "dc_solution_$(field)_error_$unit", nothing),
            "dc_solution_relative_error" => get(row, "dc_solution_$(field)_relative_error", nothing),
            "final_window" => get(row, "$(field)_$unit", nothing),
            "final_window_relative_difference" => get(row, "final_window_$(field)_relative_difference", nothing),
            "settling" => row["settling"])))
        push!(check_rows, merge(base, Dict(
            "target_dc_solution" => get(row, "target_dc_$(field)_$unit", nothing),
            "target_dc_error" => get(row, "target_dc_$(field)_error_$unit", nothing),
            "target_dc_relative_error" => get(row, "target_dc_$(field)_relative_error", nothing))))
    end
end
write_tables("dc_solutions", "Fitted DC circuit solutions and final-window means, compared with settled NLT. Fault current is receiving voltage divided by fault resistance, directed into the fault. Voltage entries refer to core 1. Errors are signed. Relative errors and window differences are fractions. Units apply to solutions and absolute errors.",
    ["cable", "model", "quantity", "unit", "dc_solution", "dc_solution_error", "dc_solution_relative_error", "final_window", "final_window_relative_difference", "settling"], solution_rows)
write_tables("dc_target_check", "Physical-target DC circuit solutions compared with settled NLT. Corrected fits use the targets kept by the package. Uncorrected fits use the independent physical DC limit. Errors are signed and relative errors are fractions. Units apply to solutions and absolute errors.",
    ["cable", "model", "quantity", "unit", "target_dc_solution", "target_dc_error", "target_dc_relative_error"], check_rows)
write_tables("dc_limits", "ULM DC errors: maximum absolute matrix difference divided by the maximum target magnitude. Every fit uses independent targets from the same physical DC limit.",
    ["cable", "model", "method", "passivity", "H0_relative_deviation", "Yc0_relative_deviation"], filter(row -> startswith(row["model"], "dc_ulm"), finals))
write_tables("dc_correction", "Short-window receiving voltage RMS deviations divided by the NLT peak. Voltage peaks are in kV. Orders are 20/20.",
    ["cable", "model", "status", "scaled_rms", "peak_receiving_voltage_kV", "wall_seconds"], matrix)
write_tables("dc_status", "Status and root messages for all models. Native messages and line-constants log paths also remain in each status archive.",
    ["cable", "model", "status", "settling", "message"], vcat(finals, matrix); layout = "lllp{0.5\\linewidth}")
export_fields = ["line_formulation", "LT1", "LT2", "LT3", "LT4", "enablf", "shntcab", "DCenab", "DCCOR", "DCCOR_meaning", "CPASS", "CPASS_meaning", "FS", "YMaxP", "AMaxP", "export_check"]
open(joinpath(output, "dc_export_variants.csv"), "w") do io
    columns = ["cable", "model", "scenario", export_fields...]
    println(io, join(csv_text.(columns), ','))
    for row in exports
        println(io, join([csv_text(get(row, key, "not found")) for key in columns], ','))
    end
end
export_values = Dict{String, Any}[]
for cable in unique([row["cable"] for row in exports]), field in export_fields
    rows = filter(row -> row["cable"] == cable, exports)
    values = sort(unique(string.([field == "DCCOR" && haskey(row, field) ? "$(row[field]) ($(dc_correction_label(row[field])))" :
        field == "CPASS" && haskey(row, field) ? "$(row[field]) ($(passivity_label(row[field])))" : get(row, field, "not found") for row in rows])))
    push!(export_values, Dict("cable" => cable, "field" => field, "distinct_values" => join(values, ", ")))
end
write_tables("dc_exports", "Distinct XML values from the package exporter. Conductance is in S/m and FS is in Hz. Loss tangents describe the exported dielectric at 50 Hz.",
    ["cable", "field", "distinct_values"], export_values; layout = "lp{0.75\\linewidth}")
write_tables("dc_pscad_fits", "Tline fitting errors in percent. Before correction uses the last Hphase attempt preceding Applying DC Correction. After correction uses Final Fit, which reports H. Missing headings give not found. The log gives no Yc error after correction.",
    ["cable", "model", "Yc_before_dc_max_error_percent", "Yc_after_dc_max_error_percent", "H_before_dc_max_error_percent", "H_after_dc_max_error_percent", "H_before_dc_rms_error_percent", "H_after_dc_rms_error_percent"], logs)
write_tables("dc_pscad_poles", "Tline pole counts and final maximum ratio of residue to pole. H counts follow delay-group order in Final Fit or the last Hmode attempt. Ratios use the final reported value, including uncorrected Hphase attempts. Yc counts come from its own fitting heading.",
    ["cable", "model", "Yc_reported_poles", "H_final_poles", "maximum_residue_pole_ratio"], logs)
passivity_rows = Dict{String, Any}[]
for row in logs
    kind = row["CPASS"] == 2 ? "enforcement" : "scan"
    push!(passivity_rows, merge(row, Dict("report" => passivity_label(row["CPASS"]),
        "violation_ranges" => row["$(kind)_violation_ranges"],
        "smallest_eigenvalue" => row["$(kind)_smallest_eigenvalue"], "outcome" => row["$(kind)_outcome"])))
end
write_tables("dc_pscad_passivity", "Native passivity reports for the selected CPASS mode. Scan ranges are in Hz. Enforcement eigenvalues come from reported point samples across iterations. Initial and last enforcement minima refer to those reports. Missing headings give not found. Raw reports preserve both scan and enforcement sections when present.",
    ["cable", "model", "report", "violation_ranges", "smallest_eigenvalue", "outcome", "enforcement_initial_smallest_eigenvalue", "enforcement_last_reported_smallest_eigenvalue", "enforcement_max_abs_error_percent"], passivity_rows;
    tex_columns = ["cable", "model", "report", "violation_ranges", "smallest_eigenvalue", "outcome"], layout = "llp{0.24\\linewidth}rp{0.27\\linewidth}")
open(joinpath(output, "dc_pscad_log_reports.toml"), "w") do io
    TOML.print(io, Dict("logs" => [Dict(key => row[key] for key in ("cable", "model", "log_path", "before_dc_report", "after_dc_report", "scan_report", "enforcement_report")) for row in logs]); sorted = true)
end
write_tables("dc_timings", "Wall times in seconds. Station simulated time sums the selected short and long windows; station wall time depends on native fitting and execution.",
    ["cable", "run_folder", "export_wall_seconds", "reference_wall_seconds", "wall_seconds", "station_runs", "station_final_runs", "station_simulated_seconds"], timings)
write(joinpath(output, "dc_sources.txt"), join(folders, '\n') * "\n")

pscad_long = filter(row -> startswith(row["model"], "dc_pscad"), finals)
write_tables("dc_pscad_long", "Consolidated PSCAD fault responses. Current is receiving voltage divided by fault resistance, into the fault. Errors and tails are fractions. Final-window errors apply to finite unsettled responses too. A settled flag requires tail variation below 0.001. Status and nonfinite diagnostics distinguish interrupted and diverged runs.",
    ["cable", "model", "status", "DCenab", "DCCOR_meaning", "CPASS_meaning", "FS", "shntcab", "fault_resistance_ohm", "fault_current_A", "receiving_voltage_V", "final_window_fault_current_relative_difference", "final_window_receiving_voltage_relative_difference", "settled", "fault_current_tail_variation", "divergence", "nonfinite_channel", "first_nonfinite_time_seconds"], pscad_long; tex_columns = ["cable", "model", "status", "final_window_receiving_voltage_relative_difference", "settled", "fault_current_tail_variation", "divergence"])
write_tables("dc_pscad_short", "Consolidated PSCAD short-window responses. RMS deviation is relative to the NLT voltage peak. Tline errors are in percent. Before and after refer to DC correction. Native passivity outcomes retain their log wording. Enforcement eigenvalues are the smallest reported across iterations.",
    ["cable", "model", "status", "DCCOR_meaning", "CPASS_meaning", "FS", "shntcab", "scaled_rms", "Yc_before_dc_max_error_percent", "Yc_after_dc_max_error_percent", "H_before_dc_max_error_percent", "H_after_dc_max_error_percent", "H_before_dc_rms_error_percent", "H_after_dc_rms_error_percent", "passivity_report", "passivity_outcome", "passivity_violation_ranges", "passivity_smallest_eigenvalue", "enforcement_initial_smallest_eigenvalue", "enforcement_last_reported_smallest_eigenvalue", "enforcement_max_abs_error_percent"], matrix; tex_columns = ["cable", "model", "status", "scaled_rms", "H_before_dc_max_error_percent", "H_after_dc_max_error_percent", "passivity_outcome"], layout = "lllllp{0.35\\linewidth}")
write_tables("dc_ulm", "Consolidated ULM fault responses. Current errors use receiving voltage divided by fault resistance. Error and tail columns are fractions. DC solution errors apply independently of transient settling. Pole time constants and reciprocals of added-pole rates are in seconds.",
    ["cable", "model", "status", "method", "passivity", "fault_current_A", "receiving_voltage_V", "dc_solution_fault_current_A", "dc_solution_receiving_voltage_V", "final_window_fault_current_relative_difference", "final_window_receiving_voltage_relative_difference", "dc_solution_fault_current_relative_error", "dc_solution_receiving_voltage_relative_error", "settled", "fault_current_tail_variation", "slowest_time_constant_seconds", "yc_added_pole_time_constant_seconds", "h_added_pole_time_constant_seconds"], ulm_rows; tex_columns = ["cable", "model", "method", "passivity", "final_window_fault_current_relative_difference", "dc_solution_fault_current_relative_error", "settled", "slowest_time_constant_seconds", "yc_added_pole_time_constant_seconds", "h_added_pole_time_constant_seconds"])
write_tables("dc_nlt", "Selected NLT references. Fault current is receiving voltage divided by fault resistance, into the fault. Final values use the last ten percent of each window. Tails are fractions. R is in ohms, voltages in V and current in A.",
    ["cable", "fault_current_A", "receiving_voltage_V", "sending_voltage_V", "fault_resistance_ohm", "settled", "fault_current_tail_variation", "source"], filter(row -> row["model"] == "dc_nlt", finals); tex_columns = ["cable", "fault_current_A", "receiving_voltage_V", "fault_resistance_ohm", "settled", "fault_current_tail_variation"])
write_tables("dc_current_consistency", "Consistency of final-window means. The native receiving-current orientation and its sign toward the fault are explicit. The signed residual is oriented measured current minus receiving voltage divided by fault resistance. Relative residuals divide its absolute value by the absolute derived current. Original current values remain available for comparison.",
    ["cable", "model", "fault_current_A", "measured_receiving_current_A", "measured_receiving_current_orientation", "measured_current_to_fault_sign", "fault_current_consistency_error_A", "fault_current_mean_consistency_relative_error", "fault_current_consistency_assessed", "archived_fault_current_A"], finals; tex_columns = ["cable", "model", "fault_current_A", "measured_receiving_current_A", "measured_receiving_current_orientation", "fault_current_consistency_error_A", "fault_current_mean_consistency_relative_error"], layout = "lllp{0.24\\linewidth}ll")
write_tables("dc_model_sources", "Archive selection for each cable and model. Later final attempts replace earlier final attempts in argument order. A non-final attempt cannot replace a final result. Without a final result, the first supplied attempt provides the metadata. Selected rows supply the consolidated tables. All paths identify read-only inputs.",
    ["cable", "model", "status", "selected", "source"], sources; layout = "lllp{0.55\\linewidth}")
write_tables("dc_incomplete_pscad", "PSCAD matrix entries without a final outcome in the supplied archives. Exported and started states are incomplete.",
    ["cable", "model", "status", "source"], filter(row -> !row["final_status"], exports); layout = "llp{0.65\\linewidth}")

println("Wrote consolidated DC behaviour tables to ", output)
