# Derive option comparisons from archived fits, responses and native logs.
using LineCableModels
using TOML
using LineCableModels.PSCAD: read_pscad_fit_report
using Statistics: mean

csv_cell(x) = "\"" * replace(x === nothing ? "" : string(x), '"' => "\"\"") * "\""
tex_text(x) = replace(string(x), '\\' => raw"\textbackslash{}", '_' => raw"\_\allowbreak{}",
    '%' => raw"\%", '&' => raw"\&", '#' => raw"\#", '$' => raw"\$", '{' => raw"\{", '}' => raw"\}")
tex_cell(::Nothing) = "{---}"
tex_cell(x::AbstractFloat) = isfinite(x) ? "\\num[round-mode=figures,round-precision=3]{$x}" : "{" * tex_text(x) * "}"
tex_cell(x) = tex_text(x)

function write_tables(destination, name, rows, columns; caption, tex_columns = columns)
    open(joinpath(destination, name * ".csv"), "w") do io
        println(io, join(csv_cell.(columns), ','))
        for row in rows
            println(io, join([csv_cell(get(row, key, nothing)) for key in columns], ','))
        end
    end
    open(joinpath(destination, name * ".tex"), "w") do io
        println(io, raw"\begingroup\tiny")
        println(io, "\\begin{longtable}{", repeat("l", length(tex_columns)), "}")
        println(io, "\\caption{", tex_text(caption), raw"}\\\toprule")
        println(io, join(tex_text.(tex_columns), " & "), raw"\\\midrule\endhead")
        for row in rows
            println(io, join([tex_cell(get(row, key, nothing)) for key in tex_columns], " & "), "\\\\")
        end
        println(io, raw"\bottomrule\end{longtable}\endgroup")
    end
end

function response_values(path)
    observed = import_data(:observed, path)
    q = only(filter(q -> q.request === voltage, observed.quantities))
    i = only(filter(q -> q.request === current, observed.quantities))
    return (; values = q.values, currents = i.values,
        times = q.coordinates.times, labels = q.coordinates.column_labels)
end

# Passivity controls do not change the unconstrained fit. Match all other native fields and the time grid.
function passivity_reference(row, runs)
    controls = ("CPASS", "ET_PE", "MER_PE", "MIT_PE", "FSP", "FEP", "NFP", "FDIS")
    matches = [run for run in runs if run["backend"] == "pscad" && run["status"] == "completed" &&
        run["scenario"] == row["scenario"] && get(run, "model", "") == row["model"] &&
        get(run, "step_seconds", nothing) == row["step_seconds"] &&
        get(run, "duration_seconds", nothing) == row["duration_seconds"] &&
        haskey(run, "exported_fit") && get(run["exported_fit"], "CPASS", 0) == 1 &&
        all(get(run["exported_fit"], field, nothing) == value
            for (field, value) in row["exported_fit"] if !(field in controls))]
    isempty(matches) && return nothing
    return matches[argmin([run["index"] for run in matches])]
end

function transient_metrics(response, reference; tmax = last(reference.times))
    response.labels == reference.labels || error("Terminal order differs")
    length(response.times) == length(reference.times) &&
        all(isapprox.(response.times, reference.times; rtol = 2e-6, atol = 1e-12)) || error("Reference time grid differs")
    core = only(findall(==("cable:1:core"), response.labels))
    a, b = response.values[2, core, :], reference.values[2, core, :]
    scale = maximum(abs, b)
    p = maximum(abs, a)
    arrival = findfirst(>=(0.1p), abs.(a))
    metrics = Dict{String, Any}("receiving_core_rms_deviation_kV" => sqrt(mean(abs2, a - b)) / 1000,
        "receiving_core_peak_deviation_kV" => maximum(abs, a - b) / 1000,
        "receiving_core_scaled_rms" => sqrt(mean(abs2, a - b)) / scale,
        "receiving_core_scaled_peak" => maximum(abs, a - b) / scale,
        "receiving_core_peak_kV" => p / 1000, "reference_core_peak_kV" => scale / 1000,
        "receiving_core_peak_error_kV" => (p - scale) / 1000,
        "receiving_core_peak_relative_error" => (p - scale) / scale,
        "arrival_time_us" => arrival === nothing || iszero(p) ? NaN : response.times[arrival] * 1e6)
    growth = Dict{String, Any}[]
    # Use the requested duration, as in the time-domain study.
    for side in 1:2, (terminal, label) in enumerate(response.labels)
        peaks = map((response, reference)) do run
            previous = findall(t -> 0.8tmax <= t < 0.9tmax, run.times)
            final = findall(t -> 0.9tmax <= t <= tmax, run.times)
            a = maximum(abs, run.values[side, terminal, previous])
            b = maximum(abs, run.values[side, terminal, final])
            ratio = iszero(a) ? (iszero(b) ? 1.0 : Inf) : b / a
            return (; ratio, late = b)
        end
        ratio_flag = peaks[1].ratio > 1 && peaks[1].ratio > 1.1peaks[2].ratio
        excess = (peaks[1].late - peaks[2].late) / scale
        push!(growth, Dict("end" => side == 1 ? "sending" : "receiving", "terminal" => label,
            "growth_ratio" => peaks[1].ratio, "reference_growth_ratio" => peaks[2].ratio,
            "late_peak_relative_to_scale" => peaks[1].late / scale,
            "reference_late_peak_relative_to_scale" => peaks[2].late / scale,
            "late_peak_excess_relative_to_scale" => excess,
            "growth_ratio_flag" => ratio_flag, "growth_flag" => ratio_flag && excess > 0.01))
        if side == 2 && endswith(label, ":sheath")
            peak = maximum(abs, response.values[side, terminal, :]) / 1000
            target = maximum(abs, reference.values[side, terminal, :]) / 1000
            metrics["$(label)_peak_kV"] = peak
            metrics["$(label)_peak_error_kV"] = peak - target
            metrics["$(label)_peak_relative_error"] = iszero(target) ? nothing : (peak - target) / target
        end
    end
    metrics["growth_ratio_flag"] = any(row["growth_ratio_flag"] for row in growth)
    metrics["growth_flag"] = any(row["growth_flag"] for row in growth)
    metrics["maximum_late_peak_relative_to_scale"] = maximum(row["late_peak_relative_to_scale"] for row in growth)
    return metrics, growth
end

# Pulse cutoff contains 99 percent of the source's continuous spectral energy.
# The ideal step has unbounded source bandwidth. Its resolved limit is the NLT band.
function spectral_limit(scenario)
    nlt_limit = (2scenario["N"] - 1) / (2scenario["Tfactor"] * scenario["duration_seconds"])
    scenario["name"] == "impulse" || return nlt_limit, Inf, "NLT upper frequency; ideal step has unbounded bandwidth"
    alpha, beta = 1.47e4, 2.47e6
    lo, hi = 0.0, 100beta
    for _ in 1:80
        w = (lo + hi) / 2
        fraction = (atan(w / alpha) / alpha - atan(w / beta) / beta) / (pi / 2 * (1 / alpha - 1 / beta))
        fraction < 0.99 ? (lo = w) : (hi = w)
    end
    cutoff = hi / (2pi)
    return min(nlt_limit, cutoff), cutoff, "99 percent source spectral energy, capped at NLT upper frequency"
end

# Convert native report data to this table schema.
function fit_report_columns(reports; fitting_band)
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
        "Yc_max_error_percent" => yc_before.maximum_error_percent,
        "Yc_after_dc_max_error_percent" => yc_after.maximum_error_percent,
        "Yc_reported_poles" => yc.poles,
        "H_max_error_percent" => h_before.maximum_error_percent,
        "H_rms_error_percent" => h_before.rms_error_percent,
        "H_after_dc_max_error_percent" => h_after.maximum_error_percent,
        "H_after_dc_rms_error_percent" => h_after.rms_error_percent,
        "H_final_poles" => (h.poles === nothing ? nothing : join(h.poles, ", ")),
        "H_poles_source" => (h.poles_source === nothing ? nothing : h.poles_source === :final_fit ? "Final Fit" : "last Hmode attempt"),
        "maximum_residue_pole_ratio" => (fitting === nothing ? nothing : fitting.maximum_residue_pole_ratio),
        "fitting_report" => raw.fitting,
        "before_dc_report" => raw.before_dc,
        "after_dc_report" => raw.after_dc,
        "scan_report" => scan.report,
        "scan_violation_ranges" => nothing,
        "scan_smallest_eigenvalue" => scan.smallest_eigenvalue,
        "scan_outcome" => scan.outcome,
        "scan_lower_Hz" => (scan.band === nothing ? nothing : scan.band.lower_Hz),
        "scan_upper_Hz" => (scan.band === nothing ? nothing : scan.band.upper_Hz),
        "enforcement_report" => enforcement.report,
        "enforcement_violation_ranges" => nothing,
        "enforcement_smallest_eigenvalue" => enforcement.smallest_eigenvalue,
        "enforcement_outcome" => enforcement.outcome,
        "enforcement_initial_smallest_eigenvalue" => enforcement.initial_smallest_eigenvalue,
        "enforcement_final_smallest_eigenvalue" => enforcement.final_smallest_eigenvalue,
        "enforcement_last_reported_smallest_eigenvalue" => enforcement.final_smallest_eigenvalue,
        "enforcement_max_abs_error_percent" => enforcement.maximum_abs_error_percent))
    for (name, ranges) in (("scan", scan.ranges), ("enforcement", enforcement.ranges))
        if ranges !== nothing && !isempty(ranges)
            data[name * "_violation_ranges"] = join(["$(row.lower_Hz) to $(row.upper_Hz) Hz: $(row.smallest_eigenvalue)" for row in ranges], "; ")
        end
    end
    data["scan_ranges"] = Dict{String, Any}[]
    if scan.ranges !== nothing
        for row in scan.ranges
            regions = String[]
            row.lower_Hz < fitting_band[1] && push!(regions, "below")
            row.upper_Hz >= fitting_band[1] && row.lower_Hz <= fitting_band[2] && push!(regions, "inside")
            row.upper_Hz > fitting_band[2] && push!(regions, "above")
            push!(data["scan_ranges"], Dict("lower_Hz" => row.lower_Hz, "upper_Hz" => row.upper_Hz,
                "smallest_eigenvalue" => row.smallest_eigenvalue, "relative_to_fitting_band" => join(regions, ", ")))
        end
    end
    return data
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) >= 2 || error("Usage: julia --project=test benchmarks/pscad_options_tables.jl <run folders>... <output folder>")
    folders = abspath.(ARGS[1:end-1])
    output = mkpath(abspath(last(ARGS)))
    summaries = [TOML.parsefile(joinpath(folder, "pscad_options.toml")) for folder in folders]
    final_status(status) = status in ("completed", "line_constants_error", "nonfinite_response", "execution_error", "julia_error")
    for cable in unique(summary["cable"] for summary in summaries)
        selected = Dict{String, Dict{String, Any}}()
        attempts = Dict{String, Any}[]
        cable_indices = findall(s -> s["cable"] == cable, summaries)
        destination = mkpath(joinpath(output, cable))
        for i in cable_indices
            folder = folders[i]
            for name in readdir(folder)
                if startswith(name, "options_fit_")
                    path = joinpath(folder, name, "status.toml")
                    isfile(path) || continue
                    fit_status = TOML.parsefile(path)
                    for scenario in summaries[i]["scenarios"]
                        run = "options_julia_$(scenario["name"])_$(fit_status["index"])_0_single"
                        isfile(joinpath(folder, run, "status.toml")) && continue
                        row = merge(fit_status, Dict("run" => run, "source" => folder, "cable" => cable,
                            "backend" => fit_status["method"], "scenario" => scenario["name"]))
                        row["status"] == "completed" && (row["status"] = "not_recorded")
                        push!(attempts, row)
                        if !haskey(selected, run) || final_status(row["status"])
                            selected[run] = row
                        end
                    end
                    continue
                end
                startswith(name, "options_pscad_") || startswith(name, "options_julia_") || continue
                path = joinpath(folder, name, "status.toml")
                identity = match(r"^options_(pscad|julia)_(step_single_point|impulse)_(\d+)", name)
                identity === nothing && continue
                row = Dict{String, Any}("status" => "not_recorded", "backend" => identity[1],
                    "scenario" => identity[2], "study" => "not_recorded", "index" => parse(Int, identity[3]))
                isfile(path) && merge!(row, TOML.parsefile(path))
                merge!(row, Dict("run" => name, "source" => folder, "cable" => cable))
                push!(attempts, row)
                if !haskey(selected, name) || final_status(row["status"])
                    selected[name] = row
                end
            end
        end
        rows, growth_rows, pole_rows, log_rows, range_rows = [Dict{String, Any}[] for _ in 1:5]
        log_cache = Dict{String, Dict{String, Any}}()
        for (name, row) in selected
            haskey(row, "exported_fit") || continue
            paths = get(row, "line_constants_log_paths", String[])
            # The copied archive remains usable after the native folder moves.
            local_paths = [joinpath(root, file) for (root, _, files) in walkdir(joinpath(row["source"], name))
                for file in files if endswith(lowercase(file), ".log")]
            paths = unique(vcat(filter(isfile, paths), local_paths))
            fit = row["exported_fit"]
            diagnostics = fit_report_columns(read_pscad_fit_report.(paths);
                fitting_band = (fit["FS"], fit["FE"]))
            log_cache[name] = merge(Dict("run" => name, "source" => row["source"], "paths" => paths), diagnostics)
        end
        reference_cache = Dict{String, Any}()
        fit_cache = Dict{String, Any}()
        baselines = Dict{Tuple{String, String}, Any}()
        for scenario in ("step_single_point", "impulse"), backend in ("ulm", "pscad")
            name = backend == "pscad" ? "options_pscad_$(scenario)_1" : "options_julia_$(scenario)_1_0_single"
            if haskey(selected, name) && selected[name]["status"] == "completed"
                baselines[(backend, scenario)] = response_values(joinpath(selected[name]["source"], name, "response.json"))
            end
        end
        for name in sort(collect(keys(selected)))
            row = copy(selected[name])
            folder = row["source"]
            if haskey(row, "exported_fit")
                merge!(row, row["exported_fit"])
                merge!(row, log_cache[name])
                row["poles_beyond_content"] = "not found: native log reports counts without pole locations"
            end
            scenario = only(s for s in summaries[last(cable_indices)]["scenarios"] if s["name"] == row["scenario"])
            cutoff, source_cutoff, definition = spectral_limit(scenario)
            row["spectral_limit_Hz"] = cutoff
            row["source_99_percent_energy_Hz"] = source_cutoff
            row["spectral_limit_definition"] = definition
            if haskey(row, "fit_archive")
                path = joinpath(folder, row["fit_archive"])
                fit = get!(fit_cache, path) do
                    import_data(:observed, path).gridpoint.options_fit
                end
                merge!(row, Dict(string(k) => v for (k, v) in pairs(fit.variant)))
                # Keep the time-step study identity when a simulation reuses its fitted baseline.
                row["study"] = row["divisor"] == 0 ? string(fit.variant.study) : "time_step"
                row["method"] = string(fit.variant.method)
                row["samples"] = length(fit.frequencies)
                row["transformation_frequency_Hz"] = fit.transformation_frequency
                row["passive_before"] = fit.passivity_before.passive
                row["passive_after"] = fit.passivity_after.passive
                row["passivity_worst_before"] = fit.passivity_before.worst
                row["passivity_worst_after"] = fit.passivity_after.worst
                row["passivity_before_report"] = string(fit.passivity_before)
                row["passivity_after_report"] = string(fit.passivity_after)
                for (key, value) in pairs(fit.fitting_error)
                    row["$(key)_rms_error"] = maximum(value.rms)
                    row["$(key)_max_error"] = maximum(value.max)
                end
                row["fitted_poles"] = sum(length, fit.poles)
                row["poles_beyond_content"] = sum(count(p -> abs(p) / (2pi) > cutoff, group) for group in fit.poles)
                propagation_ratios = Float64[]
                section_ratios = Float64[]
                for (group, (ps, rs)) in enumerate(zip(fit.poles, fit.residues)), (k, pole) in enumerate(ps)
                    ratio = fit.variant.method === :marti ? abs(rs[k]) / abs(pole) : maximum(abs, rs[:, :, k]) / abs(pole)
                    is_propagation = fit.variant.method === :marti ? group > length(fit.poles) ÷ 2 : group > 1
                    quantity = is_propagation ? "propagation" : fit.variant.method === :marti ? "impedance" : "admittance"
                    push!(is_propagation ? propagation_ratios : section_ratios, ratio)
                    push!(pole_rows, Dict("run" => name, "source" => folder, "scenario" => row["scenario"],
                        "quantity" => quantity, "group" => group, "pole_index" => k, "pole_real_per_s" => real(pole), "pole_imag_per_s" => imag(pole),
                        "pole_frequency_Hz" => abs(pole) / (2pi), "spectral_limit_Hz" => cutoff,
                        "beyond_content" => abs(pole) / (2pi) > cutoff, "residue_pole_ratio" => ratio))
                end
                row["maximum_residue_pole_ratio"] = maximum(propagation_ratios; init = 0.0)
                section_key = fit.variant.method === :marti ? "maximum_impedance_residue_pole_ratio_ohm" : "maximum_admittance_residue_pole_ratio_S"
                row[section_key] = maximum(section_ratios; init = 0.0)
                if fit.variant.enforce
                    before = get!(fit_cache, joinpath(folder, "options_fit_1/fit.json")) do
                        import_data(:observed, joinpath(folder, "options_fit_1/fit.json")).gridpoint.options_fit
                    end
                    for key in (:admittance, :propagation)
                        row["$(key)_max_error_change"] = maximum(getproperty(fit.fitting_error, key).max) -
                            maximum(getproperty(before.fitting_error, key).max)
                        row["$(key)_rms_error_change"] = maximum(getproperty(fit.fitting_error, key).rms) -
                            maximum(getproperty(before.fitting_error, key).rms)
                    end
                end
                if fit.marti !== nothing
                    row["transformation_imaginary_norm"] = maximum(fit.marti.transformation_imaginary_norm)
                    row["offdiagonal_ratio"] = string(fit.marti.offdiagonal_ratio)
                    row["phase_error_degrees"] = maximum(fit.marti.phase_error_degrees)
                end
            end
            comparison = nothing
            if row["backend"] == "pscad" && get(row, "CPASS", 0) == 2
                comparison = passivity_reference(row, values(selected))
                row["enforcement_transient_identical"] = "not recorded"
                row["enforcement_simulation"] = "not established"
                if comparison !== nothing
                    row["enforcement_reference_run"] = comparison["run"]
                    for (key, value) in log_cache[comparison["run"]]
                        startswith(key, "scan_") && (row["comparison_" * key] = value)
                    end
                end
            end
            if row["status"] == "completed"
                reference = get!(reference_cache, joinpath(folder, row["reference"])) do
                    response_values(joinpath(folder, row["reference"]))
                end
                response = response_values(joinpath(folder, name, "response.json"))
                metrics, growth = transient_metrics(response, reference; tmax = get(row, "duration_seconds", last(reference.times)))
                merge!(row, metrics)
                append!(growth_rows, [merge(Dict("run" => name, "source" => folder, "scenario" => row["scenario"]), entry) for entry in growth])
                baseline_key = (row["backend"], row["scenario"])
                baseline = if comparison !== nothing
                    get!(reference_cache, joinpath(comparison["source"], comparison["run"], "response.json")) do
                        response_values(joinpath(comparison["source"], comparison["run"], "response.json"))
                    end
                else
                    row["backend"] == "pscad" ? nothing : get(baselines, baseline_key, nothing)
                end
                if row["study"] == "passivity" && baseline !== nothing
                    changed, _ = transient_metrics(response, baseline; tmax = get(row, "duration_seconds", last(baseline.times)))
                    row["enforcement_transient_change_scaled_rms"] = changed["receiving_core_scaled_rms"]
                    row["enforcement_peak_change_kV"] = changed["receiving_core_peak_error_kV"]
                end
                if comparison !== nothing
                    # Compare both voltage and current at all terminals and times.
                    row["enforcement_transient_identical"] = isequal(response, baseline)
                    failed = occursin(r"(?i)failed", row["enforcement_outcome"])
                    if failed && row["enforcement_transient_identical"]
                        row["enforcement_simulation"] = "enforcement failed; unenforced model simulated"
                    elseif failed
                        row["enforcement_simulation"] = "enforcement failed; transient differs from matching scan run"
                    end
                end
            end
            if haskey(log_cache, name)
                log_row = log_cache[name]
                for (key, value) in row
                    (startswith(key, "enforcement_") || startswith(key, "comparison_scan_")) && (log_row[key] = value)
                end
                push!(log_rows, log_row)
                for (origin, source_run, report) in (("native", name, log_row),
                        ("matching CPASS=1", get(row, "enforcement_reference_run", ""), row))
                    prefix = origin == "native" ? "" : "comparison_"
                    for entry in get(report, prefix * "scan_ranges", Dict{String, Any}[])
                        push!(range_rows, merge(Dict("run" => name, "scenario" => row["scenario"],
                            "scan_origin" => origin, "scan_source_run" => source_run,
                            "scan_lower_Hz" => report[prefix * "scan_lower_Hz"],
                            "scan_upper_Hz" => report[prefix * "scan_upper_Hz"],
                            "FS" => row["FS"], "FE" => row["FE"]), entry))
                    end
                end
            end
            push!(rows, row)
        end
        identity = ["run", "backend", "scenario", "study", "status"]
        settings = ["FS", "FE", "Numf", "YMaxP", "AMaxP", "YMaxE", "AMaxE", "MaxRPtol", "CPASS",
            "ET_PE", "MER_PE", "MIT_PE", "FSP", "FEP", "NFP", "FDIS", "Interp1", "FT",
            "lower", "density", "order", "max_ratio", "enforce", "frequency", "interp", "divisor", "step_seconds"]
        metrics = ["receiving_core_scaled_rms", "receiving_core_scaled_peak", "receiving_core_peak_kV", "receiving_core_peak_error_kV",
            "arrival_time_us", "growth_flag", "maximum_late_peak_relative_to_scale", "poles_beyond_content", "maximum_residue_pole_ratio"]
        columns = unique(vcat(identity, settings, metrics, sort(unique([key for row in rows for key in keys(row)]))))
        write_tables(destination, "pscad_options", rows, columns;
            caption = "Option effects against exact NLT at matching sample times. Errors use the NLT receiving core peak.",
            tex_columns = vcat(identity, metrics[1:5]))
        for study in ("baseline", "band_samples", "band", "density", "residue_limit", "order_error", "order", "passivity", "time_step", "transformation")
            subset = filter(row -> string(row["study"]) == study, rows)
            write_tables(destination, "pscad_options_" * study, subset, columns;
                caption = "Option comparison: " * study * ". Blank measurements indicate an incomplete or failed run.",
                tex_columns = vcat(["run", "status"], metrics))
        end
        write_tables(destination, "pscad_options_growth", growth_rows,
            ["run", "scenario", "end", "terminal", "growth_ratio", "reference_growth_ratio",
                "late_peak_relative_to_scale", "reference_late_peak_relative_to_scale", "late_peak_excess_relative_to_scale",
                "growth_ratio_flag", "growth_flag"];
            caption = "Late-window peaks use the NLT receiving core peak as scale. Flags require excess growth and a late peak more than one percent of scale above NLT.")
        write_tables(destination, "pscad_options_passivity_ranges", range_rows,
            ["run", "scenario", "scan_origin", "scan_source_run", "scan_lower_Hz", "scan_upper_Hz",
                "FS", "FE", "lower_Hz", "upper_Hz", "smallest_eigenvalue", "relative_to_fitting_band"];
            caption = "Native violation ranges and their positions relative to the fitting band. A range can cross a band limit. Comparison scans identify their source run.")
        write_tables(destination, "pscad_options_poles", pole_rows,
            ["run", "scenario", "quantity", "group", "pole_index", "pole_real_per_s", "pole_imag_per_s", "pole_frequency_Hz", "spectral_limit_Hz", "beyond_content", "residue_pole_ratio"];
            caption = "Fitted pole frequencies are magnitudes divided by two pi. Spectral limits use the pulse's 99 percent source energy or the step's resolved NLT band.")
        write_tables(destination, "pscad_options_logs", log_rows,
            unique(vcat(["run", "source", "paths"], sort(unique([k for row in log_rows for k in keys(row)])))) ;
            caption = "Native Tline fit and passivity reports. Missing headings remain not found.",
            tex_columns = ["run", "Yc_max_error_percent", "H_max_error_percent", "H_rms_error_percent", "maximum_residue_pole_ratio", "scan_outcome", "enforcement_outcome"])
        write_tables(destination, "pscad_options_sources", attempts, ["run", "source", "status", "message", "exception_type"];
            caption = "All supplied attempts. Latest final outcomes take precedence over exports.")
        open(joinpath(destination, "pscad_options_log_reports.toml"), "w") do io
            TOML.print(io, Dict("reports" => log_rows); sorted = true)
        end
        growth_counts = Dict{String, Any}[]
        open(joinpath(destination, "pscad_options_summary.txt"), "w") do io
            for i in cable_indices
                println(io, folders[i], ": ", summaries[i])
            end
            for study in ("baseline", "band_samples", "band", "density", "residue_limit", "order_error", "order", "passivity", "time_step", "transformation"),
                    scenario in ("step_single_point", "impulse"), backend in ("pscad", "ulm", "marti")
                completed = filter(row -> row["backend"] == backend && string(row["study"]) == study && row["scenario"] == scenario &&
                    row["status"] == "completed" && haskey(row, "receiving_core_scaled_rms"), rows)
                isempty(completed) && continue
                best = completed[argmin([row["receiving_core_scaled_rms"] for row in completed])]
                before = count(row -> row["growth_ratio_flag"], completed)
                after = count(row -> row["growth_flag"], completed)
                push!(growth_counts, Dict("backend" => backend, "study" => study, "scenario" => scenario,
                    "completed" => length(completed), "growth_before" => before, "growth_after" => after))
                println(io, backend, " ", study, " ", scenario, ": smallest sampled RMS ", best["receiving_core_scaled_rms"], " at ", best["run"],
                    "; growth flagged ", before, " -> ", after, "/", length(completed))
            end
        end
        write_tables(destination, "pscad_options_growth_counts", growth_counts,
            ["backend", "study", "scenario", "completed", "growth_before", "growth_after"];
            caption = "Growth flags before and after the amplitude condition, on the same completed runs.")
        println("Option tables: ", destination, " (", length(rows), " runs)")
    end
end
