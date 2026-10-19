# Compare archived checks with archived transient outcomes.
using LineCableModels
import LineCableModels: validate
using LineCableModels.PSCAD: read_pscad_fit_report
using TOML
using EzXML
using LinearAlgebra: eigvals
using Statistics: mean

const VERIFICATION_CHECKS = (
    "fit_error_0.2pct", "fit_error_2pct", "fit_rms_2pct", "passivity", "passivity_scenario_band",
    "H0_gain", "H0_target_1pct", "Yc0_target_1pct", "ratio_1e3", "ratio_2e6", "late_growth")

csv_cell(x) = "\"" * replace(x === nothing ? "" : string(x), '"' => "\"\"") * "\""
tex_text(x) = replace(string(x), '\\' => raw"\textbackslash{}", '{' => raw"\{", '}' => raw"\}",
    '_' => raw"\_\allowbreak{}", '%' => raw"\%", '&' => raw"\&", '#' => raw"\#", '$' => raw"\$")
tex_cell(::Nothing) = "{---}"
tex_cell(x::AbstractFloat) = isfinite(x) ? "\\num[round-mode=figures,round-precision=4]{$x}" : tex_text(x)
tex_cell(x) = tex_text(x)

function write_tables(output, name, rows, columns; caption, tex_columns = columns)
    open(joinpath(output, name * ".csv"), "w") do io
        println(io, join(csv_cell.(columns), ','))
        for row in rows
            println(io, join([csv_cell(get(row, key, nothing)) for key in columns], ','))
        end
    end
    open(joinpath(output, name * ".tex"), "w") do io
        println(io, raw"\begingroup\scriptsize\begin{longtable}{", repeat("l", length(tex_columns)), "}")
        println(io, "\\caption{", tex_text(caption), raw"}\\\toprule")
        println(io, join(tex_text.(tex_columns), " & "), raw"\\\midrule\endhead")
        for row in rows
            println(io, join([tex_cell(get(row, key, nothing)) for key in tex_columns], " & "), "\\\\")
        end
        println(io, raw"\bottomrule\end{longtable}\endgroup")
    end
end

function voltage_data(observed)
    q = only(filter(q -> q.request === voltage, observed.quantities))
    return (; values = q.values, times = q.coordinates.times, labels = q.coordinates.column_labels)
end

function frequency_band(folder, coarse = false)
    path = joinpath(folder, coarse ? "coarse_baseline.json" : "baseline.json")
    observed = import_data(:observed, path)
    return extrema(first(observed.quantities).coordinates.frequencies)
end

# The impulse cutoff contains 99 percent of the source's continuous spectral energy.
# Step sources use the resolved NLT band. The ideal source has unbounded bandwidth.
function scenario_band(scenario, N, T)
    low, high = 1 / (2T), (2N - 1) / (2T)
    if occursin("impulse", scenario)
        alpha, beta = 1.47e4, 2.47e6
        lo, hi = 0.0, 100beta
        for _ in 1:80
            w = (lo + hi) / 2
            fraction = (atan(w / alpha) / alpha - atan(w / beta) / beta) / (pi / 2 * (1 / alpha - 1 / beta))
            fraction < 0.99 ? (lo = w) : (hi = w)
        end
        high = min(high, hi / (2pi))
    end
    return low, high
end

function fit_row(study, cable, fit, scenario, source, band, content)
    row = Dict{String, Any}("study" => study, "cable" => cable, "fit" => fit,
        "scenario" => scenario, "source" => source, "fitting_lower_Hz" => band[1],
        "fitting_upper_Hz" => band[2], "scenario_lower_Hz" => content[1], "scenario_upper_Hz" => content[2])
    for field in VERIFICATION_CHECKS
        row[field] = nothing
    end
    return row
end

# Classify configuration contrasts independently of transient errors.
function classify_configuration!(row; fit = String[], simulation = String[], structure = String[],
        line = String[], long_window = false)
    changes = ("fit settings" => fit, "simulation settings" => simulation,
        "model structure" => structure, "line settings" => line)
    components = [name for (name, values) in changes if !isempty(values)]
    long_window && push!(components, "long-window DC behaviour")
    row["cause_components"] = join(components, "; ")
    row["cause_class"] = isempty(components) ? "baseline" :
        length(components) == 1 ? only(components) : "mixed: " * join(components, " + ")
    row["configuration_changes"] = join(["$name: " * join(values, "; ") for (name, values) in changes if !isempty(values)], " | ")
    row["long_window"] = long_window
    return row
end

# Split reported violations at band limits. A worst value alone does not identify location.
function violation_regions!(row, ranges, range_rows)
    fs, fe = row["fitting_lower_Hz"], row["fitting_upper_Hz"]
    sl, sh = row["scenario_lower_Hz"], row["scenario_upper_Hz"]
    for (low, high, worst) in ranges
        regions = String[]
        low < fs && push!(regions, "below fitting band")
        sl > fs && low < sl && max(low, fs) <= min(high, fe, sl) && push!(regions, "inside fitting band below scenario band")
        max(low, fs, sl) <= min(high, fe, sh) && push!(regions, "inside scenario band")
        high > sh && fe > sh && max(low, fs, sh) <= min(high, fe) && push!(regions, "inside fitting band above scenario band")
        high > fe && push!(regions, "above fitting band")
        push!(range_rows, merge(Dict(key => row[key] for key in ("study", "cable", "fit", "scenario", "source")),
            Dict("lower_Hz" => low, "upper_Hz" => high, "smallest_eigenvalue" => worst,
                "location" => join(regions, "; "))))
        "inside scenario band" in regions && (row["passivity_scenario_band"] = true)
    end
end

function julia_checks!(row, fit, range_rows)
    row["method"] = string(get(fit, :method, get(get(fit, :variant, (;)), :method, :ulm)))
    errors = fit.fitting_error
    for (section, record) in pairs(errors)
        row["$(section)_max_error_percent"] = 100maximum(record.max)
        row["$(section)_rms_error_percent"] = 100maximum(record.rms)
    end
    row["reported_max_fit_error_percent"] = maximum(100maximum(record.max) for record in errors)
    report = get(fit, :passivity_after, get(fit, :passivity, nothing))
    if report !== nothing && report.checked
        row["passivity"] = !report.passive
        row["smallest_eigenvalue"] = report.worst
        row["passivity_scenario_band"] = false
        ranges = [(band[1], band[2], nothing) for band in report.bands]
        violation_regions!(row, ranges, range_rows)
        row["passivity_evidence"] = "archived fit passivity record"
        row["ratio_source"] = "archived propagation ratio summary, when present"
        row["maximum_residue_pole_ratio"] = get(fit, :largest_residue_pole_ratio,
            report.perturbation === nothing ? nothing : maximum(report.perturbation.res_pole_ratio; init = 0.0))
    end
    # Only dimensionless propagation ratios of residue to pole enter the ratio checks.
    # Admittance and impedance ratios have units and cannot use the same thresholds.
    if haskey(fit, :poles)
        ratios = Float64[]
        groups = row["method"] == "marti" ? ((length(fit.poles) ÷ 2 + 1):length(fit.poles)) : (2:length(fit.poles))
        for group in groups, (k, pole) in enumerate(fit.poles[group])
            residue = row["method"] == "marti" ? abs(fit.residues[group][k]) : maximum(abs, fit.residues[group][:, :, k])
            push!(ratios, residue / abs(pole))
        end
        row["ratio_source"] = "archived propagation poles and residues"
        row["maximum_residue_pole_ratio"] = maximum(ratios; init = 0.0)
    elseif haskey(fit, :model)
        row["ratio_source"] = "archived propagation poles and residues"
        row["maximum_residue_pole_ratio"] = maximum(maximum(abs, group.residues[:, :, k]) / abs(pole)
            for group in fit.model.propagation for (k, pole) in enumerate(group.poles))
    end
    dc = get(fit, :dc, (; targets = get(fit, :targets, nothing), achieved = get(fit, :achieved, nothing)))
    if dc.achieved !== nothing
        h = dc.achieved.propagation
        row["H0_spectral_radius"] = maximum(abs, eigvals(h))
        row["H0_maximum_entry"] = maximum(abs, h)
        row["H0_gain"] = row["H0_spectral_radius"] > 1
    end
    if dc.targets !== nothing && dc.achieved !== nothing
        for (key, field) in (("H0", :propagation), ("Yc0", :admittance))
            target, actual = getproperty(dc.targets, field), getproperty(dc.achieved, field)
            row[key * "_relative_deviation"] = maximum(abs, actual - target) / maximum(abs, target)
        end
    end
end

function native_checks!(row, directory, status, range_rows)
    controls = copy(get(status, "exported_fit", get(status, "fit", Dict{String, Any}())))
    local_paths = sort([joinpath(root, file) for (root, _, files) in walkdir(directory) for file in files])
    projects = filter(path -> endswith(path, ".pscx"), local_paths)
    for path in projects
        document = readxml(path)
        for node in findall("//User[@defn='master:Line_FrePhase_Options' or @defn='master:Line_FreDep_Options']/paramlist/param", document)
            controls[node["name"]] = parse(Float64, node["value"])
        end
    end
    row["fitting_lower_Hz"] = get(controls, "FS", row["fitting_lower_Hz"])
    row["fitting_upper_Hz"] = get(controls, "FE", row["fitting_upper_Hz"])
    row["fit_controls"] = join(["$key=$(controls[key])" for key in sort(collect(keys(controls)))], "; ")
    row["method"] = get(status, "model", "fd_phase")
    paths = filter(path -> endswith(lowercase(path), ".log"), local_paths)
    row["native_log_count"] = length(paths)
    reports = read_pscad_fit_report.(paths)
    for report in reports
        if report.reports.fitting !== nothing
            for (section, record) in (("admittance", report.admittance), ("propagation", report.propagation))
                record === nothing && continue
                final = something(record.after_dc, record.before_dc, (; maximum_error_percent = nothing, rms_error_percent = nothing))
                row[section * "_max_error_percent"] = final.maximum_error_percent
                row[section * "_rms_error_percent"] = final.rms_error_percent
            end
            row["maximum_residue_pole_ratio"] = report.maximum_residue_pole_ratio
            row["ratio_source"] = "native final fit report"
        end
        if report.modal !== nothing
            for (section, quantity) in (("impedance", :impedance), ("propagation", :propagation))
                modal_rows = filter(record -> record.quantity === quantity, report.modal)
                isempty(modal_rows) || (row[section * "_native_reported_error"] = maximum(record.reported_error for record in modal_rows))
            end
            row["native_modal_error_units"] = "not stated in log"
        end
        scan = report.scan
        if scan !== nothing
            row["passivity"] = isempty(scan.ranges) ? (scan.outcome === nothing ? nothing : false) : true
            row["smallest_eigenvalue"] = scan.smallest_eigenvalue
            row["passivity_evidence"] = "native scan"
            row["scan_outcome"] = scan.outcome
            row["passivity_scenario_band"] = row["passivity"] === nothing ? nothing : false
            if scan.band !== nothing
                row["scan_lower_Hz"], row["scan_upper_Hz"] = scan.band.lower_Hz, scan.band.upper_Hz
            end
            violation_regions!(row, [(r.lower_Hz, r.upper_Hz, r.smallest_eigenvalue) for r in scan.ranges], range_rows)
        end
        enforcement = report.enforcement
        if enforcement !== nothing
            row["enforcement_outcome"] = enforcement.outcome
            row["enforcement_initial_smallest_eigenvalue"] = enforcement.initial_smallest_eigenvalue
            row["enforcement_final_smallest_eigenvalue"] = enforcement.final_smallest_eigenvalue
            # The successful outcome supersedes negative samples from earlier iterations.
            succeeded = enforcement.outcome !== nothing && occursin("no passivity violations", enforcement.outcome)
            row["passivity"] = succeeded ? false : enforcement.initial_smallest_eigenvalue === nothing ? nothing : enforcement.initial_smallest_eigenvalue < 0
            row["smallest_eigenvalue"] = succeeded ? nothing : enforcement.initial_smallest_eigenvalue
            row["passivity_evidence"] = succeeded ? "native enforcement succeeded" :
                enforcement.outcome == "Constrained solution failed" ? "failed enforcement; initial samples of unenforced model" : "initial enforcement samples; outcome unknown"
            row["passivity_scenario_band"] = row["passivity"] === nothing ? nothing : false
            if !succeeded
                violation_regions!(row, [(r.lower_Hz, r.upper_Hz, r.smallest_eigenvalue) for r in enforcement.ranges], range_rows)
                initial = enforcement.initial_samples
                if initial !== nothing
                    violation_regions!(row, [(f, f, v) for (f, v) in zip(initial.frequencies, initial.eigenvalues) if v < 0], range_rows)
                end
            end
        end
    end
    values = [get(row, section * "_max_error_percent", nothing) for section in ("impedance", "admittance", "propagation")]
    present = filter(!isnothing, values)
    isempty(present) || (row["reported_max_fit_error_percent"] = maximum(present))
    return controls
end

function transient_outcomes!(row, observed, reference, growth_rows, sheath_rows; endpoint = false)
    response, target = voltage_data(observed), voltage_data(reference)
    response.labels == target.labels || error("Terminal order differs: $(row["source"])")
    length(response.times) == length(target.times) &&
        all(isapprox.(response.times, target.times; rtol = 2e-6, atol = 1e-10)) || error("NLT output time grid differs: $(row["source"])")
    core = only(findall(==("cable:1:core"), response.labels))
    model_core, exact_core = response.values[2, core, :], target.values[2, core, :]
    scale = maximum(abs, exact_core)
    row["reference_core_peak_kV"] = scale / 1000
    sheath_terminals = findall(label -> endswith(label, ":sheath"), response.labels)
    reference_sheath_peaks = [maximum(abs, target.values[2, terminal, :]) for terminal in sheath_terminals]
    sheath_scale = maximum(reference_sheath_peaks)
    row["reference_receiving_sheath_peak_scale_kV"] = sheath_scale / 1000
    row["nonfinite_output"] = any(q -> !all(isfinite, q.values), observed.quantities)
    row["nonfinite_output"] && return
    model_scale = maximum(abs, model_core)
    row["receiving_core_peak_kV"] = model_scale / 1000
    row["receiving_core_peak_error_percent"] = 100(model_scale - scale) / scale
    row["receiving_core_scaled_rms"] = sqrt(mean(abs2, model_core - exact_core)) / scale
    sheath_errors = Float64[]
    for (terminal, expected) in zip(sheath_terminals, reference_sheath_peaks)
        label = response.labels[terminal]
        actual = maximum(abs, response.values[2, terminal, :])
        deviation = iszero(expected) ? nothing : 100(actual - expected) / expected
        scaled = iszero(sheath_scale) ? nothing : 100abs(actual - expected) / sheath_scale
        scaled === nothing || push!(sheath_errors, scaled)
        if label == replace(response.labels[core], ":core" => ":sheath")
            row["receiving_energized_sheath_peak_error_percent"] = deviation
        end
        push!(sheath_rows, merge(Dict(k => row[k] for k in ("study", "cable", "fit", "scenario", "source")),
            Dict("terminal" => label, "peak_kV" => actual / 1000, "reference_peak_kV" => expected / 1000,
                "peak_error_percent" => deviation, "scaled_peak_error_percent" => scaled,
                "reference_peak_scale_kV" => sheath_scale / 1000)))
    end
    if !isempty(sheath_errors)
        row["receiving_sheath_scaled_peak_error_percent"] = maximum(sheath_errors)
    end
    tmax = get(row, "duration_seconds", min(last(response.times), last(target.times)))
    row["late_growth"] = false
    row["reference_qualified_growth"] = false
    for side in 1:2, (terminal, label) in enumerate(response.labels)
        peaks = map((response, target)) do run
            previous = findall(t -> 0.8tmax <= t < 0.9tmax, run.times)
            final = findall(t -> 0.9tmax <= t <= tmax, run.times)
            a, b = maximum(abs, run.values[side, terminal, previous]), maximum(abs, run.values[side, terminal, final])
            ratio = iszero(a) ? (iszero(b) ? 1.0 : Inf) : b / a
            return (; previous = a, late = b, ratio)
        end
        actual, exact = peaks
        independent = actual.ratio > 1.1 && actual.late - actual.previous > 0.01model_scale
        qualified = actual.ratio > 1 && actual.ratio > 1.1exact.ratio && actual.late - exact.late > 0.01scale
        row["late_growth"] |= independent
        row["reference_qualified_growth"] |= qualified
        push!(growth_rows, merge(Dict(k => row[k] for k in ("study", "cable", "fit", "scenario", "source")),
            Dict("end" => side == 1 ? "sending" : "receiving", "terminal" => label,
                "growth_ratio" => actual.ratio, "reference_growth_ratio" => exact.ratio,
                "late_peak_relative_to_scale" => actual.late / scale,
                "late_peak_excess_relative_to_scale" => (actual.late - exact.late) / scale,
                "late_growth" => independent, "reference_qualified_growth" => qualified)))
    end
    if endpoint
        tail = findall(>=(0.9tmax), response.times)
        row["DC_endpoint_error_percent"] = 100(mean(model_core[tail]) - mean(exact_core[tail])) / mean(exact_core[tail])
    end
end

function complete_checks!(row)
    rms = filter(!isnothing, [get(row, section * "_rms_error_percent", nothing) for section in ("admittance", "propagation")])
    isempty(rms) || (row["reported_rms_fit_error_percent"] = maximum(rms))
    for (check, field, threshold) in (("fit_error_0.2pct", "reported_max_fit_error_percent", 0.2),
            ("fit_error_2pct", "reported_max_fit_error_percent", 2),
            ("fit_rms_2pct", "reported_rms_fit_error_percent", 2),
            ("H0_target_1pct", "H0_relative_deviation", 0.01), ("Yc0_target_1pct", "Yc0_relative_deviation", 0.01),
            ("ratio_1e3", "maximum_residue_pole_ratio", 1e3), ("ratio_2e6", "maximum_residue_pole_ratio", 2e6))
        value = get(row, field, nothing)
        row[check] = value === nothing ? nothing : value > threshold
    end
    row["divergence"] = get(row, "nonfinite_output", false) || get(row, "reference_qualified_growth", false)
    for (name, core, sheath) in (("strict", 2, 10), ("loose", 5, 20))
        row["failed_" * name] = row["divergence"] || any(get(row, key, nothing) !== nothing && abs(row[key]) > limit
            for (key, limit) in (("receiving_core_peak_error_percent", core),
                ("receiving_sheath_scaled_peak_error_percent", sheath), ("DC_endpoint_error_percent", 5)))
    end
    row["any_check"] = any(get(row, key, nothing) === true for key in VERIFICATION_CHECKS)
    row["any_check_fit2pct"] = any(get(row, key, nothing) === true for key in VERIFICATION_CHECKS if key != "fit_error_0.2pct")
    row["any_check_rms2pct"] = any(get(row, key, nothing) === true for key in VERIFICATION_CHECKS
        if !(key in ("fit_error_0.2pct", "fit_error_2pct")))
    row["unchecked_checks"] = join([key for key in VERIFICATION_CHECKS if row[key] === nothing], "; ")
end

function detection_counts(rows)
    detection, compact, silent = [Dict{String, Any}[] for _ in 1:3]
    classes = sort(unique(row["cause_class"] for row in rows))
    for cause in ("pooled", classes...), definition in ("strict", "loose")
        subset = cause == "pooled" ? rows : filter(row -> row["cause_class"] == cause, rows)
        all_failed = count(row -> row["failed_" * definition], subset)
        all_fine = length(subset) - all_failed
        for check in (VERIFICATION_CHECKS..., "any_check", "any_check_fit2pct", "any_check_rms2pct")
            data = Dict{String, Any}("cause_class" => cause, "definition" => definition, "check" => check)
            for (column, flag, failed) in (("flagged_and_failed", true, true), ("flagged_and_fine", true, false),
                    ("missed_failures", false, true), ("quiet_and_fine", false, false),
                    ("unavailable_failed", nothing, true), ("unavailable_fine", nothing, false))
                data[column] = count(row -> get(row, check, nothing) === flag && row["failed_" * definition] == failed, subset)
            end
            push!(detection, data)
            push!(compact, merge(data, Dict("all_failed" => all_failed, "all_fine" => all_fine,
                "sensitivity" => iszero(all_failed) ? nothing : data["flagged_and_failed"] / all_failed,
                "specificity" => iszero(all_fine) ? nothing : data["quiet_and_fine"] / all_fine,
                "sensitivity_count" => "$(data["flagged_and_failed"])/$all_failed",
                "specificity_count" => "$(data["quiet_and_fine"])/$all_fine")))
            if cause != "pooled" && startswith(check, "any_check")
                append!(silent, [merge(row, Dict("definition" => definition, "check_set" => check)) for row in subset
                    if row["failed_" * definition] && !row[check]])
            end
        end
    end
    return (; detection, compact, silent)
end

function verification_tables(folders, output)
    validate(folders, verification_tables, output)
    rows, ranges, growth, sheaths, coverage = [Dict{String, Any}[] for _ in 1:5]
    dc_folders = String[]
    for folder in folders
        cable = basename(dirname(folder))
        println("Reading ", folder); flush(stdout)
        if isfile(joinpath(folder, "dc_behaviour.toml"))
            push!(dc_folders, folder)
        elseif isfile(joinpath(folder, "pscad_options.toml"))
            metadata = TOML.parsefile(joinpath(folder, "pscad_options.toml"))
            native_baseline = TOML.parsefile(joinpath(folder, "options_pscad_impulse_1/status.toml"))["fit"]
            julia_baseline = TOML.parsefile(joinpath(folder, "options_fit_1/status.toml"))
            refs, fits = Dict{String, Any}(), Dict{String, Any}()
            for name in readdir(folder)
                occursin(r"^options_(?:julia|pscad)_", name) || continue
                directory = joinpath(folder, name)
                status = TOML.parsefile(joinpath(directory, "status.toml"))
                state = status["status"]
                push!(coverage, Dict("source" => directory, "study" => "pscad_options", "status" => state,
                    "included" => state in ("completed", "nonfinite_response")))
                state in ("completed", "nonfinite_response") || continue
                scenario = only(s for s in metadata["scenarios"] if s["name"] == status["scenario"])
                content = scenario_band(status["scenario"], scenario["N"], scenario["Tfactor"] * scenario["duration_seconds"])
                row = fit_row("pscad_options", cable, name, status["scenario"], directory, (0.5, 1e6), content)
                row["status"] = state
                row["duration_seconds"] = status["duration_seconds"]
                fit_changes, simulation_changes, structure_changes = String[], String[], String[]
                status["divisor"] == 0 || push!(simulation_changes, "time-step divisor=$(status["divisor"])")
                row["time_step_divisor"] = status["divisor"]
                row["step_seconds"] = status["step_seconds"]
                if status["backend"] == "pscad"
                    native_checks!(row, directory, status, ranges)
                    row["interpolation"] = get(status["fit"], "Interp1", 1)
                    row["interpolation"] == 1 || push!(simulation_changes, "Interp1=$(row["interpolation"])")
                    if status["model"] == "fd_mode"
                        push!(structure_changes, "FD-Mode constant transformation; FT=$(status["fit"]["FT"])")
                    else
                        append!(fit_changes, ["$key=$value (baseline $(native_baseline[key]))"
                            for (key, value) in sort(collect(status["fit"]); by = first)
                            if key != "Interp1" && value != native_baseline[key]])
                    end
                else
                    fit = get!(fits, status["fit_archive"]) do
                        import_data(:observed, joinpath(folder, status["fit_archive"])).gridpoint.options_fit
                    end
                    row["fitting_lower_Hz"], row["fitting_upper_Hz"] = extrema(fit.frequencies)
                    row["fit_controls"] = string(fit.variant)
                    julia_checks!(row, fit, ranges)
                    row["interpolation"] = status["interp"]
                    status["interp"] == "single" || push!(simulation_changes, "interp=$(status["interp"])")
                    append!(fit_changes, ["$key=$value (baseline $(julia_baseline[string(key)]))"
                        for (key, value) in pairs(fit.variant)
                        if !(key in (:study, :method, :frequency)) && value != julia_baseline[string(key)]])
                    fit.variant.method === :marti && push!(structure_changes,
                        "J. Marti constant transformation; frequency=$(fit.variant.frequency)")
                end
                classify_configuration!(row; fit = fit_changes, simulation = simulation_changes, structure = structure_changes)
                if state == "completed"
                    reference = get!(refs, status["reference"]) do
                        import_data(:observed, joinpath(folder, status["reference"]))
                    end
                    transient_outcomes!(row, import_data(:observed, joinpath(directory, "response.json")), reference, growth, sheaths)
                else
                    row["nonfinite_output"] = true
                end
                push!(rows, row)
            end
        elseif isfile(joinpath(folder, "transient_sensitivity.json"))
            points = import_data(:observed, joinpath(folder, "transient_sensitivity.json"))
            for observed in points
                point = observed.gridpoint.sensitivity
                push!(coverage, Dict("study" => "transient_sensitivity", "source" => folder,
                    "fit" => "$(point.scenario)/$(point.study)/$(point.input)/$(point.level)", "status" => point.status,
                    "included" => false, "reason" => "physical-input NLT response; no fitted model or independent model error"))
            end
        elseif isfile(joinpath(folder, "pscad_station.toml"))
            references = import_data(:observed, joinpath(folder, "references.json"))
            for name in readdir(folder)
                occursin(r"^pscad_(?:impulse|step_single_point)_\d+$", name) || continue
                directory = joinpath(folder, name)
                status = TOML.parsefile(joinpath(directory, "status.toml"))
                state = status["status"]
                push!(coverage, Dict("study" => "time_domain", "source" => directory, "status" => state,
                    "included" => state in ("completed", "nonfinite_response")))
                state in ("completed", "nonfinite_response") || continue
                reference = only(filter(p -> p.gridpoint.study.scenario == Symbol(status["scenario"]) && p.gridpoint.study.resistance == 10, references))
                evidence = reference.gridpoint.study
                row = fit_row("time_domain", cable, name, status["scenario"], directory, (0.5, 1e6),
                    scenario_band(status["scenario"], evidence.N, evidence.T))
                row["status"] = state
                row["duration_seconds"] = evidence.window.tmax
                controls = native_checks!(row, directory, status, ranges)
                fit_changes = row["method"] == "fd_mode" ? String[] :
                    ["$key=$(controls[key]) (baseline $value)" for (key, value) in
                        ("YMaxP" => 20, "AMaxP" => 20, "DCenab" => 0, "CPASS" => 1) if controls[key] != value]
                classify_configuration!(row; fit = fit_changes,
                    structure = row["method"] == "fd_mode" ? ["FD-Mode constant transformation"] : String[])
                state == "completed" ? transient_outcomes!(row, import_data(:observed, joinpath(directory, "response.json")), reference, growth, sheaths) : (row["nonfinite_output"] = true)
                push!(rows, row)
            end
        elseif isfile(joinpath(folder, "line_models.json"))
            fits = import_data(:observed, joinpath(folder, "line_models.json"))
            settings = isfile(joinpath(folder, "model_settings.txt"))
            references = import_data(:observed, joinpath(folder, settings ? "reference.json" : "references.json"))
            for (index, observed) in enumerate(fits)
                fit = getproperty(observed.gridpoint, settings ? :model_settings : :study)
                band = frequency_band(folder, settings ? fit.method === :ulm : fit.line === :coarse)
                for scenario in (settings ? ("step_single_point",) : ("step_single_point", "step_both_ends", "impulse"))
                    reference = settings ? references : only(filter(p -> string(p.gridpoint.study.scenario) == scenario && p.gridpoint.study.resistance == 10, references))
                    evidence = settings ? nothing : reference.gridpoint.study
                    content = settings ? scenario_band(scenario, 4096, 1.25fit.window.tmax) : scenario_band(scenario, evidence.N, evidence.T)
                    for interp in (settings ? (fit.method === :marti ? ("single",) : ("single", "double")) : ("single",))
                        file = "$(scenario)_response_$(index)$(settings ? "_" * interp : "").json"
                        path = joinpath(folder, file)
                        isfile(path) || continue
                        identity = "fit_$(index)$(settings ? "_" * interp : "")"
                        row = fit_row(settings ? "model_settings" : "time_domain", cable, identity, scenario, folder, band, content)
                        row["response_file"] = file
                        row["status"] = "completed"
                        row["duration_seconds"] = settings ? fit.window.tmax : evidence.window.tmax
                        row["fit_controls"] = string(fit.options)
                        julia_checks!(row, fit, ranges)
                        fit_changes, structure_changes = String[], String[]
                        if fit.method === :ulm
                            for (key, actual, nominal) in (("Yc order", fit.options.yc.npoles, 16),
                                    ("H order", fit.options.h.npoles, 16), ("H time rows", fit.options.h.td, true),
                                    ("passivity enforcement", fit.options.passivity, true), ("DC correction", fit.options.dc, (;)),
                                    ("maximum residue/pole ratio", get(fit.options.h, :max_ratio, Inf), Inf))
                                actual == nominal || push!(fit_changes, "$key=$actual (baseline $nominal)")
                            end
                        else
                            push!(structure_changes, fit.method === :marti ?
                                "J. Marti constant transformation; $(fit.options.transformation)" : "FDCM formulation")
                        end
                        (settings ? fit.method === :ulm : fit.line === :coarse) &&
                            push!(fit_changes, "coarse sampling: 101 logarithmic points from 1 Hz to 1 MHz plus 50 Hz")
                        row["interpolation"] = interp
                        classify_configuration!(row; fit = fit_changes,
                            simulation = interp == "single" ? String[] : ["interp=$interp"], structure = structure_changes)
                        transient_outcomes!(row, import_data(:observed, path), reference, growth, sheaths)
                        push!(rows, row)
                    end
                end
            end
        else
            error("No supported archived study in $folder")
        end
        GC.gc()
    end
    # Consolidate reruns by cable and run name in argument order, as in the DC generator.
    selected = Dict{Tuple{String, String}, Any}()
    summaries = Dict(folder => TOML.parsefile(joinpath(folder, "dc_behaviour.toml")) for folder in dc_folders)
    for folder in dc_folders, name in readdir(folder)
        occursin(r"^dc_(?:nlt|ulm_\d+|pscad_\d+)$", name) || continue
        path = joinpath(folder, name, "status.toml")
        isfile(path) || continue
        status = TOML.parsefile(path)
        cable = summaries[folder]["cable"]
        key = (cable, name)
        final = status["status"] in ("completed", "nonfinite_response", "line_constants_error", "known_station_failure", "julia_error", "unsettled_reference")
        if !haskey(selected, key) || final
            selected[key] = (; folder, name, status)
        end
    end
    for cable in unique(first(key) for key in keys(selected))
        reference_source = selected[(cable, "dc_nlt")]
        long_ref = import_data(:observed, joinpath(reference_source.folder, "dc_nlt/response.json"))
        short_refs = Dict{String, Any}()
        for key in sort(collect(keys(selected)))
            first(key) == cable && last(key) != "dc_nlt" || continue
            (; folder, name, status) = selected[key]
            directory = joinpath(folder, name)
            state = status["status"]
            push!(coverage, Dict("study" => "dc_behaviour", "source" => directory, "status" => state,
                "included" => state in ("completed", "nonfinite_response")))
            state in ("completed", "nonfinite_response") || continue
            short = get(status, "scenario", "long") == "short"
            scenario = short ? "step_short" : "step_long"
            reference = short ? get!(short_refs, folder) do
                import_data(:observed, joinpath(folder, "dc_short_reference.json"))
            end : long_ref
            duration = last(voltage_data(reference).times)
            row = fit_row("dc_behaviour", cable, name, scenario, directory, (0.001, 1e6),
                short ? scenario_band(scenario, 4096, 1.25duration) : scenario_band(scenario, 8192, 2duration))
            row["status"] = state
            row["duration_seconds"] = get(status, "duration_seconds", duration)
            fit_changes, line_changes = String[], String[]
            if startswith(name, "dc_pscad")
                controls = native_checks!(row, directory, status, ranges)
                append!(fit_changes, ["$key=$(controls[key]) (baseline $value)" for (key, value) in
                    ("FS" => 0.5, "YMaxP" => 20, "AMaxP" => 20, "DCenab" => 0, "CPASS" => 0)
                    if controls[key] != value])
                get(status, "line_formulation", "lossy") == "lossy" || push!(line_changes, "lossless dielectric formulation")
                get(controls, "enablf", 1) == 1 || push!(line_changes,
                    "constant shunt conductance=$(controls["shntcab"]) S/m instead of loss tangent")
                row["line_controls"] = "$(get(status, "line_formulation", "lossy")); enablf=$(get(controls, "enablf", 1)); shntcab=$(get(controls, "shntcab", 0))"
            else
                fit = import_data(:observed, joinpath(directory, "fit.json")).gridpoint.dc_behaviour
                row["fit_controls"] = "dc_method=$(fit.method); passivity=$(fit.enforce)"
                julia_checks!(row, fit, ranges)
                row["method"] = "ulm"
                row["dc_method"] = string(fit.method)
                fit.method === :none || push!(fit_changes, "DC correction=$(fit.method)")
                fit.enforce && push!(fit_changes, "passivity enforcement=true (baseline false)")
            end
            classify_configuration!(row; fit = fit_changes, line = line_changes, long_window = !short)
            if state == "completed"
                transient_outcomes!(row, import_data(:observed, joinpath(directory, "response.json")), reference, growth, sheaths; endpoint = !short)
            else
                row["nonfinite_output"] = true
            end
            push!(rows, row)
            GC.gc(false)
        end
    end
    foreach(complete_checks!, rows)
    sort!(rows; by = row -> (row["study"], row["cable"], row["fit"], row["scenario"]))
    mkpath(output)
    fields = sort(unique(vcat([collect(keys(row)) for row in rows]...)))
    identity = ["study", "cable", "fit", "scenario", "method", "source"]
    write_tables(output, "verification_fits", rows, vcat(identity, setdiff(fields, identity));
        caption = "One archived fitted configuration and scenario per row. Empty fields mean unavailable. Errors are percent; scaled RMS is a fraction.",
        tex_columns = ["cause_class", identity[1:5]..., "receiving_core_peak_error_percent", "receiving_energized_sheath_peak_error_percent", "receiving_sheath_scaled_peak_error_percent", "DC_endpoint_error_percent", "divergence", "failed_strict", "failed_loose"])
    (; detection, compact, silent) = detection_counts(rows)
    classes = [Dict("cause_class" => r["cause_class"], "definition" => r["definition"],
        "rows" => r["all_failed"] + r["all_fine"], "failed" => r["all_failed"], "fine" => r["all_fine"],
        "silent" => r["missed_failures"]) for r in compact if r["check"] == "any_check"]
    write_tables(output, "verification_classes", classes,
        ["cause_class", "definition", "rows", "failed", "fine", "silent"];
        caption = "Configuration classes and silent failures with the default combined check set. Mixed classes are disjoint from the individual classes.")
    write_tables(output, "verification_detection", detection,
        ["cause_class", "definition", "check", "flagged_and_failed", "flagged_and_fine", "missed_failures", "quiet_and_fine", "unavailable_failed", "unavailable_fine"];
        caption = "Detection counts use only available checks. Unavailable checks have separate failure and fine counts. Fine means no observed threshold exceedance in retained outcomes.")
    write_tables(output, "verification_power", compact,
        ["cause_class", "definition", "check", "sensitivity_count", "specificity_count", "sensitivity", "specificity",
            "all_failed", "all_fine", "unavailable_failed", "unavailable_fine"];
        caption = "Sensitivity is flagged-failed / all-failed; specificity is quiet-fine / all-fine. Unavailable checks remain in the denominators and count in neither numerator. A zero denominator has no rate.",
        tex_columns = ["cause_class", "definition", "check", "sensitivity_count", "specificity_count", "unavailable_failed", "unavailable_fine"])
    write_tables(output, "verification_silent", silent,
        ["cause_class", "definition", "check_set", identity..., "reported_max_fit_error_percent", "reported_rms_fit_error_percent", "receiving_core_peak_error_percent", "receiving_core_scaled_rms", "receiving_energized_sheath_peak_error_percent", "receiving_sheath_scaled_peak_error_percent", "DC_endpoint_error_percent", "divergence", "unchecked_checks"];
        caption = "Failures with no available check flag. Missing checks are listed; silence does not imply every check passed.")
    write_tables(output, "verification_mixed", filter(row -> startswith(row["cause_class"], "mixed:"), rows),
        ["cause_class", "configuration_changes", identity..., "failed_strict", "failed_loose"];
        caption = "Configurations changing more than one cause class. These rows have their own detection counts and do not enter either single-class count.")
    for (name, data, columns) in (("verification_ranges", ranges, [identity[1:4]..., "lower_Hz", "upper_Hz", "smallest_eigenvalue", "location", "source"]),
            ("verification_growth", growth, [identity[1:4]..., "end", "terminal", "growth_ratio", "reference_growth_ratio", "late_peak_relative_to_scale", "late_peak_excess_relative_to_scale", "late_growth", "reference_qualified_growth", "source"]),
            ("verification_sheaths", sheaths, [identity[1:4]..., "terminal", "peak_kV", "reference_peak_kV", "peak_error_percent", "scaled_peak_error_percent", "reference_peak_scale_kV", "source"]),
            ("verification_coverage", coverage, ["study", "source", "fit", "status", "included", "reason"]))
        write_tables(output, name, data, columns; caption = "Archived evidence and availability for the verification comparison.")
    end
    open(joinpath(output, "verification_summary.txt"), "w") do io
        println(io, "Failure: nonfinite output or reference-qualified growth; or absolute receiving-core peak error >2%, scaled receiving-sheath peak error >10%, or DC endpoint error >5%.")
        println(io, "Loose definition: core >5%, scaled sheath >20%; divergence and DC threshold unchanged.")
        println(io, "Each sheath peak is the maximum absolute receiving-end voltage over the scenario window. Energized-sheath error is the signed model-minus-reference peak divided by that sheath's reference peak, in percent. The energized sheath belongs to the cable with the driven core.")
        println(io, "Scaled sheath error is the largest absolute model-minus-reference peak difference across receiving sheaths, divided by the largest reference receiving-sheath peak of the scenario, in percent. Only this common-scale sheath measure enters failure criteria. A zero reference denominator leaves its measure unavailable.")
        println(io, "Fit-error checks use maximum reported entry error >0.2% and >2%, and reported RMS >2%. The any_check set includes maximum >0.2%; any_check_fit2pct uses maximum >2%; any_check_rms2pct uses RMS >2%.")
        println(io, "FD-Mode raw fitting errors have no unit label in the log. The native values are retained separately; normalized percent-threshold checks remain unavailable.")
        println(io, "DC-target checks use maximum absolute matrix difference / maximum target magnitude >1%. H0 gain uses spectral radius >1. Only archived targets/achieved values are used.")
        println(io, "Ratios use dimensionless propagation residues only. J. Marti passivity remains unchecked. Missing evidence is empty, never a passing check.")
        println(io, "Standalone late growth: peak ratio >1.1 and late-minus-previous peak >1% of model receiving-core peak, at any voltage terminal.")
        println(io, "Reference-qualified growth outcome: ratio >1 and >1.1 times NLT ratio, and late peak exceeds NLT by >1% of NLT receiving-core peak. This outcome needs NLT and is not an independent check.")
        println(io, "Both growth rules compare [80%,90%) and [90%,100%] of the requested scenario duration. DC endpoint uses last-10% voltage means, including finite unsettled runs.")
        println(io, "Scenario band: NLT lower frequency to NLT upper frequency, capped for impulses at 99% source spectral energy. Violations spanning boundaries retain every intersected region; band minima are not invented.")
        println(io, "Passivity uses final archived fit records; native successful enforcement supersedes earlier negative samples. Native failed enforcement uses initial samples of the unenforced model. Eigenvalue normalizations differ by backend.")
        println(io, "Physical-input transient_sensitivity points are exact NLT responses, not fitted-model outcomes, and are listed in coverage only. Fit/configuration failures without transient outcomes are also excluded.")
        println(io, "Rows retain model settings and time-step/interpolation configurations separately. Repeated DC attempts consolidate in argument order; final attempts replace earlier attempts.")
        println(io, "Cause classes describe configuration contrasts, not established physical causes. Classification never uses an outcome or a diagnostic value.")
        println(io, "time_domain Julia baseline: dense baseline grid, ULM orders 16/16, H time rows enabled, passivity enforcement enabled, no DC correction, single interpolation, original scenario window. Order/time-row/passivity/DC/sampling changes are fit settings. Marti and FDCM are model-structure changes; FDCM is not labelled a constant transformation. Coarse Marti is mixed fit settings + model structure.")
        println(io, "time_domain PSCAD baseline: FD-Phase orders 20/20, DCenab=0, CPASS=1, unchanged fitting band, original scenario window and travel interpolation. Changed order/DC/passivity controls are fit settings. FD-Mode is model structure.")
        println(io, "pscad_options baseline: the archived variant 1 per backend, divisor=0 and native Interp1=1 or Julia interp=single. Compare fitting controls to variant 1; ignore sweep names, so duplicate default configurations remain baseline. Divisor/interpolation changes are simulation settings. FD-Mode/Marti and their transformation frequencies are model structure. Changes in multiple fitting controls remain one fit-settings class.")
        println(io, "model_settings baseline: the same dense ULM 16/16 fit as time_domain with max_ratio=Inf and single interpolation. Both ULM fits use the coarse 1 Hz grid, hence fit settings; double interpolation adds simulation settings. Marti transformation-frequency choices are model structure on the dense grid.")
        println(io, "dc_behaviour baseline: lossy dielectric formulation, native FS=0.5, orders 20/20, DCenab=0, CPASS=0, loss tangent; ULM orders 16/16, H time rows enabled, no DC correction, no passivity enforcement. Short scenarios retain their original window. Every long scenario adds long-window DC behaviour; DC/passivity/band changes add fit settings, and lossless dielectric or constant shunt-conductance overrides add line settings. These line changes are explicitly outside the five requested single classes.")
        println(io, "All simultaneous class changes have a distinct mixed class and are listed in verification_mixed.csv. No mixed row enters a single-class detection count. Baseline means no contrast against the stated study/backend baseline, not identical defaults across studies.")
        println(io, "Sensitivity = flagged-failed / all-failed; specificity = quiet-fine / all-fine. Missing checks stay in denominators but in neither numerator, so these are observed detection/clearance fractions, not estimates conditioned on availability. Zero denominators have no rate. Counts and unavailable evidence are in verification_power.csv.")
        println(io, "Total fitted configuration/scenario rows: ", length(rows))
        for study in sort(unique(row["study"] for row in rows))
            subset = filter(row -> row["study"] == study, rows)
            println(io, study, ": ", length(subset), "; failures strict/loose = ", count(row -> row["failed_strict"], subset), "/", count(row -> row["failed_loose"], subset))
        end
        println(io, "\nClass counts: class | definition | rows | failed | fine | silent (default combined checks)")
        for row in classes
            println(io, join([row[key] for key in ("cause_class", "definition", "rows", "failed", "fine", "silent")], " | "))
        end
        for cause in ("pooled", sort(unique(row["cause_class"] for row in rows))...), definition in ("strict", "loose")
            subset = filter(row -> row["cause_class"] == cause && row["definition"] == definition, detection)
            println(io, "\n", cause, "; ", definition, " detection: check | flagged-failed | flagged-fine | missed | quiet-fine | unavailable-failed | unavailable-fine")
            for row in subset
                println(io, join([row[k] for k in ("check", "flagged_and_failed", "flagged_and_fine", "missed_failures", "quiet_and_fine", "unavailable_failed", "unavailable_fine")], " | "))
            end
        end
        println(io, "\nSilent failures (full identities and outcomes in verification_silent.csv): ", length(silent))
        for row in silent
            println(io, join([get(row, k, "") for k in ("cause_class", "definition", "check_set", "study", "cable", "fit", "scenario", "receiving_core_peak_error_percent", "receiving_energized_sheath_peak_error_percent", "receiving_sheath_scaled_peak_error_percent", "DC_endpoint_error_percent", "divergence")], " | "))
        end
    end
    write(joinpath(output, "verification_sources.txt"), join(folders, '\n') * "\n")
    println("Wrote ", length(rows), " fitted configuration/scenario rows to ", output)
    return (; rows, detection, compact, silent, coverage)
end

function validate(folders::AbstractVector{<:AbstractString}, ::typeof(verification_tables), output)
    all(isdir, folders) || throw(ArgumentError("Every input must be a run folder"))
    any(folder -> output == folder || startswith(output, folder * "/"), folders) &&
        throw(ArgumentError("Output must be outside the input run folders"))
    return folders
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) >= 2 || error("Usage: julia --project=test benchmarks/verification_tables.jl <run folders>... <output folder>")
    folders = unique(abspath.(ARGS[1:end-1]))
    output = abspath(last(ARGS))
    verification_tables(folders, output)
end
