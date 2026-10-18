# Read the newest complete time-domain, PSCAD and model-settings archives.
# Run: `julia --project=test benchmarks/time_domain_tables.jl [output folder]`
# Tables use longtable, booktabs and siunitx, with three significant digits.
using LineCableModels
import LineCableModels: validate
using TOML
using XLSX

root = dirname(@__DIR__)
output = isempty(ARGS) ? joinpath(root, ".linecablemodels", "benchmarks", "tables") : abspath(ARGS[1])
length(ARGS) <= 1 || error("Expected at most one output folder")
scenarios = (:step_single_point, :step_both_ends, :impulse)
sources = String[]

tex_text(value) = replace(string(value), '\\' => raw"\textbackslash{}", '_' => raw"\_\allowbreak{}",
    '%' => raw"\%", '&' => raw"\&", '#' => raw"\#", '$' => raw"\$", '{' => raw"\{", '}' => raw"\}")
tex_number(::Nothing) = "{---}"
tex_number(value) = isnan(value) ? "{NaN}" : isinf(value) ?
    (value > 0 ? raw"$\infty$" : raw"$-\infty$") :
    string(raw"\num[round-mode=figures,round-precision=3]{", value, "}")
tex_range(values) = isempty(values) ? "{---}" :
    minimum(values) == maximum(values) ? tex_number(first(values)) :
    string(tex_number(minimum(values)), "--", tex_number(maximum(values)))
tex_flag(value) = value === nothing ? "unchecked" : value ? "yes" : "no"
dc_method(study) = get(study.options.dc, :method, :none)
dc_label(method) = method === :none ? "no correction" :
    method === :add_pole ? "added pole" : "functional form"
family(study) = study.method === :ulm ? "ULM, $(dc_label(dc_method(study)))" :
    study.method === :fdcm ? "FDCM" : "J. Marti"

function read_time_domain(folder)
    validate(folder, read_time_domain, ("segment.json", "line_models.json", "references.json",
        "coarse_baseline.json", "sampling_comparison.json"))
    fits = [point.gridpoint.study for point in import_data(:observed, joinpath(folder, "line_models.json"))]
    length(fits) == 54 || throw(ArgumentError("expected 54 baseline and coarse fits"))
    for method in (:none, :add_pole, :functional_form)
        selected = filter(s -> s.line === :baseline && s.method === :ulm && dc_method(s) === method, fits)
        length(selected) == 16 || throw(ArgumentError("expected 16 ULM variants for $method"))
        for enforce in (false, true)
            count(s -> s.options.passivity == enforce, selected) == 8 ||
                throw(ArgumentError("expected eight ULM variants per passivity setting"))
        end
    end
    for method in (:fdcm, :marti)
        count(s -> s.line === :baseline && s.method === method, fits) == 1 ||
            throw(ArgumentError("expected one baseline $method fit"))
    end
    for (index, study) in enumerate(fits)
        validate(folder, read_time_domain, ["line_model_$(index)_$(q)_1_fitting_error.xlsx" for q in (:admittance, :propagation)])
        validate(folder, read_time_domain, ["$(scenario)_response_$(index).json" for scenario in scenarios
            if scenario === :step_single_point || study.line === :baseline])
    end
    references = import_data(:observed, joinpath(folder, "references.json"))
    length(references) == 5 || throw(ArgumentError("expected five NLT references"))
    sampling = [point.gridpoint for point in import_data(:observed, joinpath(folder, "sampling_comparison.json"))]
    length(sampling) == 4 || throw(ArgumentError("expected four sampling comparisons"))
    return (; fits, sampling)
end

function read_pscad(folder)
    validate(folder, read_pscad, ("pscad_station.toml", "references.json"))
    station = TOML.parsefile(joinpath(folder, "pscad_station.toml"))
    rows = NamedTuple[]
    for scenario in (:step_single_point, :impulse), index in 1:13
        directory = joinpath(folder, "pscad_$(scenario)_$(index)")
        validate(directory, read_pscad, ("status.toml",))
        status = TOML.parsefile(joinpath(directory, "status.toml"))
        status["scenario"] == string(scenario) || throw(ArgumentError("scenario mismatch in $directory"))
        status["status"] == "completed" && validate(directory, read_pscad, ("response.json",))
        push!(rows, (; directory, scenario, index, status))
    end
    return (; station, rows)
end

function read_settings(folder)
    validate(folder, read_settings, ("model_settings.txt", "line_models.json", "reference.json"))
    fits = [point.gridpoint.model_settings for point in import_data(:observed, joinpath(folder, "line_models.json"))]
    length(fits) == 4 && [s.method for s in fits] == [:marti, :marti, :ulm, :ulm] ||
        throw(ArgumentError("expected two J. Marti and two ULM fits"))
    for (index, study) in enumerate(fits), interp in (study.method === :marti ? (:single,) : (:single, :double))
        validate(folder, read_settings, ("step_single_point_response_$(index)_$(interp).json",))
    end
    return (; fits)
end

function validate(folder::AbstractString,
        ::Union{typeof(read_time_domain),typeof(read_pscad),typeof(read_settings)}, filenames)
    absent = filter(name -> !isfile(joinpath(folder, name)), filenames)
    isempty(absent) || throw(ArgumentError("missing $(join(absent, ", "))"))
    return folder
end

function newest_complete(reader, directory, kind)
    isdir(directory) || error("No result folder: $directory")
    folders = sort(filter(name -> occursin(r"^\d{8}-\d{6}$", name) &&
        isdir(joinpath(directory, name)), readdir(directory)); rev = true)
    rejected = String[]
    for name in folders
        folder = joinpath(directory, name)
        result = try
            reader(folder)
        catch exception
            exception isa InterruptException && rethrow()
            push!(rejected, "$name: $(sprint(showerror, exception))")
            continue
        end
        push!(sources, folder)
        println("Used ", kind, ": ", folder)
        flush(stdout)
        return merge((; folder), result)
    end
    error("No complete $kind run in $directory: $(join(rejected, "; "))")
end

function write_table(body, filename, caption, columns, header, folders)
    path = joinpath(output, filename)
    open(path, "w") do io
        for folder in folders
            println(io, "% Source: ", folder)
        end
        println(io, raw"\begingroup\small\setlength{\tabcolsep}{3pt}")
        println(io, raw"\begin{longtable}{", columns, "}")
        println(io, raw"\caption{", caption, raw"}\\\\")
        println(io, raw"\toprule", '\n', header, raw"\\\\", '\n', raw"\midrule\endfirsthead")
        println(io, raw"\toprule", '\n', header, raw"\\\\", '\n', raw"\midrule\endhead")
        println(io, raw"\midrule\endfoot", '\n', raw"\bottomrule\endlastfoot")
        body(io)
        println(io, raw"\end{longtable}", '\n', raw"\endgroup")
    end
    println("Wrote ", path)
end

table_row(io, values...) = println(io, join(values, " & "), raw"\\\\")
cable_row(io, label, columns) = println(io, raw"\multicolumn{", columns, raw"}{l}{\textbf{", label, raw"}}\\\\")

function fit_groups(fits)
    groups = [findall(s -> s.line === :baseline && s.method === :ulm && dc_method(s) === method, fits)
        for method in (:none, :add_pole, :functional_form)]
    append!(groups, [findall(s -> s.line === :baseline && s.method === method, fits) for method in (:fdcm, :marti)])
    return groups
end

function workbook_error(folder, index, quantity, expected)
    path = joinpath(folder, "line_model_$(index)_$(quantity)_1_fitting_error.xlsx")
    values = XLSX.openxlsx(path) do book
        table = book["values"][:]
        table[1, 2] == "value" || error("Unexpected fitting-error workbook header: $path")
        Float64.(table[2:end, 2])
    end
    values == vec(expected.rms) || error("Workbook and archive fitting errors differ: $path")
    return maximum(values)
end

function residue_ratio(study)
    perturbation = study.passivity.perturbation
    return perturbation === nothing ? nothing : maximum(perturbation.res_pole_ratio)
end

function dc_relative(study, quantity)
    study.dc.targets === nothing && return nothing
    target = getproperty(study.dc.targets, quantity)
    achieved = getproperty(study.dc.achieved, quantity)
    # The relative scale is the maximum magnitude of the target matrix entries.
    scale = maximum(abs, target)
    difference = maximum(abs, achieved - target)
    return iszero(scale) ? (iszero(difference) ? 0.0 : Inf) : difference / scale
end

cases = NamedTuple[]
for (script, label) in (
        ("cable_320kv_no_armour_dc_bipole", "320 kV DC, unarmoured, land, bipole, stranded core"),
        ("cable_320kv_armoured_dc_bipole", "320 kV DC, armoured, land, bipole, stranded core"),
        ("cable_525kv_subsea_armoured_dc_bipole", "525 kV DC, subsea, armoured, bipole, solid equivalent core"))
    directory = joinpath(root, ".linecablemodels", "benchmarks", script)
    time_domain = newest_complete(read_time_domain, directory, "time domain")
    pscad = newest_complete(read_pscad, directory, "PSCAD")
    settings = newest_complete(read_settings, directory, "model settings")
    push!(cases, (; label, time_domain, pscad, settings))
end

mkpath(output)
write_table("fitting_errors.tex",
        "Baseline fits: maximum matrix-entry band RMS errors (fractions); ULM ranges cover 16 variants per DC method. Passive counts use the eight enforced ULM variants; unchecked fits are marked. Residue/pole ratios are H maxima from archived passivity perturbations, with coverage shown; other ratios were not retained. Wall times are seconds.",
        "@{}lrrrrrr@{}", raw"Model & $Y_c$ RMS & $H$ RMS & Passive & Max $|R|/|p|$ & Coverage & Fit (s)",
        [c.time_domain.folder for c in cases]) do io
    for case in cases
        cable_row(io, case.label, 7)
        for indices in fit_groups(case.time_domain.fits)
            studies = case.time_domain.fits[indices]
            errors = map((:admittance, :propagation)) do q
                [workbook_error(case.time_domain.folder, index, q,
                    getproperty(case.time_domain.fits[index].fitting_error, q)) for index in indices]
            end
            enforced = filter(s -> s.method !== :ulm || s.options.passivity, studies)
            checked = filter(s -> s.passivity.checked, enforced)
            passive = isempty(checked) ? "unchecked" : "$(count(s -> s.passivity.passive === true, checked))/$(length(checked))"
            ratios = filter(!isnothing, residue_ratio.(studies))
            table_row(io, family(first(studies)), tex_range(errors[1]), tex_range(errors[2]),
                passive, isempty(ratios) ? "{---}" : tex_number(maximum(ratios)),
                "$(length(ratios))/$(length(studies))", tex_range([s.wall_seconds for s in studies]))
        end
    end
end

write_table("dc_behaviour.tex",
        raw"Baseline ULM DC errors: $\max|\mathrm{achieved}-\mathrm{target}|/\max|\mathrm{target}|$, maximised over eight variants per enforcement setting. Errors are fractions. Dashes denote absent targets (no correction) or an accepted DC-limit frequency omitted from the archives and available logs.",
        "@{}lrrrrr@{}", raw"DC method & \shortstack{$Y_c$\\unenforced} & \shortstack{$H$\\unenforced} & \shortstack{$Y_c$\\enforced} & \shortstack{$H$\\enforced} & DC limit (Hz)",
        [c.time_domain.folder for c in cases]) do io
    for case in cases
        cable_row(io, case.label, 6)
        for indices in first(fit_groups(case.time_domain.fits), 3)
            studies = case.time_domain.fits[indices]
            errors = [filter(!isnothing, [dc_relative(s, q) for s in studies if s.options.passivity == enforce])
                for enforce in (false, true) for q in (:admittance, :propagation)]
            table_row(io, dc_label(dc_method(first(studies))),
                [isempty(e) ? "{---}" : tex_number(maximum(e)) for e in errors]..., "{---}")
        end
    end
end

responses = Dict()
for case in cases
    run = case.time_domain
    for scenario in scenarios
        entries = Dict{Int,NamedTuple}()
        for (index, study) in enumerate(run.fits)
            study.line === :baseline || continue
            point = import_data(:observed, joinpath(run.folder, "$(scenario)_response_$(index).json"))
            evidence = point.gridpoint.study
            evidence.scenario === scenario && evidence.fit == index || error("Response identity mismatch")
            entries[index] = (; deviation = evidence.scaled_rms, growth = evidence.stability_flag)
        end
        responses[(run.folder, scenario)] = entries
        println("Read ", scenario, " responses: ", run.folder)
        flush(stdout)
    end
end

pscad_responses = Dict()
for case in cases, row in case.pscad.rows
    push!(sources, row.directory)
    println("Used PSCAD variant: ", row.directory)
    row.status["status"] == "completed" || continue
    point = import_data(:observed, joinpath(row.directory, "response.json"))
    evidence = point.gridpoint.pscad
    evidence.scaled_rms == row.status["scaled_rms"] || error("PSCAD status and response deviations differ")
    pscad_responses[row.directory] = evidence.scaled_rms
end

for scenario in scenarios
    folders = [c.time_domain.folder for c in cases]
    scenario !== :step_both_ends && append!(folders, [c.pscad.folder for c in cases])
    write_table("transient_deviation_$(scenario).tex",
            "$(tex_text(scenario)): receiving core 1 voltage RMS deviation divided by the NLT reference peak (fraction). ULM ranges and growth counts cover 16 baseline variants per DC method. PSCAD FD-Phase ranges cover four DCenab/CPASS settings per order; completed/total counts show availability. Growth counts use archived whole-response flags; PSCAD growth flags were not retained.",
            "@{}lrrr@{}", "Model & RMS / reference peak & Growth flags & Completed / total", folders) do io
        for case in cases
            cable_row(io, case.label, 4)
            entries = responses[(case.time_domain.folder, scenario)]
            for indices in fit_groups(case.time_domain.fits)
                table_row(io, family(case.time_domain.fits[first(indices)]),
                    tex_range([entries[i].deviation for i in indices]),
                    "$(count(i -> entries[i].growth, indices))/$(length(indices))",
                    "$(length(indices))/$(length(indices))")
            end
            scenario === :step_both_ends && continue
            for order in (10, 20, 30, nothing)
                rows = filter(case.pscad.rows) do row
                    row.scenario === scenario && (order === nothing ? row.status["model"] == "fd_mode" :
                        row.status["model"] == "fd_phase" && row.status["fit"]["YMaxP"] == order && row.status["fit"]["AMaxP"] == order)
                end
                values = [pscad_responses[row.directory] for row in rows if haskey(pscad_responses, row.directory)]
                label = order === nothing ? "PSCAD FD-Mode" : "PSCAD FD-Phase, order $order"
                table_row(io, label, tex_range(values), "{---}", "$(length(values))/$(length(rows))")
            end
        end
    end
end

stations = unique(["PSCAD $(c.pscad.station["version"]), compiler $(c.pscad.station["fortran_version"])" for c in cases])
write_table("pscad_runs.tex",
        "PSCAD runs; station identity: $(tex_text(join(stations, "; "))). Deviations are receiving core 1 voltage RMS divided by the NLT reference peak (fractions). Orders are YMaxP/AMaxP. Failures retain their status kind and first message line (up to 180 characters). Dashes mark default settings or absent deviations.",
        raw"@{}llrrrp{0.34\linewidth}r@{}", "Scenario & Model & Orders & DCenab & CPASS & Status / message head & Deviation",
        [c.pscad.folder for c in cases]) do io
    for case in cases
        cable_row(io, case.label, 7)
        for row in case.pscad.rows
            status = row.status
            fit = status["fit"]
            state = status["status"] == "completed" ? "completed" :
                "$(status["status"]): $(first(first(split(status["message"], '\n')), 180))"
            orders = status["model"] == "fd_phase" ? "$(fit["YMaxP"])/$(fit["AMaxP"])" : "{---}"
            label = status["model"] == "fd_phase" ? "FD-Phase" : "FD-Mode"
            table_row(io, tex_text(row.scenario), label, orders, tex_number(get(fit, "DCenab", nothing)),
                tex_number(get(fit, "CPASS", nothing)), tex_text(state), tex_number(get(pscad_responses, row.directory, nothing)))
        end
    end
end

write_table("sampling.tex",
        "Baseline and coarse sampling: stored sample counts include 50 Hz. Fitting errors are maximum matrix-entry band RMS errors; step deviations are receiving core 1 voltage RMS divided by the NLT reference peak. All errors are fractions; coarse ULM fits use orders 16/16, time rows and passivity enforcement.",
        "@{}llrrrr@{}", raw"Model & Grid & Samples & $Y_c$ RMS & $H$ RMS & Step deviation",
        [c.time_domain.folder for c in cases]) do io
    for case in cases
        cable_row(io, case.label, 6)
        for point in case.time_domain.sampling, grid in (:baseline, :coarse)
            comparison = point.sampling
            errors = getproperty(comparison.fitting_error, grid)
            table_row(io, family(point.study), string(grid), tex_number(getproperty(comparison.samples, grid)),
                tex_number(maximum(errors.admittance.rms)), tex_number(maximum(errors.propagation.rms)),
                tex_number(getproperty(comparison.scaled_rms, grid)))
        end
    end
end

write_table("model_settings.tex",
        "Single-point step model settings. J. Marti rows give requested and actual transformation frequencies, modal Z/Y off-diagonal ratios and the range of residual imaginary norms across modes. ULM rows give the H residue/pole limit, interpolation, whole-response growth flag, passivity and largest archived H residue/pole ratio. Deviations are receiving core 1 voltage RMS divided by the reference peak (fractions).",
        "@{}llrrrrrr@{}", raw"Model & \shortstack{Requested\\frequency (Hz)} & \shortstack{Used (Hz) /\\interpolation} & Deviation & \shortstack{Z ratio /\\growth} & \shortstack{Y ratio /\\passive} & \shortstack{Imaginary norm /\\max $|R|/|p|$} & Limit",
        [c.settings.folder for c in cases]) do io
    for case in cases
        cable_row(io, case.label, 8)
        for (index, study) in enumerate(case.settings.fits)
            for interp in (study.method === :marti ? (:single,) : (:single, :double))
                point = import_data(:observed, joinpath(case.settings.folder, "step_single_point_response_$(index)_$(interp).json"))
                transient = point.gridpoint.transient
                if study.method === :marti
                    requested = study.transformation_requested === nothing ? "default" : tex_number(study.transformation_requested)
                    table_row(io, "J. Marti", requested, tex_number(study.transformation_frequency),
                        tex_number(transient.scaled_rms), tex_number(study.modal.offdiagonal_ratio[1]),
                        tex_number(study.modal.offdiagonal_ratio[2]), tex_range(study.modal.transformation_imaginary_norm), "{---}")
                else
                    table_row(io, "Coarse ULM", "residue/pole", string(interp), tex_number(transient.scaled_rms),
                        tex_flag(transient.stability_flag), tex_flag(study.passivity.passive),
                        tex_number(study.largest_residue_pole_ratio), tex_number(study.max_ratio))
                end
            end
        end
    end
end

source_file = joinpath(output, "time_domain_sources.txt")
write(source_file, join(unique(sources), '\n') * "\n")
println("Wrote ", source_file)
