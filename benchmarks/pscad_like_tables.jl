# Derive earth-model and archived PSCAD comparisons from saved observations.
using LineCableModels
using LineCableModels.Engine: compare
using LineCableModels.Commons: observation_product
using LinearAlgebra: diag
using TOML

length(ARGS) >= 2 || error("Usage: julia --project=test benchmarks/pscad_like_tables.jl <run folders>... <output folder>")
output = abspath(last(ARGS))
source_file = joinpath(dirname(dirname(@__DIR__)), "LineCableModels-benchmarks",
    "local", "report", "tables", "time_domain_sources.txt")
archived_sources = unique(strip.(readlines(source_file)))

function csv_row(io, values)
    println(io, join((v === nothing ? "" : "\"" * replace(string(v), "\"" => "\"\"") * "\"" for v in values), ','))
end

function write_csv(path, rows, columns)
    open(path, "w") do io
        csv_row(io, columns)
        for row in rows
            csv_row(io, (getproperty(row, key) for key in columns))
        end
    end
    println("Wrote ", path, ": ", length(rows), " rows")
end

tex_text(value) = replace(string(value), '_' => raw"\_\allowbreak{}", '%' => raw"\%", '&' => raw"\&", '#' => raw"\#")
tex_number(::Nothing) = "{---}"
tex_number(value::Number) = string(raw"\num[round-mode=figures,round-precision=6]{", value, "}")
tex_number(value) = tex_text(value)

function write_tex(path, rows, columns, header, caption)
    open(path, "w") do io
        println(io, raw"\begin{longtable}{", "@{}" * repeat("l", length(columns)) * "@{}", "}")
        println(io, raw"\caption{", caption, raw"}\\\\")
        println(io, raw"\toprule", '\n', join(header, " & "), raw"\\\\", '\n', raw"\midrule\endfirsthead")
        println(io, raw"\toprule", '\n', join(header, " & "), raw"\\\\", '\n', raw"\midrule\endhead")
        println(io, raw"\bottomrule\endfoot")
        for row in rows
            println(io, join((tex_number(getproperty(row, key)) for key in columns), " & "), raw"\\\\")
        end
        println(io, raw"\end{longtable}")
    end
end

function response(point)
    q = only(q for q in point.quantities if q.request === voltage)
    q.unit == LineCableModels.Units.units(:base, :volt) || error("Expected archived volts")
    q.coordinates.time_unit == LineCableModels.Units.units(:base, :second) || error("Expected archived seconds")
    return TransientResponse(q.coordinates.times, observe(point, voltage), observe(point, current),
        q.coordinates.column_labels, nothing)
end

function differences(reference, result, terminals)
    deviation = compare(reference, result, voltage)
    left, right = peak(reference, voltage), peak(result, voltage)
    return [(; terminal, rms_kV = deviation.absolute[2, i] / 1e3,
        scaled_rms = deviation.absolute[2, i] / left.value[2, i],
        reference_peak_kV = left.value[2, i] / 1e3,
        result_peak_kV = right.value[2, i] / 1e3,
        peak_error_kV = (right.value[2, i] - left.value[2, i]) / 1e3,
        peak_error_fraction = (right.value[2, i] - left.value[2, i]) / left.value[2, i],
        maximum_difference_kV = maximum(abs, result.voltage[2, i, :] - reference.voltage[2, i, :]) / 1e3)
        for (i, terminal) in enumerate(reference.coordinates) if terminal in terminals]
end

for source in ARGS[1:end-1]
    folder = abspath(source)
    script, timestamp = basename(dirname(folder)), basename(folder)
    destination = mkpath(joinpath(output, script, timestamp))
    sources = [folder, source_file]
    unified_points = import_data(:observed, joinpath(folder, "references.json"))
    pscad_like_points = import_data(:observed, joinpath(folder, "pscad_like_references.json"))
    references = Dict()
    reference_rows = NamedTuple[]
    for point in pscad_like_points
        scenario = point.gridpoint.study.scenario
        other = only(p for p in unified_points if p.gridpoint.study.scenario === scenario)
        unified, pscad_like = response(other), response(point)
        references[(scenario, :unified)] = unified
        references[(scenario, :pscad_like)] = pscad_like
        terminals = filter(t -> endswith(t, ":1:core") || endswith(t, ":sheath"), unified.coordinates)
        for row in differences(unified, pscad_like, terminals)
            push!(reference_rows, merge((; scenario), row))
        end
    end
    length(references) == 6 || error("Expected three pairs of reference scenarios in $folder")
    write_csv(joinpath(destination, "reference_differences.csv"), reference_rows,
        (:scenario, :terminal, :rms_kV, :scaled_rms, :reference_peak_kV, :result_peak_kV,
            :peak_error_kV, :peak_error_fraction, :maximum_difference_kV))
    write_tex(joinpath(destination, "reference_differences.tex"), reference_rows,
        (:scenario, :terminal, :scaled_rms, :peak_error_kV),
        ("Scenario", "Terminal", "RMS / Unified peak", "Peak difference (kV)"),
        "Pollaczek impedance with ideal earth admittance against Unified. Receiving voltage RMS uses the Unified peak. Signed peak differences subtract the Unified absolute peak.")

    earth_rows = NamedTuple[]
    for point in import_data(:observed, joinpath(folder, "pscad_like_earth.json"))
        study = point.gridpoint.study
        core = study.status == "computable" ? only(p for p in study.receiving_peaks if endswith(p.terminal, ":1:core")) : nothing
        push!(earth_rows, (; scenario = study.scenario, earth_rho_ohm_m = study.earth_rho,
            status = study.status, passivity = study.passivity.status,
            highest_frequency_Hz = study.nlt.N / study.nlt.T - 1 / (2study.nlt.T),
            receiving_peak_kV = core === nothing ? nothing : core.peak_kV,
            receiving_peak_time_us = core === nothing ? nothing : core.peak_time_us,
            minimum_Z_eigenvalue_ohm_per_m = study.passivity.minimum_Z_eigenvalue_ohm_per_m,
            minimum_Y_eigenvalue_S_per_m = study.passivity.minimum_Y_eigenvalue_S_per_m))
    end
    write_csv(joinpath(destination, "earth_resistivity.csv"), earth_rows,
        (:scenario, :earth_rho_ohm_m, :status, :passivity, :highest_frequency_Hz,
            :receiving_peak_kV, :receiving_peak_time_us, :minimum_Z_eigenvalue_ohm_per_m, :minimum_Y_eigenvalue_S_per_m))

    pscad_rows = NamedTuple[]
    variants = filter(archived_sources) do path
        basename(dirname(dirname(path))) == script && occursin(r"^pscad_(step_single_point|impulse)_\d+$", basename(path))
    end
    # Accept study folders in the source list and collect their archived variants.
    for path in archived_sources
        basename(dirname(path)) == script || continue
        append!(variants, filter(p -> occursin(r"^pscad_(step_single_point|impulse)_\d+$", basename(p)),
            readdir(path; join = true)))
    end
    isempty(variants) && error("No archived PSCAD runs for $script in $source_file")
    for directory in sort!(unique(variants))
        status = TOML.parsefile(joinpath(directory, "status.toml"))
        scenario = Symbol(status["scenario"])
        fit = status["fit"]
        push!(sources, directory)
        identity = (; scenario, run = basename(directory), model = status["model"],
            order = get(fit, "YMaxP", nothing), DCenab = get(fit, "DCenab", nothing),
            CPASS = get(fit, "CPASS", nothing), status = status["status"], source = directory)
        completed = status["status"] == "completed"
        native = completed ? response(import_data(:observed, joinpath(directory, "response.json"))) : nothing
        for model in (:unified, :pscad_like)
            reference = references[(scenario, model)]
            core1 = only(t for t in reference.coordinates if endswith(t, ":1:core"))
            metrics = completed ? only(differences(reference, native, (core1,))) :
                (; terminal = core1, rms_kV = nothing, scaled_rms = nothing, reference_peak_kV = nothing,
                    result_peak_kV = nothing, peak_error_kV = nothing, peak_error_fraction = nothing, maximum_difference_kV = nothing)
            push!(pscad_rows, merge(identity, (; reference = model), metrics))
        end
    end
    write_csv(joinpath(destination, "pscad_transients.csv"), pscad_rows,
        (:scenario, :run, :model, :order, :DCenab, :CPASS, :status, :reference, :terminal,
            :rms_kV, :scaled_rms, :reference_peak_kV, :result_peak_kV, :peak_error_kV,
            :peak_error_fraction, :maximum_difference_kV, :source))
    write_tex(joinpath(destination, "pscad_transients.tex"), pscad_rows,
        (:scenario, :run, :reference, :status, :scaled_rms, :peak_error_kV),
        ("Scenario", "Run", "Reference", "Status", "RMS / reference peak", "Peak error (kV)"),
        "Archived PSCAD receiving core 1 voltage against both NLT references. Deviations use each reference peak. Signed errors subtract each reference absolute peak. Failed runs retain empty metrics.")

    native_path = joinpath(folder, "pscad_native_line_parameters.json")
    if isfile(native_path)
        owned = import_data(:observed, joinpath(folder, "pscad_like_line_parameters.json"))
        native = import_data(:observed, native_path)
        frequency_rows, entry_rows = NamedTuple[], NamedTuple[]
        for (quantity, real_selector, imaginary_selector) in ((:Z, R, X), (:Y, G, B))
            left = observation_product(owned, real_selector)
            right = observation_product(native, real_selector)
            left.unit == right.unit && left.basis == right.basis || error("Line-parameter units or bases differ")
            left.coordinates.frequencies == right.coordinates.frequencies || error("Line-parameter frequency grids differ")
            labels = left.coordinates.labels
            labels == right.coordinates.labels || error("Line-parameter terminal orders differ")
            a = complex.(observe(owned, real_selector), observe(owned, imaginary_selector))
            b = complex.(observe(native, real_selector), observe(native, imaginary_selector))
            cables = [split(label, ':')[2] for label in labels]
            outer = [findlast(==(cable), cables) for cable in unique(cables)]
            groups = ((:all, [(i, j) for i in eachindex(labels) for j in eachindex(labels)]),
                (:self, [(i, i) for i in eachindex(labels)]),
                (:mutual, [(i, j) for i in eachindex(labels) for j in eachindex(labels) if i != j]),
                (:outer_self, [(i, i) for i in outer]),
                (:earth_return_mutual, [(i, j) for i in outer for j in outer if cables[i] != cables[j]]))
            for (k, frequency_Hz) in enumerate(left.coordinates.frequencies)
                scale = maximum(abs, diag(a[:, :, k]))
                ratio(value) = iszero(scale) ? (iszero(value) ? 0.0 : Inf) : value / scale
                for (group, entries) in groups
                    isempty(entries) && continue
                    push!(frequency_rows, (; quantity, frequency_Hz, group,
                        scaled_deviation = ratio(maximum(abs(b[i, j, k] - a[i, j, k]) for (i, j) in entries))))
                end
                for i in eachindex(labels), j in eachindex(labels)
                    push!(entry_rows, (; quantity, frequency_Hz, row = labels[i], column = labels[j],
                        group = i == j ? :self : :mutual, same_cable = cables[i] == cables[j],
                        unit = quantity === :Z ? "ohm/m" : "S/m",
                        owned_real = real(a[i, j, k]), owned_imaginary = imag(a[i, j, k]),
                        native_real = real(b[i, j, k]), native_imaginary = imag(b[i, j, k]),
                        scaled_deviation = ratio(abs(b[i, j, k] - a[i, j, k]))))
                end
            end
        end
        write_csv(joinpath(destination, "line_parameters.csv"), frequency_rows,
            (:quantity, :frequency_Hz, :group, :scaled_deviation))
        write_csv(joinpath(destination, "line_parameter_entries.csv"), entry_rows,
            (:quantity, :frequency_Hz, :row, :column, :group, :same_cable, :unit,
                :owned_real, :owned_imaginary, :native_real, :native_imaginary, :scaled_deviation))
    end
    write(joinpath(destination, "sources.txt"), join(unique(sources), '\n') * "\n")
    write(joinpath(destination, "definitions.txt"),
        "RMS deviations divide receiving voltage RMS differences by each reference absolute peak.\n" *
        "Peak errors subtract reference absolute peaks. Maximum differences use aligned waveform samples.\n" *
        "Line-parameter scaled deviations divide entry differences by the largest absolute owned diagonal at each frequency.\n" *
        "Self and mutual denote phase diagonal and off-diagonal entries. Outer self includes internal impedance.\n" *
        "Cross-cable outer Z mutual entries isolate earth return. Native isolated earth self coefficients are unavailable.\n" *
        "Native defaults also change conductor equations and dielectric export representation. Total differences include these effects.\n")
end
