# Derive transient sensitivity tables from saved response observations.
# Tables use longtable, booktabs and siunitx.
using LineCableModels

length(ARGS) >= 2 || error("Usage: julia --project=test benchmarks/transient_sensitivity_tables.jl <run folders>... <output folder>")
output = abspath(last(ARGS))

function transient_csv_row(io, values)
    println(io, join((value === nothing ? "" :
        "\"" * replace(string(value), "\"" => "\"\"") * "\"" for value in values), ','))
end

transient_tex_text(value) = replace(string(value), '_' => raw"\_\allowbreak{}",
    '%' => raw"\%", '&' => raw"\&", '#' => raw"\#")
transient_tex_number(::Nothing) = "{---}"
transient_tex_number(value) = isnan(value) ? "{NaN}" :
    string(raw"\num[round-mode=figures,round-precision=6]{", value, "}")

function transient_values(point, terminals)
    computable = point.status == "computable"
    values = [
        (observable = :receiving_peak, line_end = :receiving, terminal = "cable:1:core",
            unit = "kV", value = computable ? point.receiving_peak_kV : nothing),
        (observable = :receiving_peak_time, line_end = :receiving, terminal = "cable:1:core",
            unit = "us", value = computable ? point.receiving_peak_time_us : nothing),
        (observable = :arrival_time, line_end = :receiving, terminal = "cable:1:core",
            unit = "us", value = computable ? point.arrival_time_us : nothing),
        (observable = :sending_current_peak, line_end = :sending, terminal = "cable:1:core",
            unit = "kA", value = computable ? point.sending_current_peak_kA : nothing),
        (observable = :receiving_scaled_rms, line_end = :receiving, terminal = "cable:1:core",
            unit = "1", value = computable ? point.receiving_scaled_rms : nothing),
    ]
    metal_peaks = computable ? point.metal_peaks :
        [(; line_end, terminal, peak_kV = nothing) for line_end in (:sending, :receiving)
            for terminal in terminals if endswith(terminal, ":sheath") || endswith(terminal, ":armor")]
    append!(values, [(observable = :metal_peak, line_end = row.line_end, terminal = row.terminal,
        unit = "kV", value = row.peak_kV) for row in metal_peaks])
    return values
end

for source in ARGS[1:end-1]
    folder = abspath(source)
    metadata = import_data(:observed, joinpath(folder, "transient_sensitivity_run.json")).gridpoint
    points = [point.gridpoint.sensitivity for point in
        import_data(:observed, joinpath(folder, "transient_sensitivity.json"))]
    all(point -> haskey(point, :status) && haskey(point, :passivity), points) ||
        error("Passivity records are absent in $folder; rerun transient_sensitivity")
    cable = basename(dirname(folder))
    destination = joinpath(output, cable, basename(folder))
    mkpath(destination)
    for scenario in (:step_single_point, :impulse)
        selected = filter(point -> point.scenario === scenario, points)
        nominal = only(filter(point -> point.study === :baseline, selected))
        baseline = transient_values(nominal, metadata.terminals)
        inputs = collect(metadata.swept)
        sort!(inputs; by = input -> maximum((abs(point.receiving_peak_kV - nominal.receiving_peak_kV)
            for point in selected if point.study === :one_at_a_time && point.input === input &&
                point.status == "computable" && nominal.status == "computable"); init = -Inf), rev = true)
        ordered = NamedTuple[nominal]
        for input in inputs
            append!(ordered, sort(filter(point -> point.study === :one_at_a_time && point.input === input &&
                point.level != metadata.nominal[input], selected); by = point -> point.level))
        end
        for study in (:temperature, :earth, :corrections)
            append!(ordered, filter(point -> point.study === study, selected))
        end
        rows = NamedTuple[]
        for point in ordered, (index, observable) in enumerate(transient_values(point, metadata.terminals))
            reference = baseline[index]
            (observable.observable, observable.line_end, observable.terminal, observable.unit) ==
                (reference.observable, reference.line_end, reference.terminal, reference.unit) ||
                error("Observable coordinates changed in $folder")
            computable = point.status == "computable" && nominal.status == "computable" &&
                observable.value !== nothing && reference.value !== nothing
            change = computable ? observable.value - reference.value : nothing
            percent = !computable || iszero(reference.value) ? nothing : 100change / abs(reference.value)
            rank = computable && point.study === :one_at_a_time ? findfirst(==(point.input), inputs) : nothing
            input_nominal = point.study in (:one_at_a_time, :temperature, :earth) ? metadata.nominal[point.input] : nothing
            input_change_percent = point.study === :one_at_a_time ? 100(point.level / input_nominal - 1) : nothing
            blocking = point.status == "not_computable" ? point.passivity : nominal.passivity
            push!(rows, (; study = point.study, input = point.input, level = point.level,
                input_nominal, input_change_percent, rank,
                observable.observable, observable.line_end, observable.terminal, observable.unit,
                status = computable ? "computable" : "not_computable",
                onset_frequency_Hz = blocking.onset_frequency_Hz, offending_quantity = blocking.offending_quantity,
                nominal = computable ? reference.value : nothing, value = computable ? observable.value : nothing,
                change, change_percent = percent))
        end
        prefix = "transient_$(scenario)"
        open(joinpath(destination, prefix * ".csv"), "w") do io
            transient_csv_row(io, keys(first(rows)))
            foreach(row -> transient_csv_row(io, values(row)), rows)
        end
        open(joinpath(destination, prefix * ".tex"), "w") do io
            println(io, "% Source: ", folder)
            for observable in baseline
                println(io, raw"\begingroup\small\setlength{\tabcolsep}{3pt}")
                println(io, raw"\begin{longtable}{lp{0.21\linewidth}lrrrrr}")
                caption = "$cable, $scenario: $(observable.observable), $(observable.line_end), $(observable.terminal)"
                unit = observable.unit == "us" ? raw"$\mu$s" : observable.unit
                println(io, raw"\caption{", transient_tex_text(caption), " [", unit, raw"]}\\\\")
                header = raw"Study & Input & Level & Rank & Nominal & Value & Change & Change (\%)\\\\"
                println(io, raw"\toprule", '\n', header, '\n', raw"\midrule\endfirsthead")
                println(io, raw"\toprule", '\n', header, '\n', raw"\midrule\endhead")
                for row in rows
                    (row.observable, row.line_end, row.terminal) ==
                        (observable.observable, observable.line_end, observable.terminal) || continue
                    level = row.study === :one_at_a_time ?
                        transient_tex_number(row.input_change_percent) * raw"\%" : transient_tex_text(row.level)
                    fields = [transient_tex_text(row.study), transient_tex_text(row.input),
                        level, row.rank === nothing ? "---" : string(row.rank)]
                    if row.status == "not_computable"
                        message = "not computable (" * transient_tex_text(row.offending_quantity) *
                            " onset " * transient_tex_number(row.onset_frequency_Hz) * " Hz)"
                        push!(fields, raw"\multicolumn{4}{l}{" * message * "}")
                    else
                        append!(fields, transient_tex_number.((row.nominal, row.value, row.change, row.change_percent)))
                    end
                    println(io, join(fields, " & "), raw"\\\\")
                end
                println(io, raw"\bottomrule\end{longtable}\endgroup")
            end
        end
    end
    open(joinpath(destination, "transient_passivity.csv"), "w") do io
        rows = [(; point.scenario, point.study, point.input, point.level, point.status,
            passivity_status = point.passivity.status,
            lowest_frequency_Hz = 1 / (2point.nlt.T), highest_frequency_Hz = (2point.nlt.N - 1) / (2point.nlt.T),
            point.passivity.onset_frequency_Hz, point.passivity.onset_lower_Hz,
            point.passivity.offending_quantity, point.passivity.minimum_Z_eigenvalue_ohm_per_m,
            point.passivity.minimum_Z_frequency_Hz, point.passivity.minimum_Y_eigenvalue_S_per_m,
            point.passivity.minimum_Y_frequency_Hz, point.passivity.archive, point.passivity.point) for point in points]
        transient_csv_row(io, keys(first(rows)))
        foreach(row -> transient_csv_row(io, values(row)), rows)
    end
    open(joinpath(destination, "transient_corrections_status.csv"), "w") do io
        transient_csv_row(io, ("correction", "status"))
        foreach(row -> transient_csv_row(io, row), metadata.correction_status)
    end
    open(joinpath(destination, "transient_corrections_status.tex"), "w") do io
        println(io, raw"\begin{longtable}{ll}\toprule Correction & Status\\\midrule")
        for (name, status) in metadata.correction_status
            println(io, transient_tex_text(name), " & ", transient_tex_text(status), raw"\\\\")
        end
        println(io, raw"\bottomrule\end{longtable}")
    end
    println("Transient tables: ", destination)
end
