# Read the newest complete ranking.txt for each of the eight cable scripts under
# `.linecablemodels/benchmarks/<script>/<yyyymmdd-HHMMSS>/` in this worktree.
# Write twelve LaTeX tables to the folder given as the first argument. Without an argument,
# the tables go to `.linecablemodels/benchmarks/tables/`. The script prints the folders it
# reads and the files it writes. No cable models are computed.
# Run: `julia benchmarks/sensitivity_tables.jl [output folder]`
# The tables use longtable, booktabs and siunitx, with three significant digits.

function read_ranking(path)
    rows = NamedTuple[]
    ratios = Dict{String,Float64}()
    counts = Dict{String,Int}()
    summaries = Set{String}()
    metric_pattern = r"^(one_at_a_time|temperature|earth|representation|uncertainty) ([^ =]+)(?:=([^ ]+))? ([ZY]) RMS(?:_e)? \(all,dc,harmonic,narrow\)=\(([^)]+)\) max(?:_e)?=([^ ]+) frequency_Hz=([^ ]+)$"
    for (number, line) in enumerate(eachline(path))
        metric = match(metric_pattern, line)
        if metric !== nothing
            section, name, level, quantity, bands, peak, frequency = metric.captures
            rms = Tuple(parse.(Float64, strip.(split(bands, ','))))
            length(rms) == 4 || throw(ArgumentError("$path:$number: expected four RMS bands"))
            push!(rows, (; section, name, level = isnothing(level) ? nothing : parse(Float64, level),
                quantity, rms, peak = parse(Float64, peak), frequency = parse(Float64, frequency)))
        elseif startswith(line, "uncertainty LEP/MC")
            ratio = match(r"^uncertainty LEP/MC RMS_s all ([ZY])=([^ ]+)$", line)
            isnothing(ratio) && throw(ArgumentError("$path:$number: malformed LEP/MC ratio"))
            haskey(ratios, ratio[1]) && throw(ArgumentError("$path:$number: duplicate LEP/MC ratio"))
            ratios[ratio[1]] = parse(Float64, ratio[2])
        elseif occursin(" RMS", line)
            throw(ArgumentError("$path:$number: malformed ranking metric"))
        else
            summary = match(r"^(one_at_a_time|temperature|earth): \d+ swept inputs, (\d+) computed points,", line)
            variants = match(r"^representation: (\d+) variants,", line)
            uncertainty = match(r"^uncertainty (MonteCarlo|LinearError):", line)
            if summary !== nothing
                counts[summary[1]] = parse(Int, summary[2])
            elseif variants !== nothing
                counts["representation"] = parse(Int, variants[1])
            elseif uncertainty !== nothing
                push!(summaries, uncertainty[1])
            end
        end
    end
    incomplete = String[]
    for section in ("one_at_a_time", "temperature", "earth", "representation")
        entries = filter(row -> row.section == section, rows)
        z = Set((row.name, row.level) for row in entries if row.quantity == "Z")
        y = Set((row.name, row.level) for row in entries if row.quantity == "Y")
        expected = get(counts, section, 0)
        if expected == 0 || length(entries) != 2expected || length(z) != expected || z != y
            push!(incomplete, section)
        end
    end
    entries = filter(row -> row.section == "uncertainty", rows)
    expected = Set((name, quantity) for name in ("MonteCarlo_s", "LinearError_s", "MonteCarlo_m")
        for quantity in ("Z", "Y"))
    if length(entries) != 6 || Set((row.name, row.quantity) for row in entries) != expected ||
            Set(keys(ratios)) != Set(("Z", "Y")) || summaries != Set(("MonteCarlo", "LinearError"))
        push!(incomplete, "uncertainty")
    end
    isempty(incomplete) || throw(ArgumentError("$path: missing or incomplete sections: $(join(incomplete, ", "))"))
    return (; rows, ratios)
end

tex_text(value) = replace(string(value), '\\' => raw"\textbackslash{}", '_' => raw"\_\allowbreak{}",
    '%' => raw"\%", '&' => raw"\&", '#' => raw"\#", '$' => raw"\$", '{' => raw"\{", '}' => raw"\}")
tex_number(value) = string(raw"\num[round-mode=figures,round-precision=3]{", value, "}")

function write_table(body, filename, caption, columns, header, sources)
    path = joinpath(output, filename)
    open(path, "w") do io
        for source in sources
            println(io, "% Source: ", source)
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

root = dirname(dirname(abspath(PROGRAM_FILE)))
output = isempty(ARGS) ? joinpath(root, ".linecablemodels", "benchmarks", "tables") : abspath(ARGS[1])
cases = NamedTuple[]
for (script, label) in (
        ("cable_320kv_armoured_dc_bipole", "320 kV DC, armoured"),
        ("cable_320kv_no_armour_dc_bipole", "320 kV DC, unarmoured"),
        ("cable_380kv_armoured_ac_flat", "380 kV AC, armoured"),
        ("cable_380kv_no_armour_ac_flat", "380 kV AC, unarmoured"),
        ("cable_525kv_land_no_armour_dc_bipole", "525 kV DC, land, unarmoured"),
        ("cable_525kv_land_no_armour_ac_flat", "525 kV AC, land, unarmoured"),
        ("cable_525kv_subsea_armoured_dc_bipole", "525 kV DC, subsea, armoured"),
        ("cable_525kv_subsea_armoured_ac_flat", "525 kV AC, subsea, armoured"))
    directory = joinpath(root, ".linecablemodels", "benchmarks", script)
    isdir(directory) || error("$label ($script): no output folder; missing all five ranking sections")
    folders = sort(filter(name -> occursin(r"^\d{8}-\d{6}$", name) &&
        isdir(joinpath(directory, name)), readdir(directory)); rev = true)
    selected = false
    rejected = String[]
    for folder in folders
        path = joinpath(directory, folder, "ranking.txt")
        if !isfile(path)
            push!(rejected, "$folder: ranking.txt missing (all five sections)")
            continue
        end
        ranking = try
            read_ranking(path)
        catch exception
            exception isa ArgumentError || rethrow()
            push!(rejected, sprint(showerror, exception))
            continue
        end
        push!(cases, (; script, label, source = path, ranking.rows, ranking.ratios))
        println("Used ", label, " (", script, "): ", dirname(path))
        selected = true
        break
    end
    for reason in rejected
        println("Skipped ", label, ": ", reason)
    end
    selected || error("$label ($script): no complete ranking.txt with one_at_a_time, temperature, earth, representation, uncertainty. $(join(rejected, "; "))")
end

mkpath(output)
sources = [case.source for case in cases]
labels = join([case.label for case in cases], ", ")

for case in cases
    write_table("drivers_$(case.script).tex",
            "One-at-a-time sensitivity, $(case.label): ten largest all-band RMS deviations for each of Y and Z. Levels use the case inputs' base units; deviations are fractions.",
            raw"@{}lp{0.24\linewidth}rrrrrrr@{}",
            raw"Q & Input & Level & RMS all & RMS dc & RMS harmonic & RMS narrow & Max $e$ & $f$ (Hz) ",
            [case.source]) do io
        println(io, raw"\multicolumn{9}{l}{RMS and maximum deviations: fractions}\\\\")
        for quantity in ("Y", "Z")
            rows = sort(filter(row -> row.section == "one_at_a_time" && row.quantity == quantity,
                case.rows); by = row -> row.rms[1], rev = true)
            length(rows) >= 10 || error("$(case.label): fewer than ten one_at_a_time entries for $quantity")
            for row in first(rows, 10)
                values = (row.level, row.rms..., row.peak, row.frequency)
                println(io, quantity, " & ", tex_text(row.name), " & ",
                    join(tex_number.(values), " & "), raw"\\\\")
            end
            quantity == "Y" && println(io, raw"\midrule")
        end
    end
end

write_table("temperature_earth.tex",
        "Temperature and earth resistivity sensitivity: $labels. All computed levels are reproduced; deviations are fractions.",
        raw"@{}llrlrrr@{}",
        raw"Case & Input (unit) & Level & Q & RMS all (fraction) & Max $e$ (fraction) & $f$ (Hz)",
        sources) do io
    for case in cases, section in ("temperature", "earth")
        rows = sort(filter(row -> row.section == section, case.rows);
            by = row -> (row.level, row.quantity))
        input = section == "temperature" ? raw"Temperature (${}^\circ$C)" : raw"Earth ($\Omega\,\mathrm{m}$)"
        for row in rows
            println(io, case.label, " & ", input, " & ", tex_number(row.level), " & ", row.quantity,
                " & ", join(tex_number.((row.rms[1], row.peak, row.frequency)), " & "), raw"\\\\")
        end
    end
end

write_table("representation.tex",
        "Representation sensitivity: $labels. Every computed variant is included; deviations are fractions.",
        raw"@{}lllrrr@{}",
        raw"Case & Variant & Q & RMS all (fraction) & Max $e$ (fraction) & $f$ (Hz)",
        sources) do io
    for case in cases
        rows = sort(filter(row -> row.section == "representation", case.rows);
            by = row -> (row.name, row.quantity))
        for row in rows
            println(io, case.label, " & ", tex_text(row.name), " & ", row.quantity, " & ",
                join(tex_number.((row.rms[1], row.peak, row.frequency)), " & "), raw"\\\\")
        end
    end
end

write_table("uncertainty.tex",
        "Joint uncertainty: $labels. MC denotes Monte Carlo and LEP denotes linear error propagation. Spread and mean shift are fractions; LEP/MC is the all-band RMS spread ratio.",
        raw"@{}lllrrrr@{}",
        raw"Case & Q & Statistic & RMS all (fraction) & Max (fraction) & $f$ (Hz) & LEP/MC",
        sources) do io
    for case in cases, quantity in ("Z", "Y")
        for (name, label) in (("MonteCarlo_s", "MC spread"), ("LinearError_s", "LEP spread"),
                ("MonteCarlo_m", "MC mean shift"))
            row = only(filter(row -> row.section == "uncertainty" && row.name == name &&
                row.quantity == quantity, case.rows))
            ratio = name == "LinearError_s" ? tex_number(case.ratios[quantity]) : "{---}"
            println(io, case.label, " & ", quantity, " & ", label, " & ",
                join(tex_number.((row.rms[1], row.peak, row.frequency)), " & "),
                " & ", ratio, raw"\\\\")
        end
    end
end

write_table("time_domain_selection.tex",
        "Time domain case selection: 320 kV DC, armoured, 320 kV DC, unarmoured and 525 kV DC, subsea, armoured.",
        raw"@{}lp{0.27\linewidth}p{0.54\linewidth}@{}", "Case & Script & Reason",
        [case.source for case in cases if case.label in ("320 kV DC, armoured", "320 kV DC, unarmoured", "525 kV DC, subsea, armoured")]) do io
    for (label, reason) in (
            ("320 kV DC, armoured", "Armoured system cable: stranded core and wire armour with high permeability steel."),
            ("320 kV DC, unarmoured", "Unarmoured system cable: paired with 320 kV DC, armoured to assess the effect of armour."),
            ("525 kV DC, subsea, armoured", "Subsea cable: different solid core, armour geometry and lower permeability steel."))
        case = only(filter(case -> case.label == label, cases))
        println(io, label, " & ", tex_text(case.script), " & ", reason, raw"\\\\")
    end
end
