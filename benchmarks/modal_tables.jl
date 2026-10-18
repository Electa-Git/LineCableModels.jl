using LineCableModels
using LinearAlgebra: dot, norm
using LineCableModels.Engine.ModalAnalysis: hungarian_assignment!

function modal_csv_row(io, values)
    println(io, join((value === nothing || ismissing(value) ? "" :
        "\"" * replace(string(value), "\"" => "\"\"") * "\"" for value in values), ','))
end

function modal_class(column, terminal_order)
    cables = unique(entry.cable for entry in terminal_order)
    terminals = unique(entry.terminal for entry in terminal_order)
    currents = reshape(column, length(terminals), length(cables))
    amplitudes = sqrt.(vec(sum(abs2, currents; dims = 2)))
    amplitudes ./= maximum(amplitudes)
    active = [j for j in axes(currents, 2) if norm(currents[:, j]) >= 0.2norm(column)]
    common = all(real(dot(currents[:, i], currents[:, j])) >= 0
        for i in active for j in active if i < j)
    parity = common ? "common" : length(cables) == 2 ?
        "cable1-against-cable2" : "cable-against-cable"
    opposite(i, j) = all(abs(angle(currents[i, c] / currents[j, c])) >= 3pi / 4
        for c in active)
    core, sheath = findfirst(==(:core), terminals), findfirst(==(:sheath), terminals)
    armour = findfirst(==(:armor), terminals)
    net_fraction = norm(vec(sum(currents; dims = 1))) /
        sum(sqrt.(vec(sum(abs2, currents; dims = 2))))
    family = if armour !== nothing && min(amplitudes[sheath], amplitudes[armour]) >= 0.5 &&
            amplitudes[core] <= 0.2 && opposite(sheath, armour)
        "sheath-armour"
    elseif min(amplitudes[core], amplitudes[sheath]) >= 0.5 &&
            (armour === nothing || amplitudes[armour] <= 0.2) && opposite(core, sheath)
        "coaxial core-sheath"
    elseif net_fraction >= 0.5
        terminal = terminals[argmax(amplitudes)]
        (terminal === :armor ? "armour" : string(terminal)) * "-dominated return"
    else
        armour === nothing ? "mixed core-sheath" : "mixed core-sheath-armour"
    end
    return family * "; " * parity
end

function modal_match(result, baseline, frequency_index)
    k = only(findall(==(frequency_index), result.basis_indices))
    current, reference = result.currents[:, :, k], baseline.currents[:, :, k]
    current ./= sqrt.(sum(abs2, current; dims = 1))
    reference ./= sqrt.(sum(abs2, reference; dims = 1))
    overlaps = abs.(current' * reference)
    n = size(overlaps, 1)
    buffers = (u = zeros(n + 1), v = zeros(n + 1), matching = zeros(Int, n + 1),
        way = zeros(Int, n + 1), minimums = zeros(n + 1), used = falses(n + 1),
        assignment = zeros(Int, n))
    assignment = hungarian_assignment!(-overlaps, buffers)
    best = vec(maximum(overlaps; dims = 2))
    second = [sort(overlaps[i, :]; rev = true)[2] for i in 1:n]
    assigned = [overlaps[i, assignment[i]] for i in 1:n]
    ambiguous = (best .< 0.9) .| (best .- second .<= 0.05) .| (assigned .< 0.9)
    return (; assignment, assigned, best, second, ambiguous)
end

function modal_quantities(result, length)
    speed = observe(result.curves, velocity) ./ 1e6
    return (attenuation_dB_per_km = observe(result.curves, (gamma, real)) .* (1000 * 20 / log(10)),
        velocity_m_per_us = speed, abs_Zc_ohm = abs.(complex.(observe(result.curves, (Zc, real)), observe(result.curves, (Zc, imag)))), travel_time_us = length ./ speed)
end

function modal_relative_jump(before, after)
    change = abs(after - before)
    return iszero(before) ? (iszero(change) ? 0.0 : Inf) : change / abs(before)
end

function modal_write_continuity(tables, summary, result, baseline, quantities, reference,
        point, study, frequencies, sampled)
    fallback_indices = result.diagnostics.fallback_frequencies
    exceeded = false
    maximum_ratio = 0.0
    for k in fallback_indices, mode in axes(quantities.velocity_m_per_us, 1),
            quantity in (:velocity_m_per_us, :attenuation_dB_per_km)
        first_index, last_index = max(1, k - 2), min(length(frequencies), k + 2)
        baseline_max = maximum(modal_relative_jump(reference[quantity][mode, j],
            reference[quantity][mode, j + 1]) for j in first_index:last_index-1)
        for (left, right) in ((k - 1, k), (k, k + 1))
            available = 1 <= left < right <= length(frequencies)
            jump = available ? modal_relative_jump(quantities[quantity][mode, left],
                quantities[quantity][mode, right]) : nothing
            ratio = !available ? nothing : iszero(baseline_max) ?
                (iszero(jump) ? 0.0 : Inf) : jump / baseline_max
            above = available && jump > 3baseline_max
            exceeded |= above
            available && (maximum_ratio = max(maximum_ratio, ratio))
            modal_csv_row(tables.fallback, (study, point, mode, mode, k, frequencies[k],
                left, right, available ? frequencies[left] : nothing,
                available ? frequencies[right] : nothing, quantity, jump,
                first_index, last_index, baseline_max, ratio, above))
            if above
                println(summary, "fallback jump point=", point, " mode=", mode,
                    " index=", k, " edge=", left, "->", right, " quantity=", quantity,
                    " relative_jump=", jump, " baseline_max=", baseline_max, " ratio=", ratio)
            end
        end
    end
    println(summary, "fallback continuity study=", study, " point=", point,
        " fallback_indices=", fallback_indices, " exceeds_threefold=", exceeded,
        " maximum_ratio=", maximum_ratio)
    for k in sampled[1:3]
        matched = modal_match(result, baseline, k)
        differs = any(mode -> matched.assignment[mode] != mode, eachindex(matched.assignment))
        low_overlap = any(<(0.9), matched.assigned)
        for mode in eachindex(matched.assignment)
            modal_csv_row(tables.order, (study, point, frequencies[k], mode,
                matched.assignment[mode], matched.assigned[mode], matched.best[mode],
                matched.second[mode], matched.assignment[mode] != mode,
                matched.assigned[mode] < 0.9, matched.ambiguous[mode]))
        end
        println(summary, "mode order study=", study, " point=", point,
            " frequency_Hz=", frequencies[k], " assignment=", matched.assignment,
            " overlap=", matched.assigned, " differs_from_identity=", differs,
            " overlap_below_0_9=", low_overlap)
    end
    return nothing
end

function modal_skin_transitions(p)
    permeability = LineCableModels.Commons.vacuum_permeability(Float64)
    sheath_rho = p.sheath_rho * (1 + p.sheath_temperature_coefficient *
        (p.temperature - p.reference_temperature))
    sheath = sheath_rho / (pi * permeability * p.sheath_mu_r * p.sheath_thickness^2)
    armour_length = nothing
    armour = nothing
    if haskey(p, :armour_wires)
        armour_length = p.armour_wire_diameter
        if p.armour_as_tube
            core_radius = get(p, :stranded_core, false) ? 13p.strand_diameter / 2 : p.core_diameter / 2
            inner = core_radius + p.inner_semicon_thickness + p.insulation_thickness +
                p.outer_semicon_thickness + p.water_blocking_thickness +
                p.sheath_thickness + p.inner_sheath_thickness + p.bedding_thickness
            armour_length = sqrt(inner^2 + p.armour_wires * (p.armour_wire_diameter / 2)^2) - inner
        end
        rho = p.armour_rho * (1 + p.armour_temperature_coefficient *
            (p.temperature - p.reference_temperature))
        armour = rho / (pi * permeability * p.armour_mu_r * armour_length^2)
    end
    return (; sheath, armour, armour_length)
end

function modal_write_point(table, conductors, summary, result, baseline, p, point,
        terminal_order, frequencies, sampled, match_index;
        continuity_tables, study, nonlinearities = nothing)
    matched = modal_match(result, baseline, match_index)
    diagnostics = result.diagnostics
    if !isempty(diagnostics.missed_frequencies) || !isempty(diagnostics.fallback_frequencies)
        println(summary, "modal diagnostics point=", point, " missed=", diagnostics.missed_frequencies,
            " fallback=", diagnostics.fallback_frequencies)
    end
    quantities = modal_quantities(result, p.line_length)
    reference = modal_quantities(baseline, p.line_length)
    modal_write_continuity(continuity_tables, summary, result, baseline, quantities, reference,
        point, study, frequencies, sampled)
    skin = modal_skin_transitions(p)
    for mode in eachindex(matched.assignment)
        if matched.ambiguous[mode]
            println(summary, "ambiguous point=", point, " mode=", mode,
                " baseline_mode=", matched.assignment[mode], " overlap=", matched.assigned[mode],
                " best=", matched.best[mode], " second=", matched.second[mode])
        end
        for (j, k) in enumerate(sampled), quantity in keys(quantities)
            value = quantities[quantity][mode, k]
            base = reference[quantity][matched.assignment[mode], k]
            linearity = nonlinearities === nothing ? nothing :
                nonlinearities[quantity][matched.assignment[mode], findfirst(==(k), sampled)]
            modal_csv_row(table, (point, mode, matched.assignment[mode],
                matched.assigned[mode], matched.best[mode], matched.second[mode], matched.ambiguous[mode],
                modal_class(result.currents[:, mode, j], terminal_order), frequencies[k], quantity,
                value, base, value - base, iszero(base) ? nothing : 100 * (value - base) / base,
                linearity, skin.sheath, skin.armour, skin.armour_length))
        end
    end
    resistance = observe(result.coefficients, R)[:, :, 1] .* 1000
    conductance = observe(result.coefficients, G)[:, :, 1] .* 1000
    for (i, row) in enumerate(terminal_order)
        modal_csv_row(conductors, (point, frequencies[1], row.cable, row.terminal,
            row.cable, row.terminal, "R_ohm_per_km", resistance[i, i]))
        for (j, column) in enumerate(terminal_order)
            modal_csv_row(conductors, (point, frequencies[1], row.cable, row.terminal,
                column.cable, column.terminal, "G_S_per_km", conductance[i, j]))
        end
    end
    observed = result.curves
    metadata = (; point, inputs = p, diagnostics, baseline_mode = matched.assignment,
        overlap = matched.assigned, best_overlap = matched.best, second_overlap = matched.second,
        ambiguous = matched.ambiguous, match_frequency_Hz = frequencies[match_index],
        mode_class_frequencies_Hz = frequencies[sampled],
        mode_classes = [modal_class(result.currents[:, mode, k], terminal_order)
            for mode in axes(result.currents, 2), k in axes(result.currents, 3)],
        sheath_transition_Hz = skin.sheath, armour_transition_Hz = skin.armour,
        armour_skin_length_m = skin.armour_length)
    return ObservedResult(merge(observed.gridpoint, (; study = metadata)),
        observed.quantities, observed.errors, observed.timings)
end

function modal_output_tables(folder, name)
    table = open(joinpath(folder, name * ".csv"), "w")
    conductors = open(joinpath(folder, name * "_conductors.csv"), "w")
    modal_csv_row(table, ("point", "mode", "baseline_mode", "overlap", "best_overlap",
        "second_overlap", "ambiguous", "mode_class", "frequency_Hz", "quantity", "value",
        "baseline", "change", "change_percent", "nonlinearity_fraction",
        "sheath_transition_Hz", "armour_transition_Hz", "armour_skin_length_m"))
    modal_csv_row(conductors, ("point", "frequency_Hz", "row_cable", "row_terminal",
        "column_cable", "column_terminal", "quantity", "value"))
    return table, conductors
end

function modal_nonlinearity(results, baseline, length, sampled, match_index)
    quantities = map(result -> modal_quantities(result, length), results)
    matches = map(result -> modal_match(result, baseline, match_index), results)
    output = map(values -> zeros(size(values, 1), Base.length(sampled)), quantities[1])
    for quantity in keys(output), mode in axes(output[quantity], 1), (j, k) in enumerate(sampled)
        curve = [quantities[i][quantity][findfirst(==(mode), matches[i].assignment), k]
            for i in 1:Base.length(results)]
        straight = range(first(curve), last(curve); length = Base.length(curve))
        deviation = maximum(abs.(curve .- straight))
        span = abs(last(curve) - first(curve))
        output[quantity][mode, j] = iszero(span) ? (iszero(deviation) ? 0.0 : NaN) : deviation / span
    end
    return output
end

function read_modal(path)
    curves, bases, coefficients = mktempdir() do temporary
        expanded = joinpath(temporary, "modal.json")
        run(pipeline(`gzip -dc -- $path`; stdout = expanded))
        import_data(:observed, expanded)
    end
    metadata = curves.gridpoint.modal
    currents = complex.(observe(bases, (Ti, real)), observe(bases, (Ti, imag)))
    result = (; curves, coefficients, currents, metadata.basis_indices, metadata.diagnostics)
    return (; result, metadata.point, metadata.inputs, metadata.study)
end

function modal_tables(run_folder, output_root)
    started = time_ns()
    metadata = import_data(:observed, joinpath(run_folder, "modal_run.json")).gridpoint
    folder = joinpath(output_root, basename(dirname(run_folder)), basename(run_folder))
    mkpath(folder)
    frequencies = metadata.nominal.frequencies
    sampled = [only(findall(f -> isapprox(f, requested; rtol = 1e-12), frequencies))
        for requested in (1.0, 50.0, 1e3, 1e4, 1e5, 1e6)]
    match_index = sampled[4]
    baseline = read_modal(joinpath(run_folder, "modal_baseline_0001.json.gz")).result
    summary = open(joinpath(folder, "modal_summary.txt"), "w")
    continuity_tables = (fallback = open(joinpath(folder, "modal_fallback_continuity.csv"), "w"),
        order = open(joinpath(folder, "modal_mode_order.csv"), "w"))
    modal_csv_row(continuity_tables.fallback, ("study", "point", "mode", "baseline_mode",
        "fallback_index", "fallback_frequency_Hz", "from_index", "to_index", "from_Hz", "to_Hz",
        "quantity", "relative_jump", "baseline_first_index", "baseline_last_index",
        "baseline_max_relative_jump", "baseline_ratio", "exceeds_threefold"))
    modal_csv_row(continuity_tables.order, ("study", "point", "frequency_Hz", "mode", "baseline_mode",
        "overlap", "best_overlap", "second_overlap", "differs_from_identity", "overlap_below_0_9", "ambiguous"))
    println(summary, "Source: ", run_folder)
    println(summary, "Modal cable study; ", length(frequencies), " frequencies; length_m=", metadata.nominal.line_length)
    println(summary, "Classification: root sum square terminal amplitudes, opposite phase >=135 degrees,",
        " local pair >=0.5 and remaining amplitude <=0.2; net fraction >=0.5 for dominant return.")
    println(summary, "Cable parity: sign of real current overlap; three cables use all active pairs.")
    println(summary, "Matches: normalized current overlap at 10 kHz, Hungarian assignment;",
        " ambiguous when best or assigned <0.9, or best minus second <=0.05.")
    println(summary, "Skin depth: rho(T)/(pi*mu*thickness^2); armour uses wire diameter or tube thickness.")
    println(summary, "R and G use phase coefficients at the lowest frequency; R is the DC estimate.")
    println(summary, "Travel time is phase travel time; zero baseline percentages are unavailable.")
    println(summary, "Nonlinearity is maximum deviation from the endpoint line divided by endpoint span;",
        " zero span with nonzero deviation is NaN.")
    println(summary, "Fallback jumps: abs(after-before)/abs(before), divided by the largest adjacent",
        " jump of the same numbered baseline mode within indices k-2:k+2, clipped to the scan.")
    println(summary, "The threshold is strictly greater than three times that baseline maximum;",
        " zero divided by zero is reported as zero, positive divided by zero as Inf.")
    println(summary, "Mode order: Hungarian current-basis assignment at 1 Hz, 50 Hz and 1 kHz;",
        " report nonidentity assignments and assigned overlaps below 0.9.")
    try
        for point in metadata.infeasible
            println(summary, "infeasible ", point, " cable exteriors overlap")
        end
        open(joinpath(folder, "modal_corrections_status.csv"), "w") do io
            modal_csv_row(io, ("correction", "status"))
            for (name, status) in metadata.correction_status
                modal_csv_row(io, (name, status))
                println(summary, "correction ", name, ": ", status)
            end
        end
        open(joinpath(folder, "modal_earth_formulas_registry.csv"), "w") do io
            modal_csv_row(io, ("family", "registered", "resolved", "description", "availability"))
            for row in metadata.registry
                modal_csv_row(io, row)
            end
        end
        for timing in metadata.timings
            study = timing.study
            section_start = time_ns()
            paths = sort(filter(path -> startswith(basename(path), "modal_" * study * "_") &&
                endswith(path, ".json.gz"), readdir(run_folder; join = true)))
            length(paths) == timing.points || error("Archive count differs from the declared grid: $study")
            table, conductors = modal_output_tables(folder, "modal_" * study)
            observations = ObservedResult[]
            try
                batches = study == "one_at_a_time" ? Iterators.partition(paths, 9) : (paths,)
                for batch in batches
                    records = read_modal.(batch)
                    nonlinearities = study == "one_at_a_time" ? modal_nonlinearity(
                        getproperty.(records, :result), baseline, metadata.nominal.line_length, sampled, match_index) : nothing
                    for record in records
                        push!(observations, modal_write_point(table, conductors, summary,
                            record.result, baseline, record.inputs, record.point,
                            metadata.terminal_order, frequencies, sampled, match_index;
                            continuity_tables, study, nonlinearities))
                    end
                end
                archive = joinpath(folder, "modal_" * study * ".json")
                save(study == "baseline" ? only(observations) : observations, archive)
                run(`gzip -n -9 -- $archive`)
            finally
                close(table)
                close(conductors)
            end
            elapsed = (time_ns() - section_start) / 1e9
            println(summary, study, " wall_seconds=", timing.wall_seconds, " points=", timing.points,
                " generation_seconds=", elapsed)
            println(basename(dirname(run_folder)), " ", study, ": generation ", elapsed, " s")
            flush(summary)
            flush(stdout)
        end
        println(summary, "total wall_seconds=", metadata.wall_seconds)
        println(summary, "generation wall_seconds=", (time_ns() - started) / 1e9)
    finally
        close(continuity_tables.fallback)
        close(continuity_tables.order)
        close(summary)
    end
    println("Generated ", folder, " in ", (time_ns() - started) / 1e9, " s")
    flush(stdout)
    return folder
end

length(ARGS) >= 2 || error("Usage: julia --project=test benchmarks/modal_tables.jl <run folder>... <output folder>")
for run_folder in ARGS[1:end-1]
    modal_tables(abspath(run_folder), abspath(ARGS[end]))
end
