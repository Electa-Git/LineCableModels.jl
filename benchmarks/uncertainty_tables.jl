using LineCableModels
using Measurements
using Statistics: mean, std, quantile
using Random: Xoshiro, rand
using LinearAlgebra: BLAS, norm
using Printf: @sprintf

function csv_row(io, values)
    println(io, join((ismissing(value) || value === nothing ? "" :
        "\"" * replace(string(value), "\"" => "\"\"") * "\"" for value in values), ','))
end

function physical_products(products, frequencies)
    haskey(products, :Z) && return products
    n, _, nf, trials = size(products.R)
    impedance = Array{ComplexF64}(undef, n, n, nf, trials)
    admittance = similar(impedance)
    for trial in 1:trials, k in 1:nf
        voltage_basis = products.Tv[:, :, k, trial]
        current_basis = products.Ti[:, :, k, trial]
        angular = 2pi * frequencies[k]
        impedance[:, :, k, trial] = voltage_basis *
            complex.(products.R[:, :, k, trial], angular .* products.L[:, :, k, trial]) / current_basis
        admittance[:, :, k, trial] = current_basis *
            complex.(products.G[:, :, k, trial], angular .* products.C[:, :, k, trial]) / voltage_basis
    end
    return merge(products, (; Z = impedance, Y = admittance))
end

function engineering_products(products, line_length)
    result = NamedTuple[]
    if haskey(products, :Z)
        for (name, quantity, unit) in ((:Z, products.Z, "ohm/m"), (:Y, products.Y, "S/m")),
                (component, operation) in (("real", real), ("imag", imag), ("magnitude", abs))
            push!(result, (; name = "$(name)_$(component)", values = operation.(quantity), unit, matrix = true))
        end
    end
    if haskey(products, :velocity)
        for (name, values, unit) in (
                ("attenuation", products.alpha .* (1000 * 20 / log(10)), "dB/km"),
                ("velocity", products.velocity ./ 1e6, "m/us"),
                ("surge_impedance", abs.(products.Zc), "ohm"),
                ("travel_time", 1e6line_length ./ products.velocity, "us"))
            push!(result, (; name, values, unit, matrix = false))
        end
    end
    return result
end

function mode_description(current_basis, terminal_order, mode, k)
    column = current_basis[:, mode, k]
    net = abs(sum(column)) / sum(abs, column)
    amplitudes = Dict(terminal => norm(column[findall(entry -> entry.terminal == terminal, terminal_order)])
        for terminal in unique(entry.terminal for entry in terminal_order))
    terminal = first(sort!(collect(keys(amplitudes)); by = t -> amplitudes[t], rev = true))
    return net >= 0.5 ? "$(terminal) return" : "balanced $(terminal)"
end

function comparison_classes(sigma, spread, nominal, ratio_low, ratio_high, shift, interval)
    scale = abs(nominal)
    sigma_class = if ismissing(sigma)
        "unsupported"
    elseif !isfinite(sigma) || !isfinite(spread)
        "unavailable"
    elseif iszero(spread) || spread < 0.001scale
        "insensitive"
    elseif ismissing(ratio_low) || ismissing(ratio_high)
        "inconclusive"
    elseif 0.9 <= ratio_low <= ratio_high <= 1.1
        "equivalent"
    elseif ratio_high < 0.9 || ratio_low > 1.1
        "different"
    else
        "inconclusive"
    end
    mean_shift_flag = !(interval.shift[1] <= 0 <= interval.shift[2]) && abs(shift) > 0.01scale
    mean_class = mean_shift_flag ? "mean shift" :
        maximum(abs, interval.shift) <= 0.01scale ? "within 1%" : "inconclusive"
    return (; sigma_class, mean_class, mean_shift_flag)
end

function sample_rows(baseline, sampled, methods, declaration)
    rows = NamedTuple[]
    trials = size(sampled.Z, 4)
    data = Vector{Vector{Float64}}()
    base_products = engineering_products(baseline, declaration.nominal.line_length)
    sampled_products = engineering_products(sampled, declaration.nominal.line_length)
    method_products = map(methods) do (method, point)
        method => engineering_products(point.gridpoint.products, declaration.nominal.line_length)
    end
    for (base, product) in zip(base_products, sampled_products)
        for coordinate in CartesianIndices(base.values)
            i = coordinate.I
            k = last(i)
            # Modal engineering tables use the six report frequencies.
            !base.matrix && k == 1 && continue
            selected = vec(product.values[i..., :])
            sigmas = Dict(method => begin
                matching = findfirst(item -> item.name == base.name, quantities)
                matching === nothing ? missing : Measurements.uncertainty(quantities[matching].values[i...])
            end for (method, quantities) in method_products)
            mode = base.matrix ? nothing : first(i)
            family = base.matrix ? "" : mode_description(baseline.Ti, declaration.terminal_order, mode, k)
            entry_class, conductor_pair, physical_class, entry_magnitude, matrix_scale = if base.matrix
                left, right = declaration.terminal_order[i[1]], declaration.terminal_order[i[2]]
                entry_class = left.cable != right.cable ? "between-cable" :
                    i[1] == i[2] ? "within-cable self" : "within-cable mutual"
                terminals = unique(entry.terminal for entry in declaration.terminal_order)
                pair = sort!([left.terminal, right.terminal]; by = t -> findfirst(==(t), terminals))
                conductor_pair = replace(join(string.(pair), '/'), "armor" => "armour")
                matrix = startswith(base.name, "Z_") ? baseline.Z : baseline.Y
                (entry_class, conductor_pair, "$entry_class $conductor_pair",
                    abs(matrix[i...]), maximum(abs, view(matrix, :, :, k)))
            else
                ("modal", "", "modal mode $mode", nothing, nothing)
            end
            structural_zero = base.matrix && (iszero(matrix_scale) || entry_magnitude < 1e-9matrix_scale)
            push!(rows, (; observable = base.name, frequency_hz = declaration.frequencies[k],
                row = first(i), column = base.matrix ? i[2] : nothing, mode, family,
                entry_class, conductor_pair, physical_class, entry_magnitude, matrix_scale,
                structural_zero, unit = base.unit, nominal = Float64(base.values[i...]), sigmas))
            push!(data, Float64.(selected))
        end
    end
    return rows, permutedims(reduce(hcat, data))
end

function transient_rows(reference, sampled, central)
    baseline = reference.baseline
    core_index = findfirst(==(reference.core1), reference.coordinates)
    rows = NamedTuple[]
    data = Vector{Vector{Float64}}()
    for (name, selector, unit) in ((:peak, :peak, "V"), (:peak_time, :peak_time, "s"),
            (:arrival_time, :arrival_time, "s"))
        for side in 1:2, (terminal, label) in enumerate(reference.coordinates)
            selected = name == :peak ? (endswith(label, ":sheath") || endswith(label, ":armor") ||
                (side == 2 && terminal == core_index)) : side == 2 && terminal == core_index
            selected || continue
            push!(rows, (; observable = string(name), frequency_hz = nothing, row = side,
                column = terminal, mode = nothing,
                family = (side == 1 ? "sending " : "receiving ") * label, unit,
                entry_class = "transient", conductor_pair = "",
                physical_class = "transient " * (side == 1 ? "sending " : "receiving ") * label,
                entry_magnitude = nothing, matrix_scale = nothing, structural_zero = false,
                nominal = Float64(getproperty(baseline, selector)[side, terminal]),
                sigmas = Dict(:central => begin
                    value = getproperty(central, selector)[side, terminal]
                    ismissing(value) ? missing : Measurements.uncertainty(value)
                end)))
            push!(data, Float64.(vec(getproperty(sampled, selector)[side, terminal, :])))
        end
    end
    return rows, permutedims(reduce(hcat, data))
end

function bootstrap_intervals(data, nominal, n; repetitions = 2000)
    # One resampled trial index serves each observable and preserves output dependence.
    rng = Xoshiro(20261018 + n)
    weights = zeros(n, repetitions)
    for repetition in 1:repetitions, _ in 1:n
        weights[rand(rng, 1:n), repetition] += 1
    end
    centered = data[:, 1:n] .- nominal
    means = centered * weights ./ n
    second_moments = abs2.(centered) * weights
    variances = max.(0.0, (second_moments .- n .* abs2.(means)) ./ (n - 1))
    deviations = sqrt.(variances)
    return [(shift = quantile(vec(means[j, :]), [0.025, 0.975]),
        sigma = quantile(vec(deviations[j, :]), [0.025, 0.975])) for j in axes(data, 1)]
end

function latex_escape(value)
    replace(string(value), '\\' => "\\textbackslash{}", '_' => "\\_", '%' => "\\%", '&' => "\\&", '#' => "\\#")
end
format_number(value) = ismissing(value) || value === nothing || !isfinite(value) ? "--" : @sprintf("%.3g", value)

function write_table(path, rows, caption)
    isempty(rows) && return
    headers = keys(first(rows))
    open(path * ".csv", "w") do io
        csv_row(io, headers)
        foreach(row -> csv_row(io, values(row)), rows)
    end
    open(path * ".tex", "w") do io
        println(io, "\\begin{longtable}{lllllrrrrl}")
        println(io, "\\caption{", caption, "}\\\\")
        println(io, "\\toprule")
        println(io, "Case & Physical class & Observable & Coordinate & Method & N & Shift (\\%) & Spread & LEP/MC [95\\% CI] & Class / mean \\\\")
        println(io, "\\midrule\\endhead")
        for row in rows
            coordinate = row.mode === nothing ? "$(row.row),$(row.column)" : "mode $(row.mode)"
            ratio = "$(format_number(row.ratio)) [$(format_number(row.ratio_low)), $(format_number(row.ratio_high))]"
            shift = "$(format_number(row.shift_percent)) [$(format_number(row.shift_low_percent)), $(format_number(row.shift_high_percent))]"
            spread = "$(format_number(row.mc_sigma)) [$(format_number(row.sigma_low)), $(format_number(row.sigma_high))]"
            frequency = row.frequency_hz === nothing ? "" : " @ $(format_number(row.frequency_hz)) Hz"
            observable = "$(row.observable) ($(row.unit))$(frequency)"
            println(io, join(latex_escape.(("$(row.nominal_id): $(row.study) $(row.percent)%", row.physical_class, observable,
                coordinate, row.method, row.n, shift, spread, ratio,
                row.structural_zero ? "structural zero" : "$(row.sigma_class) / $(row.mean_class)")), " & "), " \\\\")
        end
        println(io, "\\bottomrule\n\\end{longtable}")
    end
end

function class_summaries(rows; by_observable = false)
    groups = Dict{Tuple,Vector{NamedTuple}}()
    for row in rows
        row.n == row.accepted || continue
        key = (row.study, row.percent, row.nominal_id, string(row.method), row.physical_class,
            by_observable ? row.observable : "all")
        push!(get!(groups, key, NamedTuple[]), row)
    end
    map(sort!(collect(keys(groups)); by = string)) do key
        group = groups[key]
        active = filter(row -> !row.structural_zero && !(row.sigma_class in ("unsupported", "unavailable")), group)
        first_row = first(group)
        (; study = key[1], percent = key[2], nominal_id = key[3], method = key[4],
            first_row.mc_source_folder, first_row.lep_source_folder, first_row.baseline_source_folder,
            physical_class = key[5], observable = key[6], first_row.entry_class,
            first_row.conductor_pair, first_row.mode, family = join(unique(row.family for row in group), "; "),
            n = first_row.accepted, comparisons = length(active),
            equivalent = count(row -> row.sigma_class == "equivalent", active),
            different = count(row -> row.sigma_class == "different", active),
            inconclusive = count(row -> row.sigma_class == "inconclusive", active),
            insensitive = count(row -> row.sigma_class == "insensitive", active),
            structural_zeros = count(row -> row.structural_zero, group),
            unsupported = count(row -> !row.structural_zero && row.sigma_class == "unsupported", group),
            unavailable = count(row -> !row.structural_zero && row.sigma_class == "unavailable", group),
            mean_shift = count(row -> row.mean_shift_flag, active),
            mean_within_one_percent = count(row -> row.mean_class == "within 1%", active),
            mean_inconclusive = count(row -> row.mean_class == "inconclusive", active),
            warm_speedup = first_row.speedup,
            sigma_verdict = isempty(active) ? "excluded" :
                any(row -> row.sigma_class == "different", active) ? "differences present" :
                any(row -> row.sigma_class == "inconclusive", active) ? "inconclusive" :
                all(row -> row.sigma_class == "insensitive", active) ? "insensitive" : "equivalent or insensitive",
            mean_verdict = isempty(active) ? "excluded" : any(row -> row.mean_shift_flag, active) ? "mean shifts present" :
                all(row -> row.mean_class == "within 1%", active) ? "within 1%" : "inconclusive")
    end
end

function lep_comparisons(rows)
    groups = Dict{Tuple,Dict{Symbol,NamedTuple}}()
    for row in rows
        key = (row.study, row.percent, row.nominal_id, row.n, row.observable, row.frequency_hz,
            row.row, row.column, row.mode)
        get!(groups, key, Dict{Symbol,NamedTuple}())[row.method] = row
    end
    comparisons = NamedTuple[]
    for key in sort!(collect(keys(groups)); by = string)
        group = groups[key]
        haskey(group, :measurements) && haskey(group, :central) || continue
        measured, central = group[:measurements], group[:central]
        ratio = if ismissing(measured.lep_sigma) || ismissing(central.lep_sigma) ||
                !isfinite(measured.lep_sigma) || !isfinite(central.lep_sigma)
            missing
        elseif iszero(measured.lep_sigma)
            iszero(central.lep_sigma) ? 1.0 : Inf
        else
            central.lep_sigma / measured.lep_sigma
        end
        disagree = !ismissing(ratio) && !(0.9 <= ratio <= 1.1)
        m, c = measured.sigma_class, central.sigma_class
        support = if measured.structural_zero
            "structural zero"
        elseif m == "insensitive" || c == "insensitive"
            "insensitive"
        elseif m in ("unsupported", "unavailable") || c in ("unsupported", "unavailable")
            "unavailable"
        elseif m == "equivalent" && c == "different"
            "measurements"
        elseif c == "equivalent" && m == "different"
            "central"
        elseif m == c == "equivalent"
            "both"
        elseif m == c == "different"
            "neither"
        else
            "inconclusive"
        end
        push!(comparisons, (; measured.study, measured.percent, measured.nominal_id, measured.observable,
            measured.mc_source_folder, measured.baseline_source_folder,
            measurements_source_folder = measured.lep_source_folder, central_source_folder = central.lep_source_folder,
            measured.frequency_hz, measured.row, measured.column, measured.mode, measured.family,
            measured.entry_class, measured.conductor_pair, measured.physical_class, measured.structural_zero,
            measured.n, measured.accepted, measurements_sigma = measured.lep_sigma, central_sigma = central.lep_sigma,
            central_over_measurements = ratio, methods_disagree = disagree,
            measurements_class = m, central_class = c, mc_support = support,
            measurements_speedup = measured.speedup, central_speedup = central.speedup))
    end
    return comparisons
end

function lep_method_summaries(rows)
    groups = Dict{Tuple,Vector{NamedTuple}}()
    for row in rows
        row.n == row.accepted || continue
        key = (row.study, row.percent, row.nominal_id, row.physical_class)
        push!(get!(groups, key, NamedTuple[]), row)
    end
    map(sort!(collect(keys(groups)); by = string)) do key
        group = groups[key]
        eligible = filter(row -> !row.structural_zero && !(row.mc_support in ("insensitive", "unavailable")), group)
        differences = filter(row -> row.methods_disagree, eligible)
        first_row = first(group)
        (; study = key[1], percent = key[2], nominal_id = key[3], physical_class = key[4],
            first_row.mc_source_folder, first_row.baseline_source_folder,
            first_row.measurements_source_folder, first_row.central_source_folder,
            n = first_row.accepted, comparisons = length(eligible), methods_disagree = length(differences),
            mc_support_measurements = count(row -> row.mc_support == "measurements", differences),
            mc_support_central = count(row -> row.mc_support == "central", differences),
            mc_support_both = count(row -> row.mc_support == "both", differences),
            mc_support_neither = count(row -> row.mc_support == "neither", differences),
            mc_support_inconclusive = count(row -> row.mc_support == "inconclusive", differences),
            first_row.measurements_speedup, first_row.central_speedup)
    end
end

function write_class_summary(path, rows)
    isempty(rows) && return
    open(path * ".csv", "w") do io
        csv_row(io, keys(first(rows)))
        foreach(row -> csv_row(io, values(row)), rows)
    end
    open(path * ".tex", "w") do io
        println(io, "\\begin{longtable}{llllrrrrrlr}")
        println(io, "\\caption{Classes at the complete accepted sample. E/D/U/I denote equivalent, " *
            "different, inconclusive and insensitive spreads. The counts exclude structural zeros and unavailable methods.}\\\\")
        println(io, "\\toprule Study & Class & Observable & Method & N & E/D/U/I & Zeros & Mean shift & Mean within 1\\% & Mean verdict & Speed-up \\\\")
        println(io, "\\midrule\\endhead")
        for row in rows
            counts = "$(row.equivalent)/$(row.different)/$(row.inconclusive)/$(row.insensitive)"
            println(io, join(latex_escape.(("$(row.study) $(row.percent)%", row.physical_class, row.observable,
                row.method, row.n, counts, row.structural_zeros, row.mean_shift, row.mean_within_one_percent,
                row.mean_verdict, format_number(row.warm_speedup))), " & "), " \\\\")
        end
        println(io, "\\bottomrule\n\\end{longtable}")
    end
end

function generate_uncertainty_tables(source, destination; transient_source = source)
    declaration = import_data(:observed, joinpath(source, "uncertainty_run.json")).gridpoint
    if transient_source != source
        transient_declaration = import_data(:observed, joinpath(transient_source, "uncertainty_run.json")).gridpoint
        transient_declaration.script == declaration.script && transient_declaration.nominal == declaration.nominal ||
            error("Transient and frequency-domain archives must describe the same cable and nominal inputs")
    end
    default_baseline = import_data(:observed, joinpath(source, "uncertainty_baseline.json")).gridpoint.baseline
    output = mkpath(joinpath(destination, declaration.script, basename(source)))
    tables = Dict("uncertainty" => NamedTuple[], "uncertainty_kink" => NamedTuple[],
        "uncertainty_transient" => NamedTuple[])
    timing_rows = NamedTuple[]
    diagnostic_rows = NamedTuple[]
    status_rows = NamedTuple[]
    kink_samples = Dict{String,Any}()
    archives = [(source, file) for file in sort(readdir(source))]
    if transient_source != source
        filter!(entry -> last(entry) != "uncertainty_transient_10_monte_carlo.json", archives)
        push!(archives, (transient_source, "uncertainty_transient_10_monte_carlo.json"))
    end
    for (origin, file) in archives
        if endswith(file, "_status.json")
            point = import_data(:observed, joinpath(origin, file)).gridpoint
            push!(status_rows, (; source_folder = origin, point.study, point.percent, point.method, point.status, point.message))
            continue
        end
        endswith(file, "_monte_carlo.json") || continue
        mc = import_data(:observed, joinpath(origin, file))
        meta = mc.gridpoint
        methods = Pair{Symbol,Any}[]
        method_sources = Dict{Symbol,String}()
        for method in (:measurements, :central)
            method_file = replace(file, "monte_carlo" => string(method))
            method_source = isfile(joinpath(source, method_file)) ? source : origin
            path = joinpath(method_source, method_file)
            if isfile(path)
                push!(methods, method => import_data(:observed, path))
                method_sources[method] = method_source
            end
        end
        baseline_file = meta.study == "transient" ? "uncertainty_transient_baseline.json" :
            get(meta, :baseline_file, "uncertainty_baseline.json")
        baseline_source = isfile(joinpath(source, baseline_file)) ? source : origin
        rows, data = if meta.study == "transient"
            reference = import_data(:observed, joinpath(baseline_source, baseline_file)).gridpoint
            transient_rows(reference, meta.products, only(last(point).gridpoint.products
                for point in methods if first(point) == :central))
        else
            baseline = baseline_file == "uncertainty_baseline.json" ? default_baseline :
                import_data(:observed, joinpath(baseline_source, baseline_file)).gridpoint.baseline
            products = physical_products(meta.products, declaration.frequencies)
            sample_rows(baseline, products, methods, declaration)
        end
        accepted = size(data, 2)
        if startswith(meta.study, "kink")
            kink_samples[meta.study] = (; rows, data, inputs = meta.inputs,
                complete = accepted == meta.trials && meta.diagnostics.failure_summary.failed == 0)
        end
        counts = sort!(unique([filter(n -> n <= accepted, [100, 200, 500, 1000, 2000]); accepted]))
        target = meta.study == "transient" ? "uncertainty_transient" :
            startswith(meta.study, "kink") ? "uncertainty_kink" : "uncertainty"
        for n in counts
            intervals = bootstrap_intervals(data, [row.nominal for row in rows], n)
            for (j, row) in enumerate(rows)
                spread = std(view(data, j, 1:n))
                shift = mean(view(data, j, 1:n)) - row.nominal
                interval = intervals[j]
                scale = abs(row.nominal)
                shift_ok = interval.shift[1] <= 0 <= interval.shift[2] ||
                    maximum(abs, interval.shift) < 0.01scale
                for (method, lep) in methods
                    sigma = get(row.sigmas, method, missing)
                    ratio = ismissing(sigma) || !isfinite(sigma) || !isfinite(spread) || spread <= 0 ? missing : sigma / spread
                    ratio_low = ismissing(ratio) || interval.sigma[2] <= 0 ? missing : sigma / interval.sigma[2]
                    ratio_high = ismissing(ratio) || interval.sigma[1] <= 0 ? missing : sigma / interval.sigma[1]
                    agrees = !ismissing(ratio_low) && !ismissing(ratio_high) &&
                        0.9 <= ratio_low <= ratio_high <= 1.1 && shift_ok
                    classes = comparison_classes(sigma, spread, row.nominal, ratio_low, ratio_high, shift, interval)
                    percent(value) = iszero(scale) ? missing : 100value / scale
                    first_pass = lep.timings.compilation_seconds
                    warm = lep.timings.warm_seconds
                    mc_warm = mc.timings.warm_seconds
                    push!(tables[target], (; study = meta.study, percent = meta.percent,
                        mc_source_folder = origin, lep_source_folder = method_sources[method],
                        baseline_source_folder = baseline_source,
                        nominal_id = get(meta, :nominal_id, "original"),
                        bedding_thickness_m = get(meta, :bedding_thickness,
                            get(declaration.nominal, :bedding_thickness, nothing)),
                        observable = row.observable, row.frequency_hz, row.row, row.column,
                        row.mode, row.family, row.unit, row.nominal, n, accepted, requested = meta.trials,
                        row.entry_class, row.conductor_pair, row.physical_class, row.entry_magnitude,
                        row.matrix_scale, row.structural_zero, classes...,
                        mean_shift = shift, shift_low = interval.shift[1], shift_high = interval.shift[2],
                        shift_percent = percent(shift), shift_low_percent = percent(interval.shift[1]),
                        shift_high_percent = percent(interval.shift[2]), mc_sigma = spread,
                        sigma_low = interval.sigma[1], sigma_high = interval.sigma[2],
                        method, lep_sigma = sigma, ratio, ratio_low, ratio_high, agrees,
                        status = ismissing(sigma) ? "unsupported" : !isfinite(sigma) ? "nonfinite LEP" :
                            ismissing(ratio) ? "undefined spread" : "available",
                        mc_first_seconds = mc.timings.compilation_seconds, mc_warm_seconds = mc_warm,
                        lep_first_seconds = first_pass, lep_warm_seconds = warm,
                        speedup = mc_warm / warm,
                        subset_mc_seconds_estimate = mc_warm * n / meta.trials,
                        subset_speedup_estimate = mc_warm * n / meta.trials / warm))
                end
            end
        end
        for (method, point) in [:monte_carlo => mc; methods]
            timing = point.timings
            target_n = meta.study == "joint" && meta.percent == 10 ? 2000 :
                startswith(meta.study, "kink") ? 1000 : 500
            multiplier = method == :monte_carlo ? target_n / meta.trials : 1.0
            overhead = max(0, timing.compilation_seconds - timing.warm_seconds)
            push!(timing_rows, (; study = meta.study, percent = meta.percent, method,
                source_folder = method == :monte_carlo ? origin : method_sources[method],
                requested = meta.trials, accepted, first_seconds = timing.compilation_seconds,
                warm_seconds = timing.warm_seconds, compile_overhead_estimate = overhead,
                warm_speedup = mc.timings.warm_seconds / timing.warm_seconds,
                full_trials = target_n, full_first_plus_warm_estimate = overhead + 2multiplier * timing.warm_seconds))
        end
        diagnostic = meta.diagnostics
        assignment = get(diagnostic, :mode_assignment, nothing)
        clustered = assignment !== nothing && haskey(assignment, :cluster_sizes) ?
            vec(any(>(1), assignment.cluster_sizes; dims = 1)) : nothing
        push!(diagnostic_rows, (; source_folder = origin,
            study = meta.study, percent = meta.percent, requested = meta.trials, accepted,
            failures = haskey(diagnostic, :failure_summary) ? diagnostic.failure_summary.failed : 0,
            minimum_overlap = haskey(diagnostic, :mode_assignment) ? minimum(diagnostic.mode_assignment.overlap) : missing,
            assignment_frequencies = clustered === nothing ? missing : length(assignment.frequencies),
            clustered_frequency_count = clustered === nothing ? missing : count(clustered),
            clustered_frequencies_hz = clustered === nothing ? missing : join(assignment.frequencies[clustered], ';'),
            modal_fallback_trials = haskey(diagnostic, :modal) ? count(point -> !isempty(point.fallback_frequencies), diagnostic.modal) : missing,
            modal_missed_trials = haskey(diagnostic, :modal) ? count(point -> !isempty(point.missed_frequencies), diagnostic.modal) : missing,
            fixed_radius = get(meta, :fixed_radius, nothing), inputs = join(get(meta, :inputs, ()), ';')))
        println("tables: ", meta.study, " ±", meta.percent, "%; ", accepted, "/", meta.trials, " trials")
        flush(stdout)
    end
    effects = NamedTuple[]
    for subset in ("all", "under", "remaining")
        rule_key, fixed_key = "kink_$(subset)_rule", "kink_$(subset)_fixed"
        haskey(kink_samples, rule_key) && haskey(kink_samples, fixed_key) || continue
        rule, fixed = kink_samples[rule_key], kink_samples[fixed_key]
        # Pair identical seeded input draws only when neither run excludes trials.
        rule.complete && fixed.complete && rule.inputs == fixed.inputs &&
            size(rule.data) == size(fixed.data) || continue
        trials = size(rule.data, 2)
        counts = sort!(unique([filter(n -> n <= trials, [100, 200, 500, 1000, 2000]); trials]))
        difference = fixed.data - rule.data
        for n in counts
            intervals = bootstrap_intervals(difference, zeros(size(difference, 1)), n)
            for (j, row) in enumerate(rule.rows)
                scale = abs(row.nominal)
                percent(value) = iszero(scale) ? missing : 100value / scale
                change = mean(view(difference, j, 1:n))
                push!(effects, (; source_folder = source, subset, row.observable, row.frequency_hz, row.row, row.column,
                    row.mode, row.family, row.unit, row.nominal, n,
                    rule_shift = mean(view(rule.data, j, 1:n)) - row.nominal,
                    fixed_shift = mean(view(fixed.data, j, 1:n)) - row.nominal,
                    fixed_minus_rule = change, low = intervals[j].shift[1], high = intervals[j].shift[2],
                    change_percent = percent(change), low_percent = percent(intervals[j].shift[1]),
                    high_percent = percent(intervals[j].shift[2])))
            end
        end
    end
    if !isempty(effects)
        open(joinpath(output, "uncertainty_kink_effects.csv"), "w") do io
            csv_row(io, keys(first(effects)))
            foreach(row -> csv_row(io, values(row)), effects)
        end
        open(joinpath(output, "uncertainty_kink_effects.tex"), "w") do io
            println(io, "\\begin{longtable}{llllrr}")
            println(io, "\\caption{Paired effect of a fixed armour radius on the Monte Carlo mean. " *
                "Change is fixed minus ordinary rule, as a percentage of the absolute nominal value.}\\\\")
            println(io, "\\toprule Inputs & Observable & Coordinate & N & Change (\\%) & 95\\% interval \\\\")
            println(io, "\\midrule\\endhead")
            for row in effects
                coordinate = row.mode === nothing ? "$(row.row),$(row.column)" : "mode $(row.mode)"
                frequency = format_number(row.frequency_hz)
                interval = "[$(format_number(row.low_percent)), $(format_number(row.high_percent))]"
                println(io, join(latex_escape.((row.subset, "$(row.observable) @ $(frequency) Hz", coordinate,
                    row.n, format_number(row.change_percent), interval)), " & "), " \\\\")
            end
            println(io, "\\bottomrule\n\\end{longtable}")
        end
    end
    for (name, rows) in tables
        write_table(joinpath(output, name), rows,
            "Monte Carlo and first-order uncertainty for $(latex_escape(declaration.script)). " *
            "Intervals use 2000 percentile bootstrap replicates. Spread uses the stated native unit; " *
            "the CSV includes both uncertainty methods, timing, mode descriptions and status.")
    end
    all_rows = reduce(vcat, values(tables))
    summaries = class_summaries(all_rows)
    observable_summaries = class_summaries(all_rows; by_observable = true)
    write_class_summary(joinpath(output, "uncertainty_classes"), summaries)
    write_class_summary(joinpath(output, "uncertainty_observable_classes"), observable_summaries)
    method_comparisons = lep_comparisons(all_rows)
    method_summaries = lep_method_summaries(method_comparisons)
    structural_zeros = filter(row -> row.structural_zero && row.n == row.accepted, all_rows)
    for (name, rows) in (("uncertainty_lep_methods", method_comparisons),
            ("uncertainty_lep_method_classes", method_summaries),
            ("uncertainty_structural_zeros", structural_zeros))
        isempty(rows) && continue
        open(joinpath(output, name * ".csv"), "w") do io
            csv_row(io, keys(first(rows)))
            foreach(row -> csv_row(io, values(row)), rows)
        end
    end
    open(joinpath(output, "uncertainty_lep_method_classes.tex"), "w") do io
        println(io, "\\begin{longtable}{llrrrrrr}")
        println(io, "\\caption{LEP method differences above 10\\%, with Monte Carlo support. " *
            "The counts exclude structural zeros and insensitive or unavailable comparisons.}\\\\")
        println(io, "\\toprule Study & Class & Differences & Measurements & Central & Both & Neither & Inconclusive \\\\")
        println(io, "\\midrule\\endhead")
        for row in method_summaries
            println(io, join(latex_escape.(("$(row.study) $(row.percent)%", row.physical_class,
                row.methods_disagree, row.mc_support_measurements, row.mc_support_central,
                row.mc_support_both, row.mc_support_neither, row.mc_support_inconclusive)), " & "), " \\\\")
        end
        println(io, "\\bottomrule\n\\end{longtable}")
    end
    for (name, rows) in (("uncertainty_timings", timing_rows), ("uncertainty_diagnostics", diagnostic_rows),
            ("uncertainty_status", status_rows))
        isempty(rows) && continue
        open(joinpath(output, name * ".csv"), "w") do io
            csv_row(io, keys(first(rows)))
            foreach(row -> csv_row(io, values(row)), rows)
        end
    end
    open(joinpath(output, "uncertainty_timings.tex"), "w") do io
        println(io, "\\begin{longtable}{lllrrrrr}")
        println(io, "\\caption{Uncertainty computation cost in seconds. The first pass includes compilation. " *
            "Full estimates include both passes and scale the sample cost with trial count.}\\\\")
        println(io, "\\toprule Study & Level & Method & First & Warm & Compile estimate & Speed-up & Full estimate \\\\")
        println(io, "\\midrule\\endhead")
        for row in timing_rows
            println(io, join(latex_escape.((row.study, row.percent, row.method,
                format_number(row.first_seconds), format_number(row.warm_seconds),
                format_number(row.compile_overhead_estimate), format_number(row.warm_speedup),
                format_number(row.full_first_plus_warm_estimate))), " & "), " \\\\")
        end
        println(io, "\\bottomrule\n\\end{longtable}")
    end
    geometry_file = joinpath(source, "uncertainty_kink_geometry.json")
    geometry = isfile(geometry_file) ? import_data(:observed, geometry_file).gridpoint.geometry : ()
    geometry_rows = NamedTuple[]
    if !isempty(geometry)
        geometry_rows = [(; source_folder = source, entry.nominal_id, bedding_mm = 1000entry.bedding_thickness,
            nominal_margin_mm = 1000entry.nominal_margin, inner_margin_mm = 1000entry.inner_margin,
            all_margin_mm = 1000entry.all_margin) for entry in geometry]
        open(joinpath(output, "uncertainty_kink_geometry.csv"), "w") do io
            csv_row(io, keys(first(geometry_rows)))
            foreach(row -> csv_row(io, values(row)), geometry_rows)
        end
        open(joinpath(output, "uncertainty_kink_geometry.tex"), "w") do io
            println(io, "\\begin{longtable}{lrrrr}")
            println(io, "\\caption{Bedding and clearance above the wire-fit radius in millimetres. " *
                "Support margins include the bedding tolerance.}\\\\")
            println(io, "\\toprule Nominal & Bedding & Nominal margin & Inner-input margin & All-input margin \\\\")
            println(io, "\\midrule\\endhead")
            for row in geometry_rows
                println(io, join(latex_escape.((row.nominal_id,
                    (format_number(getproperty(row, key)) for key in
                        (:bedding_mm, :nominal_margin_mm, :inner_margin_mm, :all_margin_mm))...)), " & "), " \\\\")
            end
            println(io, "\\bottomrule\n\\end{longtable}")
        end
    end
    source_rows = NamedTuple[]
    for (name, rows) in [collect(pairs(tables));
            "uncertainty_classes" => summaries;
            "uncertainty_observable_classes" => observable_summaries;
            "uncertainty_lep_methods" => method_comparisons;
            "uncertainty_lep_method_classes" => method_summaries;
            "uncertainty_structural_zeros" => structural_zeros;
            "uncertainty_timings" => timing_rows; "uncertainty_diagnostics" => diagnostic_rows;
            "uncertainty_status" => status_rows; "uncertainty_kink_effects" => effects;
            "uncertainty_kink_geometry" => geometry_rows]
        origins = sort!(unique([getproperty(row, key) for row in rows for key in keys(row)
            if endswith(string(key), "source_folder")]))
        for extension in (".csv", ".tex"), origin in origins
            isfile(joinpath(output, name * extension)) || continue
            push!(source_rows, (; table = name * extension, source_folder = origin))
        end
    end
    open(joinpath(output, "uncertainty_sources.csv"), "w") do io
        csv_row(io, (:table, :source_folder))
        foreach(row -> csv_row(io, values(row)), source_rows)
    end
    open(joinpath(output, "uncertainty_summary.txt"), "w") do io
        println(io, "Frequency-domain source: ", abspath(source))
        for row in source_rows
            println(io, "Table ", row.table, ": ", row.source_folder)
        end
        println(io, "Bootstrap: 2000 paired trial resamples; seed 20261018 + subset count.")
        println(io, "Full compute estimate, first and warm passes: ",
            sum(row.full_first_plus_warm_estimate for row in timing_rows), " s")
        println(io, "Complete samples only. Counts include directed entries, cable copies and observable components. " *
            "They are not independent experiments.")
        println(io, "Spread counts use E for equivalent and D for different. U denotes inconclusive and I denotes insensitive. " *
            "Insensitivity means MC sigma below 0.1% of the absolute nominal, or exactly zero.")
        println(io, "Structural zeros use the complex entry magnitude below 1e-9 of the largest complex matrix entry. " *
            "The comparison counts exclude these zeros and unsupported or nonfinite rows.")
        println(io, "Mean shift: mean-shift interval excludes zero and absolute point shift exceeds 1%. " *
            "Mean within 1% requires the entire interval within that tolerance. Other means are inconclusive.")
        for row in summaries
            println(io, row.study, " ±", row.percent, "% ", row.nominal_id, " ", row.method,
                " | ", row.physical_class, " | N=", row.n, " E/D/U/I=",
                (row.equivalent, row.different, row.inconclusive, row.insensitive),
                "; structural=", row.structural_zeros, "; unsupported/unavailable=", (row.unsupported, row.unavailable),
                "; mean shift/within 1%/inconclusive=", (row.mean_shift, row.mean_within_one_percent, row.mean_inconclusive),
                "; warm speed-up=", format_number(row.warm_speedup),
                "; sigma: ", row.sigma_verdict, "; mean: ", row.mean_verdict)
        end
        for row in method_summaries
            row.methods_disagree > 0 || continue
            println(io, "LEP methods differ | ", (row.study, row.percent, row.nominal_id, row.physical_class),
                " | MC support measurements/central/both/neither/inconclusive=",
                (row.mc_support_measurements, row.mc_support_central, row.mc_support_both,
                    row.mc_support_neither, row.mc_support_inconclusive))
        end
        for entry in geometry
            println(io, entry.nominal_id, ": bedding ", 1000entry.bedding_thickness,
                " mm; nominal/inner/all margins ",
                1000 .* (entry.nominal_margin, entry.inner_margin, entry.all_margin), " mm")
        end
        for status in status_rows
            println(io, status.study, " ", status.method, ": ", status.status, ": ", status.message)
        end
        for timing in declaration.timings
            haskey(timing, :wall_seconds) && println(io, timing.study, " ±", timing.percent,
                "%: ", timing.wall_seconds, " s including archives.")
        end
    end
    println("wrote ", output)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("Usage: julia --project=test benchmarks/uncertainty_tables.jl <run folder>... <output folder>")
    BLAS.set_num_threads(1)
    sources = abspath.(ARGS[1:end-1])
    declarations = Dict(source => import_data(:observed, joinpath(source, "uncertainty_run.json")).gridpoint
        for source in sources)
    frequency_sources = filter(source -> get(declarations[source], :frequency_domain, false), sources)
    allunique(declarations[source].script for source in frequency_sources) ||
        error("Supply at most one frequency-domain folder per cable")
    for source in sources
        any(frequency -> frequency != source && declarations[frequency].script == declarations[source].script,
            frequency_sources) && continue
        candidates = filter(candidate -> declarations[candidate].script == declarations[source].script &&
            isfile(joinpath(candidate, "uncertainty_transient_10_monte_carlo.json")), sources)
        transient_source = get(declarations[source], :frequency_domain, false) && !isempty(candidates) ?
            first(candidates) : source
        generate_uncertainty_tables(source, abspath(last(ARGS)); transient_source)
    end
end
