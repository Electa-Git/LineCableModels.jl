# Explicit stranded core, lead sheath, unarmoured.
# Two DC poles, 150 km long.

using LineCableModels
using LineCableModels.Engine: compare
using LineCableModels.ReportBuilder: BenchmarkTableDefinition
using LinearAlgebra: diag, eigmin, Hermitian
using Dates
using Measurements, Distributions
using Statistics: mean, std

ARGS in (String[], ["time_domain"], ["pscad"], ["time_domain", "pscad"], ["model_settings"], ["modal"], ["transient_sensitivity"], ["dc_behaviour"], ["dc_behaviour", "pscad"], ["dc_behaviour", "pscad", "supplementary"]) ||
    (length(ARGS) == 4 && ARGS[1:3] == ["dc_behaviour", "pscad", "supplementary"]) ||
    (length(ARGS) >= 5 && ARGS[1:3] == ["dc_behaviour", "pscad", "runs"]) ||
    error("Usage: julia --project=test benchmarks/cable_320kv_no_armour_dc_bipole.jl [time_domain] [pscad] | model_settings | modal | transient_sensitivity | dc_behaviour [pscad [supplementary [base_folder]]] | dc_behaviour pscad runs <indices>... <base_folder>")
run_time_domain = "time_domain" in ARGS
run_pscad = "pscad" in ARGS
run_model_settings = ARGS == ["model_settings"]
run_transient_sensitivity = ARGS == ["transient_sensitivity"]
transient_started = time_ns()
run_dc_behaviour = "dc_behaviour" in ARGS
run_dc_supplementary = "supplementary" in ARGS
run_dc_selected = length(ARGS) >= 3 && ARGS[1:3] == ["dc_behaviour", "pscad", "runs"]
dc_run_indices = run_dc_selected ? parse.(Int, ARGS[4:end-1]) : Int[]
allunique(dc_run_indices) && all(i -> 1 <= i <= 55, dc_run_indices) || error("Run indices must be unique integers from 1 to 55")
dc_base_folder = run_dc_selected || length(ARGS) == 4 ? abspath(last(ARGS)) : nothing

frequencies = sort!(vcat(exp10.(range(-3, 6; length = 226)), 50.0))

"""
Construct the cable line-parameter problem from physical data in SI units and temperatures in °C.
"""
function line_parameters_problem(
        strand_diameter, core_lay_ratio, inner_semicon_thickness,
        insulation_thickness, outer_semicon_thickness, water_blocking_thickness,
        sheath_thickness, inner_sheath_thickness, jacket_thickness, core_rho, core_mu_r,
        core_temperature_coefficient, inner_semicon_rho, inner_semicon_eps_r,
        outer_semicon_rho, outer_semicon_eps_r, insulation_rho, insulation_eps_r,
        insulation_mu_r, sheath_rho, sheath_mu_r, sheath_temperature_coefficient, pe_rho,
        pe_eps_r, pp_rho, pp_eps_r, pp_mu_r, reference_temperature, cable_x_left,
        cable_x_right, cable_y, line_length, voltage, conductor_cross_section, temperature,
        earth_rho, earth_eps_r, earth_mu_r, frequencies,
        stranded_core, core_diameter,
)
    core = Material(:conductor, core_rho, 1.0, core_mu_r,
        reference_temperature, core_temperature_coefficient)
    semicon_inner = Material(:semicon, inner_semicon_rho, inner_semicon_eps_r)
    semicon_outer = Material(:semicon, outer_semicon_rho, outer_semicon_eps_r)
    xlpe = Material(:insulator, insulation_rho, insulation_eps_r, insulation_mu_r)
    sheath = Material(:conductor, sheath_rho, 1.0, sheath_mu_r,
        reference_temperature, sheath_temperature_coefficient)
    pe = Material(:insulator, pe_rho, pe_eps_r)
    pp = Material(:insulator, pp_rho, pp_eps_r, pp_mu_r)

    # The centre wire and its six courses give an envelope of 13 strand diameters.
    radius = stranded_core ? 13strand_diameter / 2 : core_diameter / 2
    conductor = if stranded_core
        stranded(core;
            shape = Disk(strand_diameter / 2), lay = LayRatio(core_lay_ratio),
            boundary = Disk(radius), fill = Material(kind = :insulator, rho = Inf))
    else
        Group(:core, Region(:core_conductor, Disk(radius), core))
    end
    parts = AbstractCablePart[conductor]
    for (tag, layer_thickness, material) in (
        (:core_semicon_inner, inner_semicon_thickness, semicon_inner),
        (:core_insulation, insulation_thickness, xlpe),
        (:core_semicon_outer, outer_semicon_thickness, semicon_outer),
        (:core_water_blocking, water_blocking_thickness, semicon_outer),
    )
        iszero(layer_thickness) && continue
        push!(parts, Region(tag, Shell(layer_thickness), material))
        radius += layer_thickness
    end
    sheath_outer = radius + sheath_thickness
    if !iszero(sheath_thickness)
        push!(parts, Group(:sheath, Region(:sheath_lead, Annulus(radius, sheath_outer), sheath)))
    end
    if !iszero(inner_sheath_thickness)
        push!(parts, Region(:sheath_inner_pe, Shell(inner_sheath_thickness), pe))
    end
    if !iszero(jacket_thickness)
        push!(parts, Region(:sheath_outer_jacket, Shell(jacket_thickness), pp))
    end

    design = build(CableDesign, "cable_320kv_no_armour", Stack(parts);
        nominal_data = (designation_code = "cable_320kv_no_armour",
            U0 = voltage / 1000,
            conductor_cross_section = conductor_cross_section / 1e-6))
    earth = homogeneous(rho = earth_rho, eps_r = earth_eps_r, mu_r = earth_mu_r)
    cable_x = (cable_x_left, cable_x_right)
    connections = [Dict(terminal => length(design.terminal_order) * (index - 1) + offset
                        for (offset, terminal) in enumerate(design.terminal_order))
                   for index in 1:2]
    system = build(LineCableSystem, fill(design, 2),
        [Pose2(cable_x[index], cable_y) for index in 1:2];
        connections, environment = earth, system_id = "cable_320kv_no_armour_dc_bipole",
        line_length = line_length)
    return LineParametersProblem(system; temperature, earth_props = earth, frequencies)
end

nominal = (;
    strand_diameter = 3.66e-3, core_lay_ratio = 15.0,
    inner_semicon_thickness = 2e-3, insulation_thickness = 20e-3,
    outer_semicon_thickness = 1.5e-3, water_blocking_thickness = 0.3e-3,
    sheath_thickness = 3.3e-3, inner_sheath_thickness = 3e-3, jacket_thickness = 10e-3,
    core_rho = 1.724e-8, core_mu_r = 1.0, core_temperature_coefficient = 0.00393,
    inner_semicon_rho = 1000.0, inner_semicon_eps_r = 1000.0, outer_semicon_rho = 500.0,
    outer_semicon_eps_r = 1000.0, insulation_rho = 1e14, insulation_eps_r = 2.3,
    insulation_mu_r = 1.0, sheath_rho = 2.14e-7, sheath_mu_r = 1.0,
    sheath_temperature_coefficient = 0.004, pe_rho = 1.97e14, pe_eps_r = 2.3, pp_rho = 1e15,
    pp_eps_r = 2.8, pp_mu_r = 1.0, reference_temperature = 20.0, cable_x_left = -0.5,
    cable_x_right = 0.5, cable_y = -1.0, line_length = 150000.0, voltage = 320e3,
    conductor_cross_section = 1600e-6, temperature = 70.0, earth_rho = 1000.0,
    earth_eps_r = 10.0, earth_mu_r = 1.0, frequencies = frequencies,
    stranded_core = true, core_diameter = 46.3e-3,
)

swept = (
    :strand_diameter, :core_lay_ratio, :inner_semicon_thickness, :insulation_thickness,
    :outer_semicon_thickness, :water_blocking_thickness, :sheath_thickness,
    :inner_sheath_thickness, :jacket_thickness, :core_rho, :core_mu_r,
    :core_temperature_coefficient, :inner_semicon_rho, :inner_semicon_eps_r,
    :outer_semicon_rho, :outer_semicon_eps_r, :insulation_rho, :insulation_eps_r,
    :insulation_mu_r, :sheath_rho, :sheath_mu_r, :sheath_temperature_coefficient, :pe_rho,
    :pe_eps_r, :pp_rho, :pp_eps_r, :pp_mu_r,
)

corrected_materials = (; core_rho = 2.3853e-8,
    core_mu_r = 1.01663, insulation_eps_r = 2.77781, insulation_mu_r = 1.59891,
    sheath_mu_r = 0.999983, pe_eps_r = 2.51862, pp_mu_r = 1.1263)

formulation = Formulation(earth_impedance = :unified, earth_admittance = :unified,
    shunt_model = :equivalent, insulation_admittance = :lossy, semicon_admittance = :lossy;
    options = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
folder = joinpath(pwd(), ".linecablemodels", "benchmarks",
    "cable_320kv_no_armour_dc_bipole", Dates.format(now(), "yyyymmdd-HHMMSS"))
mkpath(folder)

folded = (inner_semicon_thickness = 0.0, outer_semicon_thickness = 0.0, water_blocking_thickness = 0.0,
    insulation_thickness = nominal.insulation_thickness + nominal.inner_semicon_thickness +
        nominal.outer_semicon_thickness + nominal.water_blocking_thickness)
correction_steps = (;
    solid_core = (; stranded_core = false, corrected_materials[(:core_rho, :core_mu_r)]...),
    core_lay = (; stranded_core = false, corrected_materials[(:core_rho, :core_mu_r, :insulation_mu_r)]...),
    folded_layers = (; folded..., insulation_eps_r = corrected_materials.insulation_eps_r),
    pe = (; pe_eps_r = corrected_materials.pe_eps_r),
    all = (; corrected_materials..., folded..., stranded_core = false),
)
correction_status = [
    "solid_core" => "computed",
    "core_lay" => "computed",
    "folded_layers" => "computed",
    "armour_tube" => "not applicable",
    "armour_lay" => "not applicable",
    "pe" => "computed",
    "all" => "computed",
]

if ARGS == ["modal"]
    let
        started = time_ns()
        timings = NamedTuple[]
        registry = [
            ("impedance", "unified", "unified", "Unified circumferential earth impedance with complete enclosed-current normalization", "supported"),
            ("impedance", "ametani2009", "ametani2009", "Ametani mixed-pair homogeneous-earth impedance (2009)", "registered; coaxial implementation unavailable"),
            ("impedance", "carson1926", "carson1926", "Carson homogeneous-earth overhead impedance (1926)", "registered; coaxial implementation unavailable"),
            ("impedance", "default", "unified", "Unified circumferential earth impedance with complete enclosed-current normalization", "supported"),
            ("impedance", "gary1976", "gary1976", "Gary complex-depth approximation (1976)", "ArgumentError: formula :gary1976 has no expression for :self, 2, 2"),
            ("impedance", "lucca1994", "lucca1994", "Lucca mixed-pair homogeneous-earth impedance (1994)", "ArgumentError: formula :lucca1994 has no expression for :self, 2, 2"),
            ("impedance", "pollaczek1926", "pollaczek1926", "Pollaczek homogeneous-earth underground impedance (1926)", "registered; coaxial implementation unavailable"),
            ("impedance", "saad1996", "saad1996", "Saad underground closed form (1996)", "supported"),
            ("impedance", "wedepohl1973", "wedepohl1973", "Wedepohl-Wilcox low-frequency underground approximation (1973)", "supported"),
            ("impedance", "wise1934", "wise1934", "Wise homogeneous-earth overhead impedance (1934) — not yet implemented", "registered; coaxial implementation unavailable"),
            ("impedance", "xue2018", "xue2018", "Xue homogeneous-earth underground impedance (2018) — not yet implemented", "registered; coaxial implementation unavailable"),
            ("admittance", "unified", "unified", "Unified circumferential earth potential with complete enclosed-current normalization", "supported"),
            ("admittance", "default", "unified", "Unified circumferential earth potential with complete enclosed-current normalization", "supported"),
            ("admittance", "ideal", "ideal", "Ideal-earth electrostatic image potential coefficients", "supported"),
            ("admittance", "pollaczek1926", "pollaczek1926", "Pollaczek underground potential coefficients (1926) — not yet implemented", "registered; coaxial implementation unavailable"),
            ("admittance", "wise1948", "wise1948", "Wise homogeneous-earth overhead potential coefficient (1948) — not yet implemented", "registered; coaxial implementation unavailable"),
            ("admittance", "xue2018", "xue2018", "Xue homogeneous-earth underground potential coefficient (2018) — not yet implemented", "registered; coaxial implementation unavailable"),
        ]
        basis_indices = [only(findall(f -> isapprox(f, requested; rtol = 1e-12), frequencies))
            for requested in (1.0, 50.0, 1e3, 1e4, 1e5, 1e6)]
        infeasible = String[]
        problem = line_parameters_problem(nominal...)
        println("modal run: ", folder)
        flush(stdout)
        for study in ("baseline", "one_at_a_time", "temperature_earth", "layout", "corrections", "earth_formulas")
            section_start = time_ns()
            point_count = 0
            for input in (study == "one_at_a_time" ? swept : (study,))
                selected = formulation
                inputs = NamedTuple[]
                names = String[]
                if study == "baseline"
                    inputs = [nominal]
                    names = ["baseline"]
                    space = Gridspace{LineParametersProblem}(p -> line_parameters_problem(p...), (Grid(inputs),))
                elseif study == "one_at_a_time"
                    levels = nominal[input] .* range(0.9, 1.1; length = 9)
                    inputs = [merge(nominal, NamedTuple{(input,)}((level,))) for level in levels]
                    space = Gridspace{LineParametersProblem}(p -> line_parameters_problem(p...), (Grid(inputs),))
                    names = [string(input, "=", level) for level in levels]
                elseif study == "temperature_earth"
                    temperatures, resistivities = 20.0:10.0:90.0, exp10.(1:0.25:4)
                    space = Gridspace{LineParametersProblem}(
                        (temperature, rho) -> line_parameters_problem(merge(nominal, (; temperature, earth_rho = rho))...),
                        (Grid(temperatures), Grid(resistivities)))
                    inputs = [merge(nominal, (; temperature, earth_rho = rho)) for rho in resistivities for temperature in temperatures]
                    names = ["temperature=$(p.temperature);earth_rho=$(p.earth_rho)" for p in inputs]
                elseif study == "layout"
                    spacings, depths = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0), (0.8, 1.0, 1.25, 1.5, 2.0)
                    radii = [LineCableModels.support(design.geometry.outer) for design in problem.system.designs]
                    clearance = maximum(radii[1:end-1] .+ radii[2:end])
                    append!(infeasible, ["spacing=$(spacing);depth=$(depth)" for depth in depths for spacing in spacings if spacing <= clearance])
                    spacings = filter(>(clearance), spacings)
                    space = Gridspace{LineParametersProblem}((spacing, depth) -> line_parameters_problem(merge(nominal,
                        (cable_x_left = -spacing / 2, cable_x_right = spacing / 2, cable_y = -depth))...),
                        (Grid(spacings), Grid(depths)))
                    inputs = [merge(nominal, (cable_x_left = -spacing / 2, cable_x_right = spacing / 2, cable_y = -depth)) for depth in depths for spacing in spacings]
                    names = ["spacing=$(spacing);depth=$(depth)" for depth in depths for spacing in spacings]
                elseif study == "corrections"
                    inputs = [merge(nominal, changes) for changes in values(correction_steps)]
                    names = string.(collect(keys(correction_steps)))
                    space = Gridspace{LineParametersProblem}(p -> line_parameters_problem(p...), (Grid(inputs),))
                else
                    earth_pairs = ((:unified, :unified), (:saad1996, :unified),
                        (:wedepohl1973, :unified), (:unified, :ideal), (:saad1996, :ideal),
                        (:wedepohl1973, :ideal))
                    selected = Formulation(earth_impedance = Grid(first.(earth_pairs)),
                        earth_admittance = Grid(last.(earth_pairs)), combine = :zip,
                        shunt_model = :equivalent, insulation_admittance = :lossy, semicon_admittance = :lossy;
                        options = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
                    inputs = fill(nominal, length(earth_pairs))
                    names = ["impedance=$(z);admittance=$(y)" for (z, y) in earth_pairs]
                    space = Gridspace{LineParametersProblem}(p -> line_parameters_problem(p...), (Grid((nominal,)),))
                end
                if !isempty(inputs)
                    results = compute(space, selected; modal = :default)
                    for (i, result) in enumerate(results)
                        point_count += 1
                        observed = ObservedResult(PropagationParameters(result; line_length = nominal.line_length),
                            ((gamma, real), (gamma, imag), velocity, Zc);
                            gridpoint = (; id = details(result).data.gridpoint), length_unit = :base, clip = false)
                        metadata = (; study, point = names[i], inputs = inputs[i],
                            basis_indices, diagnostics = details(result).data.modal.diagnostics)
                        archive = joinpath(folder, "modal_$(study)_$(lpad(point_count, 4, '0')).json")
                        save([ObservedResult(merge(observed.gridpoint, (; modal = metadata)), observed.quantities,
                                observed.errors, observed.timings),
                            ObservedResult(result[basis_indices], (Ti, Tv);
                                gridpoint = (; id = observed.gridpoint.id), length_unit = :base, clip = false),
                            ObservedResult(transform(PhaseDomain, result[1]), (Z, Y);
                                gridpoint = (; id = observed.gridpoint.id), length_unit = :base, clip = false)],
                            archive)
                        run(`gzip -n -9 -- $archive`)
                    end
                end
                study == "one_at_a_time" && println("modal one_at_a_time: ", input, " complete")
                flush(stdout)
            end
            elapsed = (time_ns() - section_start) / 1e9
            push!(timings, (; study, wall_seconds = elapsed, points = point_count))
            println("modal ", study, ": ", elapsed, " s, ", point_count, " points")
            flush(stdout)
        end
        metadata = (; id = nothing, nominal, swept, correction_status = Tuple.(correction_status), registry, infeasible,
            terminal_order = problem.system.terminal_order, timings,
            wall_seconds = (time_ns() - started) / 1e9)
        save(ObservedResult(metadata, NamedTuple[], NamedTuple[], (;)), joinpath(folder, "modal_run.json"))
    end
else
    problem = line_parameters_problem(nominal...)
    compute_time = @elapsed result = compute(problem, formulation)
    export_data(:json, problem; file_name = joinpath(folder, "problem.json"))
    save(ObservedResult(result, (R, X, G, B); length_unit = :base, clip = false),
        joinpath(folder, "baseline.json"))
    println("cable_320kv_no_armour_dc_bipole: ", length(problem.system.terminal_order),
        " terminals, compute ", compute_time, " s, ", folder)
end

if isempty(ARGS)
# One at a time, temperature, and earth resistivity.
open(joinpath(folder, "ranking.txt"), "w") do io
    for (section, names) in ((:one_at_a_time, swept),
            (:temperature, (:temperature,)), (:earth, (:earth_rho,)))
        points = ObservedResult[]
        ranking = NamedTuple[]
        section_time = @elapsed begin
            for name in names
                levels = section === :one_at_a_time ? nominal[name] .* (0.9, 1.1) :
                    section === :temperature ? (20.0, 40.0, 90.0) : (100.0, 10000.0)
                if section === :one_at_a_time
                    run = [compute(line_parameters_problem(
                        map(k -> k === name ? level : nominal[k], keys(nominal))...),
                        formulation) for level in levels]
                else
                    grids = map(k -> k === name ? Grid(levels) : Grid((nominal[k],)), keys(nominal))
                    space = Gridspace{LineParametersProblem}(line_parameters_problem, grids)
                    run = compute(space, formulation)
                end
                completed = compare(result, run, [Z, Y]; bands = (:all, :dc, :harmonic, :narrow))
                observed = observables(run; comparisons = completed, length_unit = :base, clip = false)
                append!(points, observed)
                for (index, level) in enumerate(levels)
                    deviations = map((Y, Z)) do quantity
                        base = observe(result, quantity)
                        point = observe(run[index], quantity)
                        curve = [maximum(abs, point[:, :, k] - base[:, :, k]) /
                                 maximum(abs, diag(base[:, :, k])) for k in eachindex(frequencies)]
                        value, k = findmax(curve)
                        rms = map((:all, :dc, :harmonic, :narrow)) do band
                            indices = only(record for record in observed[index].errors
                                if record.request === quantity && record.band === band).settings.indices
                            sqrt(sum(abs2, curve[indices]) / length(indices))
                        end
                        (; rms, value, frequency = frequencies[k])
                    end
                    push!(ranking, (; name, level, deviations))
                end
            end
            artifact = report(BenchmarkTableDefinition(), points;
                reference = ObservedResult(result; length_unit = :base, clip = false))
            save(artifact, joinpath(folder, string(section, ".json")))
        end
        sort!(ranking; by = entry -> Tuple(deviation.rms[1] for deviation in entry.deviations),
            rev = true)
        printed = IOBuffer()
        println(printed, section, ": ", length(names), " swept inputs, ", length(points),
            " computed points, wall ", section_time, " s")
        for entry in ranking, (index, quantity) in enumerate((Y, Z))
            println(printed, section, " ", entry.name, "=", entry.level, " ", nameof(quantity),
                " RMS_e (all,dc,harmonic,narrow)=", entry.deviations[index].rms,
                " max_e=", entry.deviations[index].value,
                " frequency_Hz=", entry.deviations[index].frequency)
        end
        text = String(take!(printed))
        print(io, text)
        print(text)
        flush(io)
    end

    # Representation variants against the same baseline.
    names = (:conductor, :semicons, :dielectric, :shunt)
    ranking = NamedTuple[]
    section_time = @elapsed begin
        run = [
            compute(line_parameters_problem(merge(nominal, (; stranded_core = false, corrected_materials...))...), formulation),
            compute(line_parameters_problem(merge(nominal, (;
                inner_semicon_thickness = 0.0, outer_semicon_thickness = 0.0,
                water_blocking_thickness = 0.0,
                insulation_thickness = nominal.insulation_thickness + nominal.inner_semicon_thickness +
                    nominal.outer_semicon_thickness + nominal.water_blocking_thickness))...), formulation),
            compute(problem, Formulation(earth_impedance = :unified, earth_admittance = :unified,
                shunt_model = :equivalent, insulation_admittance = :lossless,
                semicon_admittance = :lossless;
                options = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false))),
            compute(problem, Formulation(earth_impedance = :unified, earth_admittance = :unified,
                shunt_model = :boundary, insulation_admittance = :lossy, semicon_admittance = :lossy;
                options = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false))),
        ]
        completed = compare(result, run, [Z, Y]; bands = (:all, :dc, :harmonic, :narrow))
        points = observables(run; comparisons = completed, length_unit = :base, clip = false)
        for (index, name) in enumerate(names)
            deviations = map((Y, Z)) do quantity
                base = observe(result, quantity)
                point = observe(run[index], quantity)
                curve = [maximum(abs, point[:, :, k] - base[:, :, k]) /
                         maximum(abs, diag(base[:, :, k])) for k in eachindex(frequencies)]
                value, k = findmax(curve)
                rms = map((:all, :dc, :harmonic, :narrow)) do band
                    indices = only(record for record in points[index].errors
                        if record.request === quantity && record.band === band).settings.indices
                    sqrt(sum(abs2, curve[indices]) / length(indices))
                end
                (; rms, value, frequency = frequencies[k])
            end
            push!(ranking, (; name, deviations))
        end
        artifact = report(BenchmarkTableDefinition(), points;
            reference = ObservedResult(result; length_unit = :base, clip = false))
        save(artifact, joinpath(folder, "representation.json"))
    end
    sort!(ranking; by = entry -> Tuple(deviation.rms[1] for deviation in entry.deviations),
        rev = true)
    printed = IOBuffer()
    println(printed, "representation: ", length(names), " variants, wall ", section_time, " s")
    for entry in ranking, (index, quantity) in enumerate((Y, Z))
        println(printed, "representation ", entry.name, " ", nameof(quantity),
            " RMS_e (all,dc,harmonic,narrow)=", entry.deviations[index].rms,
            " max_e=", entry.deviations[index].value,
            " frequency_Hz=", entry.deviations[index].frequency)
    end
    text = String(take!(printed))
    print(io, text)
    print(text)
    flush(io)
end

# Joint uncertainty: uniform support of ±10 % for the swept inputs.
grids = map(k -> k in swept ? Grid(nominal[k], (10 / sqrt(3),)) :
    Grid((nominal[k],)), keys(nominal))
space = Gridspace{LineParametersProblem}(line_parameters_problem, grids)
monte_carlo_time = @elapsed monte_carlo = compute(ParametricProblem(space),
    MonteCarlo(formulation; trials = 200, seed = 20261016, distribution = :uniform,
        options = (retain_details = true, on_error = :resample, max_failures = 20)))
linear_error_time = @elapsed linear_error = compute(ParametricProblem(space),
    LinearError(formulation))

open(joinpath(folder, "ranking.txt"), "a") do io
    printed = IOBuffer()
    summary = only(details(monte_carlo).data.failure_summary)
    println(printed, "uncertainty MonteCarlo: trials=", trial_count(monte_carlo, 1),
        " failures=", summary.failed, " wall_s=", monte_carlo_time)
    println(printed, "uncertainty LinearError: computations=1 failures=0 wall_s=", linear_error_time)
    for failure in only(details(monte_carlo).data.failures)
        println(printed, "uncertainty failed realization: ", failure)
    end
    points = ObservedResult[]
    requests = (R, X, G, B, (statistics, Z, mean), (statistics, Z, std),
        (statistics, Y, mean), (statistics, Y, std))
    for (run, elapsed) in ((monte_carlo, monte_carlo_time), (linear_error, linear_error_time))
        completed = compare(result, collect(run), [Z, Y]; bands = (:all, :dc, :harmonic, :narrow))
        append!(points, observables(run, requests; comparisons = completed,
            timings = (wall_seconds = elapsed,), length_unit = :base, clip = false))
    end
    save(points, joinpath(folder, "uncertainty.json"))
    for quantity in (Z, Y)
        base = observe(result, quantity)
        scale = [maximum(abs, diag(base[:, :, k])) for k in eachindex(frequencies)]
        spread = map((monte_carlo, linear_error)) do run
            deviation = observe(run, statistics, quantity, std, 1)
            [maximum(deviation[:, :, k]) / scale[k] for k in eachindex(frequencies)]
        end
        center = observe(monte_carlo, statistics, quantity, mean, 1)
        shift = [maximum(abs, center[:, :, k] - base[:, :, k]) / scale[k]
                 for k in eachindex(frequencies)]
        metrics = map(((:MonteCarlo_s, spread[1]), (:LinearError_s, spread[2]),
                (:MonteCarlo_m, shift))) do (label, curve)
            value, k = findmax(curve)
            rms = map((:all, :dc, :harmonic, :narrow)) do band
                indices = only(record for record in first(points).errors
                    if record.request === quantity && record.band === band).settings.indices
                sqrt(sum(abs2, curve[indices]) / length(indices))
            end
            println(printed, "uncertainty ", label, " ", nameof(quantity),
                " RMS (all,dc,harmonic,narrow)=", rms, " max=", value,
                " frequency_Hz=", frequencies[k])
            (; rms, value, frequency = frequencies[k])
        end
        println(printed, "uncertainty LEP/MC RMS_s all ", nameof(quantity), "=",
            metrics[2].rms[1] / metrics[1].rms[1])
    end
    text = String(take!(printed))
    print(io, text)
    print(text)
    flush(io)
end
end

if run_time_domain || (run_pscad && !run_dc_behaviour) || run_transient_sensitivity
    using UniversalLineModel, XLSX, TOML
    using LineCableModels: DoubleExponential
    using LineCableModels.PSCAD: PSCADTransientFormulation, LineConstantsError, RemoteConfig, identify

    let
        println("baseline: 227 samples from 1 mHz, including 50 Hz")
        section_start = time_ns()
        modal = compute(ModalAnalysisProblem(result), ModalAnalysisFormulation(:default))
        segment = PropagationParameters(modal; line_length = nominal.line_length)
        # Travel times use the modal phase velocities at the highest sampled frequency.
        travel_times = nominal.line_length ./ velocity(segment)[:, end]
        tmax = 20minimum(travel_times)
        window = (tmax = tmax, step = tmax / 4000)
        save(ObservedResult(segment; timings = (wall_seconds = (time_ns() - section_start) / 1e9,),
            length_unit = :base, clip = false), joinpath(folder, "segment.json"))
        println("modal: wall ", (time_ns() - section_start) / 1e9,
            " s, travel times ", travel_times, " s, window ", window)

        if run_time_domain
            fitting_problem = LineModelProblem(segment; line = problem)
            section_start = time_ns()
            sampling_frequencies = sort!(vcat(exp10.(range(0, 6; length = 101)), 50.0))
            sampling_problem = line_parameters_problem(merge(nominal,
                (; frequencies = sampling_frequencies))...)
            sampling_result = compute(sampling_problem, formulation)
            sampling_modal = compute(ModalAnalysisProblem(sampling_result), ModalAnalysisFormulation(:default))
            sampling_segment = PropagationParameters(sampling_modal; line_length = nominal.line_length)
            sampling_fitting_problem = LineModelProblem(sampling_segment; line = sampling_problem)
            save(ObservedResult(sampling_result; clip = false), joinpath(folder, "coarse_baseline.json"))
            println("coarse: 101 samples from 1 Hz, plus 50 Hz; wall ",
                (time_ns() - section_start) / 1e9, " s")

            formulations = LineModelFormulation(:ulm; options = (
                yc = (npoles = Grid((8, 16)),),
                h = (npoles = Grid((8, 16)), td = Grid((false, true))),
                passivity = Grid((false, true)),
                dc = Grid(((;), (method = :add_pole,), (method = :functional_form,))),
                line = formulation))
            fitted_models = LineModel[]
            fit_points = ObservedResult[]
            fit_times = Float64[]
            for (sampling, selected_problem, selections) in (
                    (:baseline, fitting_problem, Iterators.flatten((formulations,
                        (LineModelFormulation(:fdcm), LineModelFormulation(:marti))))),
                    (:coarse, sampling_fitting_problem, Iterators.flatten((
                        LineModelFormulation(:ulm; options = (yc = (npoles = 16,),
                            h = (npoles = 16, td = true), passivity = true,
                            dc = Grid(((;), (method = :add_pole,), (method = :functional_form,))),
                            line = formulation)), (LineModelFormulation(:marti),)))))
                section_start = time_ns()
                for selected in selections
                    elapsed = @elapsed fitted = compute(selected_problem, selected; options = (timing = true,))
                    push!(fitted_models, fitted)
                    push!(fit_times, elapsed)
                    diagnostic = details(fitted).data
                    evidence = (; method = formula_id(selected.formula),
                        options = (; (key => value for (key, value) in pairs(diagnostic.options)
                            if key !== :line)...), line = sampling,
                        samples = sampling === :baseline ? length(frequencies) : length(sampling_frequencies),
                        orders = length.(poles(fitted)), passivity = passivity(fitted),
                        fitting_error = fitting_error(fitted), dc = fitted.dc,
                        step_check = get(diagnostic, :step_check, false),
                        step_err = get(diagnostic, :step_err, NaN),
                        step_err_late = get(diagnostic, :step_err_late, NaN),
                        delays = delays(fitted), travel_times, window, wall_seconds = elapsed,
                        action_timing = diagnostic.timing)
                    println("time domain fits compute ", length(fitted_models), ": wall ", elapsed, " s")
                    flush(stdout)
                    point = ObservedResult(fitted; timings = (wall_seconds = elapsed,
                        action = diagnostic.timing), clip = false)
                    # Keep fit diagnostics as ordinary metadata in the portable archive.
                    push!(fit_points, ObservedResult(merge(point.gridpoint, (study = evidence,)),
                        NamedTuple[], point.errors, point.timings))
                    for quantity in (:admittance, :propagation)
                        observed = ObservedResult(fitted, ((fitting_error, quantity),);
                            gridpoint = fit_points[end].gridpoint, clip = false)
                        export_data(:xlsx, observed;
                            file_name = joinpath(folder, "line_model_$(length(fitted_models))_$(quantity).xlsx"))
                    end
                    println("time domain fits ", length(fitted_models), " ", evidence.method,
                        ": orders ", evidence.orders, ", wall ", elapsed, " s")
                    flush(stdout)
                end
                save(fit_points, joinpath(folder, "line_models.json"))
                println("time domain fits ", sampling, ": wall ", (time_ns() - section_start) / 1e9, " s")
            end
        end

        section_start = time_ns()
        labels = terminal_labels(problem.system)
        core1 = first(label for label in labels if endswith(label, ":1:core"))
        # Zip the sheath grids so both bonds have the same R at each reference point.
        bonds = [TerminalImpedance(label; R = Grid((1.0, 10.0, 20.0)))
            for label in labels if endswith(label, ":sheath")]
        reference_problems = TransientAnalysisProblem[TransientAnalysisProblem(problem;
            sending = (Step(core1, nominal.voltage), assignments...), receiving = (), window)
            for assignments in zip(bonds...)]
        both_end_bonds = Tuple(TerminalImpedance(label; R = 10.0)
            for label in labels if endswith(label, ":sheath"))
        push!(reference_problems, TransientAnalysisProblem(problem;
            sending = (Step(core1, nominal.voltage), both_end_bonds...),
            receiving = both_end_bonds, window))
        alpha, beta = 1.47e4, 2.47e6
        pulse_peak_time = log(beta / alpha) / (beta - alpha)
        pulse_amplitude = nominal.voltage / (exp(-alpha * pulse_peak_time) - exp(-beta * pulse_peak_time))
        pulse_window = (tmax = 3e-3, step = 0.25e-6)
        push!(reference_problems, TransientAnalysisProblem(problem;
            sending = (DoubleExponential(core1, pulse_amplitude, alpha, beta), both_end_bonds...),
            receiving = (), window = pulse_window))
        if run_transient_sensitivity
            println("transient sensitivity run: ", folder)
            flush(stdout)
            observations = ObservedResult[]
            timings = NamedTuple[]
            retained = NamedTuple[]
            passivity_grids = NamedTuple[]
            for (scenario, reference_index) in ((:step_single_point, 2), (:impulse, 5))
                scenario_start = time_ns()
                reference_problem = reference_problems[reference_index]
                N = scenario === :impulse ? 8192 : 4096
                selected = TransientFormulation(:nlt; options = (line = formulation, N, Tfactor = 1.25))
                controls = selected.formula.options.data
                T = controls.Tfactor * reference_problem.window.tmax
                band = (lowest_Hz = 1 / (2T), highest_Hz = (2controls.N - 1) / (2T))
                count = ceil(Int, 400log10(band.highest_Hz / band.lowest_Hz)) + 1
                passivity_frequencies = exp10.(range(log10(band.lowest_Hz), log10(band.highest_Hz); length = count))
                passivity_frequencies[1], passivity_frequencies[end] = band.lowest_Hz, band.highest_Hz
                nlt_settings = (; N = controls.N, T, Tfactor = controls.Tfactor, taper = controls.taper)
                push!(passivity_grids, (; scenario, band_source = "nlt_settings", band,
                    points_per_decade = 400, frequencies_Hz = passivity_frequencies, nlt = nlt_settings))
                nominal_response = nothing
                nominal_peak = NaN
                extremes = NamedTuple[]
                for study in (:baseline, :one_at_a_time, :temperature, :earth, :corrections)
                    for input in (study === :one_at_a_time ? swept : (study,))
                        inputs = if study === :baseline
                            [nominal]
                        elseif study === :corrections
                            [merge(nominal, changes) for changes in values(correction_steps)]
                        else
                            key = study === :temperature ? :temperature : study === :earth ? :earth_rho : input
                            levels = study === :temperature ? (20.0, 40.0, 70.0, 90.0) :
                                study === :earth ? (100.0, 1000.0, 10000.0) : nominal[key] .* (0.9, 1.0, 1.1)
                            [merge(nominal, NamedTuple{(key,)}((level,))) for level in levels]
                        end
                        isempty(inputs) && continue
                        group_start = time_ns()
                        space = Gridspace{LineParametersProblem}(p -> line_parameters_problem(p...), (Grid(inputs),))
                        passivity_space = Gridspace{LineParametersProblem}(line ->
                            LineParametersProblem(line.system; temperature = line.temperature,
                                earth_props = line.earth_props, frequencies = passivity_frequencies, sigma = 0.0), (space,))
                        passivity_elapsed = @elapsed scans = compute(passivity_space, formulation)
                        checks = NamedTuple[]
                        curves = NamedTuple[]
                        for scan in scans
                            z, y = observe(scan, Z), observe(scan, Y)
                            z_min = [eigmin(Hermitian((z[:, :, k] + z[:, :, k]') / 2)) for k in eachindex(passivity_frequencies)]
                            y_min = [eigmin(Hermitian((y[:, :, k] + y[:, :, k]') / 2)) for k in eachindex(passivity_frequencies)]
                            onset = findfirst(k -> z_min[k] < 0 || y_min[k] < 0, eachindex(passivity_frequencies))
                            offending_quantity = onset === nothing ? nothing : join(
                                (string(quantity) for (quantity, values) in ((:Z, z_min), (:Y, y_min)) if values[onset] < 0), ',')
                            check = (; status = onset === nothing ? "passive" : "non_passive",
                                onset_frequency_Hz = onset === nothing ? nothing : passivity_frequencies[onset],
                                onset_lower_Hz = onset === nothing || onset == 1 ? nothing : passivity_frequencies[onset - 1],
                                offending_quantity,
                                minimum_Z_eigenvalue_ohm_per_m = minimum(z_min),
                                minimum_Z_frequency_Hz = passivity_frequencies[argmin(z_min)],
                                minimum_Y_eigenvalue_S_per_m = minimum(y_min),
                                minimum_Y_frequency_Hz = passivity_frequencies[argmin(y_min)])
                            push!(checks, check)
                            push!(curves, (; id = details(scan).data.gridpoint, check,
                                Z_min_eigenvalues_ohm_per_m = z_min, Y_min_eigenvalues_S_per_m = y_min))
                        end
                        passivity_file = "transient_$(scenario)_$(study)_$(input)_passivity.json"
                        save(ObservedResult((; id = nothing, scenario, study, input,
                            band_source = "nlt_settings", band, nlt = nlt_settings,
                            points_per_decade = 400, frequencies_Hz = passivity_frequencies, curves),
                            NamedTuple[], NamedTuple[], (; wall_seconds = passivity_elapsed)), joinpath(folder, passivity_file))
                        computable = findall(check -> check.status == "passive", checks)
                        responses = TransientResponse[]
                        transient_elapsed = 0.0
                        if !isempty(computable)
                            accepted_space = Gridspace{LineParametersProblem}(p -> line_parameters_problem(p...), (Grid(inputs[computable]),))
                            problems = TransientAnalysisProblem(accepted_space; sending = reference_problem.sending,
                                receiving = reference_problem.receiving, window = reference_problem.window)
                            transient_elapsed = @elapsed responses = compute(problems, selected; options = (timing = true,))
                        end
                        if study === :baseline && !isempty(responses)
                            nominal_response = only(responses)
                            core_index = only(findall(==(core1), nominal_response.coordinates))
                            nominal_peak = peak(nominal_response, voltage).value[2, core_index]
                        end
                        group_points = ObservedResult[]
                        for index in eachindex(inputs)
                            input_name = study === :earth ? :earth_rho : input
                            level = study === :baseline ? "nominal" : study === :corrections ?
                                string(collect(keys(correction_steps))[index]) : inputs[index][input_name]
                            passivity = merge(checks[index], (; archive = passivity_file, point = index))
                            evidence = (; scenario, study, input = input_name, level,
                                status = passivity.status == "passive" ? "computable" : "not_computable",
                                passivity, window = reference_problem.window, nlt = nlt_settings)
                            response_index = findfirst(==(index), computable)
                            id = details(scans[index]).data.gridpoint
                            if response_index !== nothing
                                response = responses[response_index]
                                core_index = only(findall(==(core1), response.coordinates))
                                voltage_peak = peak(response, voltage)
                                current_peak = peak(response, current)
                                metal_indices = findall(label -> endswith(label, ":sheath") ||
                                    endswith(label, ":armor"), response.coordinates)
                                metal_peaks = [(; line_end, terminal = response.coordinates[terminal],
                                    peak_kV = voltage_peak.value[side, terminal] / 1e3)
                                    for (side, line_end) in enumerate((:sending, :receiving)) for terminal in metal_indices]
                                deviation = nominal_response === nothing ? nothing :
                                    compare(nominal_response, response, voltage).absolute[2, core_index] / nominal_peak
                                evidence = merge(evidence, (;
                                    receiving_peak_kV = voltage_peak.value[2, core_index] / 1e3,
                                    receiving_peak_time_us = voltage_peak.time[2, core_index] * 1e6,
                                    arrival_time_us = arrival_time(response, voltage)[2, core_index] * 1e6,
                                    sending_current_peak_kA = current_peak.value[1, core_index] / 1e3,
                                    receiving_scaled_rms = deviation, metal_peaks, nlt = details(response).data.nlt))
                                id = details(response).data.gridpoint
                            else
                                println("transient ", scenario, " ", input_name, " ", level,
                                    " not computable: ", passivity.offending_quantity, " onset ", passivity.onset_frequency_Hz, " Hz")
                            end
                            point = ObservedResult((; id, sensitivity = evidence), NamedTuple[], NamedTuple[],
                                (; passivity_wall_seconds = passivity_elapsed, transient_wall_seconds = transient_elapsed))
                            push!(group_points, point)
                            push!(observations, point)
                        end
                        save(observations, joinpath(folder, "transient_sensitivity.json"))
                        if study === :baseline && nominal_response !== nothing
                            save(ObservedResult(nominal_response; gridpoint = only(group_points).gridpoint, clip = false),
                                joinpath(folder, "transient_$(scenario)_nominal.json"))
                            println("transient ", scenario, " nominal receiving peak: ", nominal_peak / 1e3, " kV")
                        elseif study in (:one_at_a_time, :temperature, :earth) && !isempty(computable) && isfinite(nominal_peak)
                            change_kV = maximum(abs(group_points[index].gridpoint.sensitivity.receiving_peak_kV - nominal_peak / 1e3)
                                for index in computable)
                            endpoints = NamedTuple[]
                            for (bound, index) in ((:lower, 1), (:upper, length(inputs)))
                                response_index = findfirst(==(index), computable)
                                response_index === nothing && continue
                                push!(endpoints, (; bound, response = responses[response_index], point = group_points[index]))
                            end
                            push!(extremes, (; input = group_points[1].gridpoint.sensitivity.input, change_kV, endpoints))
                            sort!(extremes; by = entry -> entry.change_kV, rev = true)
                            length(extremes) > 5 && pop!(extremes)
                        end
                        elapsed = (time_ns() - group_start) / 1e9
                        push!(timings, (; scenario, study, input, points = length(inputs), computed_points = length(responses),
                            passivity_wall_seconds = passivity_elapsed, transient_wall_seconds = transient_elapsed, wall_seconds = elapsed))
                        println("transient ", scenario, " ", study, " ", input, ": ", elapsed, " s")
                        flush(stdout)
                    end
                end
                for entry in extremes, endpoint in entry.endpoints
                    filename = "transient_$(scenario)_$(entry.input)_$(endpoint.bound).json"
                    save(ObservedResult(endpoint.response; gridpoint = endpoint.point.gridpoint, clip = false), joinpath(folder, filename))
                    push!(retained, (; scenario, input = entry.input, bound = endpoint.bound, change_kV = entry.change_kV, filename))
                end
                println("transient ", scenario, " complete: ", (time_ns() - scenario_start) / 1e9, " s")
                flush(stdout)
            end
            metadata = (; id = nothing, nominal, swept, correction_steps, terminals = labels, passivity_grids,
                correction_status = Tuple.(correction_status), travel_times, timings, retained,
                wall_seconds = (time_ns() - transient_started) / 1e9)
            save(ObservedResult(metadata, NamedTuple[], NamedTuple[], (;)), joinpath(folder, "transient_sensitivity_run.json"))
            println("transient sensitivity complete: ", metadata.wall_seconds, " s, ", folder)
            exit()
        end

        references = Vector{TransientResponse}(undef, length(reference_problems))
        reference_points = ObservedResult[]
        for (reference_index, reference_problem) in enumerate(reference_problems)
            run_time_domain || reference_index in (2, 5) || continue
            resistance = (1.0, 10.0, 20.0, 10.0, 10.0)[reference_index]
            scenario = reference_index <= 3 ? :step_single_point : reference_index == 4 ? :step_both_ends : :impulse
            N = scenario === :impulse ? 8192 : 4096
            elapsed = @elapsed reference = compute(reference_problem,
                TransientFormulation(:nlt; options = (line = formulation, N, Tfactor = 1.25));
                options = (timing = true,))
            references[reference_index] = reference
            point = ObservedResult(reference; timings = (wall_seconds = elapsed,), clip = false)
            source = first(reference_problem.sending)
            nlt = details(reference).data.nlt
            evidence = (; scenario, resistance, N = nlt.N, T = nlt.T,
                highest_frequency = (2nlt.N - 1) / (2nlt.T),
                window = reference_problem.window, travel_times,
                sending = (source = core1, waveform = nameof(typeof(source)), amplitude = source.amplitude,
                    alpha = scenario === :impulse ? alpha : nothing, beta = scenario === :impulse ? beta : nothing,
                    sheaths = [bond.terminal for bond in reference_problem.sending[2:end]]),
                receiving = [bond.terminal for bond in reference_problem.receiving],
                pulse_peak_time = scenario === :impulse ? pulse_peak_time : nothing,
                before_first_arrival = reference_problem.window.tmax < minimum(travel_times))
            push!(reference_points, ObservedResult(merge(point.gridpoint, (study = evidence,)),
                point.quantities, point.errors, point.timings))
            save(reference_points, joinpath(folder, "references.json"))
            println("reference ", scenario, " R=", resistance, ": wall ", elapsed,
                " s, N=", nlt.N, " T=", nlt.T, " s, highest frequency=", evidence.highest_frequency, " Hz")
            flush(stdout)
        end
        println("references: wall ", (time_ns() - section_start) / 1e9, " s")

        if run_time_domain
            model_deviations = fill(NaN, 3, length(fitted_models))
            for (scenario_index, reference_index) in enumerate((2, 4, 5))
                section_start = time_ns()
                scenario = (:step_single_point, :step_both_ends, :impulse)[scenario_index]
                reference = references[reference_index]
                sending = reference_problems[reference_index].sending
                receiving = reference_problems[reference_index].receiving
                window = reference_problems[reference_index].window
                tmax = window.tmax
                core_index = findfirst(==(core1), reference.coordinates)
                growth_factor = 1.1
                for (index, fitted) in enumerate(fitted_models)
                    scenario !== :step_single_point && fit_points[index].gridpoint.study.line === :coarse && continue
                    response_file = "$(scenario)_response_$(index).json"
                    elapsed = @elapsed response = compute(TransientAnalysisProblem(fitted;
                        sending, receiving, window), TransientFormulation(:td);
                        options = (timing = true,))
                    save(ObservedResult(response; timings = (wall_seconds = elapsed,), clip = false),
                        joinpath(folder, response_file))
                    errors = map((voltage, current)) do selector
                        deviation = compare(reference, response, selector)
                        (; quantity = nameof(selector), absolute = deviation.absolute,
                            relative = deviation.relative, settings = details(deviation).data)
                    end
                    kpis = map((voltage, current)) do selector
                        extrema = peak(response, selector)
                        arrival = arrival_time(response, selector)
                        decrement, damping_error = try
                            (damping(response, selector)[:, core_index], nothing)
                        catch exception
                            exception isa DomainError || rethrow()
                            (fill(missing, 2), sprint(showerror, exception))
                        end
                        expected = selector === voltage && scenario !== :impulse ? nominal.voltage : 0.0
                        offset = dc_offset(response, selector; expected)
                        (; quantity = nameof(selector), peak = extrema.value[:, core_index],
                            peak_time = extrema.time[:, core_index], arrival_time = arrival[:, core_index],
                            damping = decrement, damping_error, dc_offset = offset[:, core_index],
                            dc_expected = expected)
                    end
                    # Compare per-channel |v| peaks in [0.8, 0.9) and [0.9, 1.0] of tmax.
                    last_windows = (findall(t -> 0.8tmax <= t < 0.9tmax, response.t),
                        findall(t -> 0.9tmax <= t <= tmax, response.t))
                    window_peaks = map((model = response, reference = reference)) do run
                        [dropdims(maximum(abs, run.voltage[:, :, indices]; dims = 3);
                            dims = 3) for indices in last_windows]
                    end
                    growth_ratios = map(window_peaks) do peaks
                        map(peaks[1], peaks[2]) do previous_peak, last_peak
                            # A zero signal in both windows has unit growth.
                            # A rise from zero has infinite growth.
                            iszero(previous_peak) ? (iszero(last_peak) ? 1.0 : Inf) :
                                last_peak / previous_peak
                        end
                    end
                    growth = (growth_ratios.model .> 1) .& (growth_ratios.model .> growth_factor .* growth_ratios.reference)
                    scaled_rms = errors[1].absolute[2, core_index] /
                        maximum(abs, reference.voltage[2, core_index, :])
                    model_deviations[scenario_index, index] = scaled_rms
                    evidence = (; scenario, fit = index, comparisons = errors, kpis, scaled_rms,
                        reference_peak = maximum(abs, reference.voltage[2, core_index, :]),
                        before_first_arrival = tmax < minimum(travel_times),
                        scale = :reference_receiving_core1_peak, growth, growth_ratios, growth_factor,
                        voltage_window_peaks = window_peaks, stability_flag = any(growth),
                        stability_windows = ((0.8tmax, 0.9tmax), (0.9tmax, tmax)))
                    point = ObservedResult(response; timings = (wall_seconds = elapsed,), clip = false)
                    point = ObservedResult(merge(point.gridpoint, (study = evidence,)),
                        point.quantities, point.errors, point.timings)
                    save(point, joinpath(folder, response_file))
                    fit_evidence = fit_points[index].gridpoint.study
                    println("model ", scenario, " ", index, " ", fit_evidence.method,
                        ": Yc_RMS=", maximum(fitting_error(fitted, :admittance)),
                        " H_RMS=", maximum(fitting_error(fitted, :propagation)),
                        " passive=", passivity(fitted).passive, " step_err=", fit_evidence.step_err,
                        " scaled_RMS=", scaled_rms, " growth=", any(growth),
                        " fit_s=", fit_times[index], " td_s=", elapsed,
                        " damping_errors=", [entry.damping_error for entry in kpis])
                    for (side, line_end) in enumerate(("sending", "receiving")),
                            (terminal, coordinate) in enumerate(response.coordinates)
                        println("growth ", index, " ", line_end, " ", coordinate,
                            ": model=", growth_ratios.model[side, terminal],
                            " reference=", growth_ratios.reference[side, terminal],
                            " flag=", growth[side, terminal])
                    end
                    flush(stdout)
                end
                println("models ", scenario, ": wall ", (time_ns() - section_start) / 1e9, " s, ", folder)
            end

            section_start = time_ns()
            sampling_points = ObservedResult[]
            for index in 51:length(fitted_models)
                coarse = fit_points[index].gridpoint.study
                baseline_index = findfirst(fit_points[1:50]) do point
                    full = point.gridpoint.study
                    full.method === coarse.method && (coarse.method === :marti ||
                        (full.options.yc.npoles == 16 && full.options.h.npoles == 16 &&
                         full.options.h.td && full.options.passivity && full.options.dc == coarse.options.dc))
                end
                full = fit_points[baseline_index].gridpoint.study
                evidence = (; coarse_fit = index, baseline_fit = baseline_index,
                    samples = (baseline = length(frequencies), coarse = length(sampling_frequencies)),
                    frequency_bounds = (baseline = extrema(frequencies), coarse = extrema(sampling_frequencies)),
                    fitting_error = (baseline = full.fitting_error, coarse = coarse.fitting_error),
                    scaled_rms = (baseline = model_deviations[1, baseline_index], coarse = model_deviations[1, index]))
                point = fit_points[index]
                push!(sampling_points, ObservedResult(merge(point.gridpoint, (; sampling = evidence)),
                    point.quantities, point.errors, point.timings))
                println("sampling ", coarse.method, " ", get(coarse.options, :dc, (;)),
                    ": fitting_error=", evidence.fitting_error, " scaled_RMS=", evidence.scaled_rms)
            end
            save(sampling_points, joinpath(folder, "sampling_comparison.json"))
            println("sampling comparison: wall ", (time_ns() - section_start) / 1e9, " s")
        end

        if run_pscad
            section_start = time_ns()
            path = ENV["LINECABLEMODELS_PSCAD_CONFIG"]
            options = (; remote = RemoteConfig(path))
            station = identify(options.remote)
            open(joinpath(folder, "pscad_station.toml"), "w") do io
                TOML.print(io, station; sorted = true)
            end
            native_folder = mkpath(joinpath(options.remote.local_root, "cable_320kv_no_armour_dc_bipole_$(basename(folder))"))
            for (scenario, reference_index) in ((:step_single_point, 2), (:impulse, 5))
                reference = references[reference_index]
                reference_problem = reference_problems[reference_index]
                core_index = findfirst(==(core1), reference.coordinates)
                pscad_formulations = PSCADTransientFormulation[PSCADTransientFormulation(model = :fd_phase; line = formulation,
                    fit = (; YMaxP = order, AMaxP = order, DCenab = dc, CPASS = passivity_mode))
                    for order in (10, 20, 30) for dc in (0, 1) for passivity_mode in (0, 1)]
                push!(pscad_formulations, PSCADTransientFormulation(model = :fd_mode, line = formulation))
                for (index, selected) in enumerate(pscad_formulations)
                    output_folder = mkpath(joinpath(folder, "pscad_$(scenario)_$(index)"))
                    before = Set(readdir(native_folder))
                    started = time_ns()
                    # Keep construction and comparison outside the station-error catch.
                    response = try
                        compute(reference_problem, selected;
                            options = merge(options, (; work_root = native_folder)))
                    catch exception
                        failure = if exception isa LineConstantsError
                            Dict("status" => "line_constants_error", "message" => exception.message,
                                "log_path" => exception.log_path)
                        elseif exception isa ErrorException
                            Dict("status" => "known_station_failure", "message" => sprint(showerror, exception))
                        else
                            rethrow()
                        end
                        open(joinpath(output_folder, "status.toml"), "w") do io
                            TOML.print(io, merge(failure, Dict("scenario" => string(scenario),
                                "model" => string(selected.model), "fit" => Dict(string(k) => v for (k, v) in pairs(selected.fit)),
                                "station" => station, "wall_seconds" => (time_ns() - started) / 1e9)); sorted = true)
                        end
                        nothing
                    end
                    # Archive the native project, logs and outputs for successful and failed runs.
                    for name in setdiff(Set(readdir(native_folder)), before)
                        cp(joinpath(native_folder, name), joinpath(output_folder, name))
                    end
                    response === nothing && continue
                    deviation = compare(reference, response, voltage)
                    reference_peak = maximum(abs, reference.voltage[2, core_index, :])
                    scaled_rms = deviation.absolute[2, core_index] / reference_peak
                    evidence = (; scenario, model = selected.model, fit = selected.fit,
                        station = details(response).data.execution.solver_identity,
                        scaled_rms, reference_peak, scale = :reference_receiving_core1_peak,
                        absolute = deviation.absolute, relative = deviation.relative,
                        before_first_arrival = reference_problem.window.tmax < minimum(travel_times))
                    point = ObservedResult(response;
                        timings = (wall_seconds = (time_ns() - started) / 1e9,), clip = false)
                    save(ObservedResult(merge(point.gridpoint, (; pscad = evidence)),
                        point.quantities, point.errors, point.timings), joinpath(output_folder, "response.json"))
                    open(joinpath(output_folder, "status.toml"), "w") do io
                        TOML.print(io, Dict("status" => "completed", "scenario" => string(scenario),
                            "model" => string(selected.model), "fit" => Dict(string(k) => v for (k, v) in pairs(selected.fit)),
                            "station" => evidence.station, "scaled_rms" => scaled_rms,
                            "reference_peak" => reference_peak,
                            "wall_seconds" => (time_ns() - started) / 1e9); sorted = true)
                    end
                end
            end
            println("pscad: wall ", (time_ns() - section_start) / 1e9, " s, ", folder)
        end
    end
end

if run_model_settings
    using UniversalLineModel

    let
        println("baseline: 227 samples from 1 mHz, including 50 Hz")
        section_start = time_ns()
        modal = compute(ModalAnalysisProblem(result), ModalAnalysisFormulation(:default))
        segment = PropagationParameters(modal; line_length = nominal.line_length)
        travel_times = nominal.line_length ./ velocity(segment)[:, end]
        tmax = 20minimum(travel_times)
        window = (tmax = tmax, step = tmax / 4000)
        fitting_problem = LineModelProblem(segment; line = problem)
        sampling_frequencies = sort!(vcat(exp10.(range(0, 6; length = 101)), 50.0))
        sampling_problem = line_parameters_problem(merge(nominal,
            (; frequencies = sampling_frequencies))...)
        sampling_result = compute(sampling_problem, formulation)
        sampling_modal = compute(ModalAnalysisProblem(sampling_result), ModalAnalysisFormulation(:default))
        sampling_segment = PropagationParameters(sampling_modal; line_length = nominal.line_length)
        sampling_fitting_problem = LineModelProblem(sampling_segment; line = sampling_problem)
        save(ObservedResult(sampling_result; clip = false),
            joinpath(folder, "coarse_baseline.json"))
        println("coarse: 101 samples from 1 Hz, plus 50 Hz")
        labels = terminal_labels(problem.system)
        core1 = first(label for label in labels if endswith(label, ":1:core"))
        bonds = Tuple(TerminalImpedance(label; R = 10.0)
            for label in labels if endswith(label, ":sheath"))
        sending = (Step(core1, nominal.voltage), bonds...)
        receiving = ()
        reference_problem = TransientAnalysisProblem(problem; sending, receiving, window)
        reference_time = @elapsed reference = compute(reference_problem,
            TransientFormulation(:nlt; options = (line = formulation, N = 4096, Tfactor = 1.25));
            options = (timing = true,))
        save(ObservedResult(reference; timings = (wall_seconds = reference_time,), clip = false),
            joinpath(folder, "reference.json"))
        core_index = findfirst(==(core1), reference.coordinates)
        reference_peak = maximum(abs, reference.voltage[2, core_index, :])
        println("model settings reference: wall ", reference_time, " s, peak ", reference_peak, " V")
        flush(stdout)

        fit_points = ObservedResult[]
        open(joinpath(folder, "model_settings.txt"), "w") do io
            for (method, selected_problem, controls) in (
                    (:marti, fitting_problem, (nothing, 2000.0)),
                    (:ulm, sampling_fitting_problem, (Inf, 100.0)))
                for control in controls
                    selected = method === :marti ? LineModelFormulation(:marti;
                        options = (transformation = (frequency = control,),)) :
                        LineModelFormulation(:ulm; options = (yc = (npoles = 16,),
                            h = (npoles = 16, td = true, max_ratio = control),
                            passivity = true, dc = (;), line = formulation))
                    fit_time = @elapsed fitted = compute(selected_problem, selected;
                        options = (timing = true,))
                    diagnostic = details(fitted).data
                    # H ratios use the largest matrix entry of each residue divided by |pole|.
                    ratios = method === :ulm ?
                        UniversalLineModel.residue_pole_ratios(UniversalLineModel.ULMModel(fitted)) : nothing
                    admittance_ratio = method === :ulm ? maximum(
                        maximum(abs, residue) / abs(pole) for (residue, pole) in
                        zip(eachslice(fitted.model.admittance.residues; dims = 3),
                            fitted.model.admittance.poles)) : nothing
                    evidence = (; method, samples = length(fitted.frequencies),
                        transformation_requested = method === :marti ? control : nothing,
                        transformation_frequency = method === :marti ? fitted.model.transformation_frequency : nothing,
                        modal = method === :marti ? diagnostic.marti : nothing,
                        max_ratio = method === :ulm ? control : nothing,
                        residue_pole_ratios = ratios,
                        largest_residue_pole_ratio = method === :ulm ? maximum(maximum, ratios) : nothing,
                        largest_admittance_residue_pole_ratio = admittance_ratio,
                        fitting_error = fitting_error(fitted), passivity = passivity(fitted),
                        orders = length.(poles(fitted)), dc = fitted.dc,
                        options = (; (key => value for (key, value) in pairs(diagnostic.options)
                            if key !== :line)...),
                        travel_times, window, wall_seconds = fit_time)
                    point = ObservedResult(fitted; timings = (wall_seconds = fit_time,), clip = false)
                    push!(fit_points, ObservedResult(merge(point.gridpoint, (; model_settings = evidence)),
                        NamedTuple[], point.errors, point.timings))
                    save(fit_points, joinpath(folder, "line_models.json"))
                    println("model settings fit ", method, " control=", control, ": wall ", fit_time, " s")
                    flush(stdout)
                    for interp in (method === :marti ? (:single,) : (:single, :double))
                        elapsed = @elapsed response = compute(TransientAnalysisProblem(fitted;
                            sending, receiving, window), TransientFormulation(:td; options = (; interp));
                            options = (timing = true,))
                        deviation = compare(reference, response, voltage)
                        scaled_rms = deviation.absolute[2, core_index] / reference_peak
                        last_windows = (findall(t -> 0.8tmax <= t < 0.9tmax, response.t),
                            findall(t -> 0.9tmax <= t <= tmax, response.t))
                        window_peaks = map((model = response, reference = reference)) do run
                            [dropdims(maximum(abs, run.voltage[:, :, indices]; dims = 3);
                                dims = 3) for indices in last_windows]
                        end
                        growth_ratios = map(window_peaks) do peaks
                            map(peaks[1], peaks[2]) do previous_peak, last_peak
                                iszero(previous_peak) ? (iszero(last_peak) ? 1.0 : Inf) :
                                    last_peak / previous_peak
                            end
                        end
                        growth_factor = 1.1
                        growth = (growth_ratios.model .> 1) .&
                            (growth_ratios.model .> growth_factor .* growth_ratios.reference)
                        transient = (; scenario = :step_single_point, resistance = 10.0, interp, scaled_rms,
                            reference_peak, scale = :reference_receiving_core1_peak,
                            absolute = deviation.absolute, relative = deviation.relative,
                            growth, growth_ratios, growth_factor, voltage_window_peaks = window_peaks,
                            stability_flag = any(growth),
                            stability_windows = ((0.8tmax, 0.9tmax), (0.9tmax, tmax)))
                        point = ObservedResult(response; timings = (wall_seconds = elapsed,), clip = false)
                        save(ObservedResult(merge(point.gridpoint, (; model_settings = evidence, transient)),
                            point.quantities, point.errors, point.timings),
                            joinpath(folder, "step_single_point_response_$(length(fit_points))_$(interp).json"))
                        text = string(method, " control=", control, " interp=", interp,
                            " frequency_Hz=", evidence.transformation_frequency,
                            " modal=", evidence.modal, " fitting_error=", evidence.fitting_error,
                            " passive=", evidence.passivity.passive,
                            " largest_residue_pole_ratio=", evidence.largest_residue_pole_ratio,
                            " scaled_RMS=", scaled_rms, " growth=", any(growth),
                            " fit_s=", fit_time, " td_s=", elapsed)
                        println(io, text)
                        println(text)
                        flush(io)
                        flush(stdout)
                    end
                end
            end
            elapsed = (time_ns() - section_start) / 1e9
            println(io, "model settings wall_seconds=", elapsed, " baseline_compute_seconds=", compute_time)
            println("model settings: wall ", elapsed, " s, ", folder)
        end
    end
end

if run_dc_behaviour
    using UniversalLineModel, TOML, EzXML
    using LinearAlgebra: Symmetric, eigvals, I
    using LineCableModels.PSCAD: PSCADTransientFormulation, LineConstantsError, NonfiniteOutputError, RemoteConfig, identify

    let
        study_start = time_ns()
        labels = terminal_labels(problem.system)
        core1 = first(label for label in labels if endswith(label, ":1:core"))
        core_index = only(findall(==(core1), labels))
        G_dc = real(result.Y.values[core_index, core_index, 1])
        core_radius = (nominal.stranded_core ? 13nominal.strand_diameter / 2 : nominal.core_diameter / 2)
        r_in = core_radius + nominal.inner_semicon_thickness
        G_analytic = 2pi / (nominal.insulation_rho * log1p(nominal.insulation_thickness / r_in))
        modal = compute(ModalAnalysisProblem(result), ModalAnalysisFormulation(:default))
        segment = PropagationParameters(modal; line_length = nominal.line_length)
        travel_times = nominal.line_length ./ velocity(segment)[:, end]
        time_constants = [maximum(eigvals(Symmetric(imag.(result.Z.values[:, :, k]) / (2pi * result.f[k])),
            Symmetric(real.(result.Z.values[:, :, k])))) for k in 1:3]
        tau = maximum(time_constants)
        step = minimum(travel_times) / 10
        tmax = ceil(max(20tau, 50maximum(travel_times)) / step) * step
        window = (; tmax, step)
        short_window = (tmax = 20minimum(travel_times), step = minimum(travel_times) / 200)
        bonds = Tuple(TerminalImpedance(label; R = 10.0) for label in labels
            if endswith(label, ":sheath") || endswith(label, ":armor"))
        sending = (Step(core1, nominal.voltage), bonds...)
        R_fault = 0.01
        receiving = (TerminalImpedance(core1; R = R_fault),)
        long_problem = TransientAnalysisProblem(problem; sending, receiving, window)
        short_problem = TransientAnalysisProblem(problem; sending, receiving = (), window = short_window)
        summary = Dict{String, Any}("cable" => problem.system.system_id,
            "G_dc_S_per_m" => G_dc, "G_analytic_S_per_m" => G_analytic,
            "G_frequency_Hz" => first(result.f), "tau_seconds" => tau,
            "time_constants_seconds" => time_constants, "time_constant_frequencies_Hz" => result.f[1:3],
            "travel_times_seconds" => travel_times, "duration_seconds" => tmax, "step_seconds" => step,
            "window_choice" => "max(20 times slowest L/R mode at first three baseline frequencies, 50 times longest travel time), rounded up to an output step",
            "short_duration_seconds" => short_window.tmax, "short_step_seconds" => short_window.step,
            "supplementary" => run_dc_supplementary, "selected_runs" => run_dc_selected,
            "fault_resistance_ohm" => R_fault, "fault_current_definition" => "receiving voltage divided by fault resistance; into fault")
        println("DC conductance: numeric ", G_dc, " S/m; main-insulation analytic ", G_analytic, " S/m")
        println("DC window: ", window, "; L/R modes ", time_constants, " s; travel times ", travel_times, " s")
        flush(stdout)

        # Export selected variants locally before optional station execution.
        variants = NamedTuple[]
        for (line_name, line) in ((:lossless, Formulation()), (:lossy, formulation)), dc in (0, 1), passivity_mode in (0, 1)
            push!(variants, (; scenario = :short, line_name, line,
                fit = (enablf = 1, shntcab = 0.0, DCenab = dc, DCCOR = 1,
                    CPASS = passivity_mode, FS = 0.5, YMaxP = 20, AMaxP = 20)))
        end
        for conductance in (G_dc, 10G_dc, 100G_dc, 1000G_dc, 1e-11), dc_correction in (1, 2), passivity_mode in (0, 1)
            push!(variants, (; scenario = :short, line_name = :lossy, line = formulation,
                fit = (enablf = 0, shntcab = conductance, DCenab = 1, DCCOR = dc_correction,
                    CPASS = passivity_mode, FS = 0.5, YMaxP = 20, AMaxP = 20)))
        end
        for dc_correction in (1, 2)
            push!(variants, (; scenario = :short, line_name = :lossy, line = formulation,
                fit = (enablf = 0, shntcab = G_dc, DCenab = 1, DCCOR = dc_correction,
                    CPASS = 1, FS = 0.01, YMaxP = 20, AMaxP = 20)))
        end
        push!(variants, (; scenario = :long, line_name = :lossy, line = formulation,
            fit = (enablf = 1, shntcab = 0.0, DCenab = 0, DCCOR = 1,
                CPASS = 0, FS = 0.5, YMaxP = 20, AMaxP = 20)))
        for dc_correction in (1, 2), passivity_mode in (0, 1)
            push!(variants, (; scenario = :long, line_name = :lossy, line = formulation,
                fit = (enablf = 0, shntcab = G_dc, DCenab = 1, DCCOR = dc_correction,
                    CPASS = passivity_mode, FS = 0.5, YMaxP = 20, AMaxP = 20)))
        end
        # Keep the original 35 indices to preserve run identity.
        for conductance in (G_dc, 100G_dc), dc_correction in (1, 2)
            push!(variants, (; scenario = :short, line_name = :lossy, line = formulation,
                fit = (enablf = 0, shntcab = conductance, DCenab = 1, DCCOR = dc_correction,
                    CPASS = 2, FS = 0.5, YMaxP = 20, AMaxP = 20)))
        end
        for start_frequency in (0.01, 0.001), dc_correction in (1, 2), conductance in (G_dc, 100G_dc), passivity_mode in (1, 2)
            push!(variants, (; scenario = :long, line_name = :lossy, line = formulation,
                fit = (enablf = 0, shntcab = conductance, DCenab = 1, DCCOR = dc_correction,
                    CPASS = passivity_mode, FS = start_frequency, YMaxP = 20, AMaxP = 20)))
        end
        if dc_base_folder !== nothing
            isdir(dc_base_folder) || error("Base run folder is missing: $dc_base_folder")
            base_summary = TOML.parsefile(joinpath(dc_base_folder, "dc_behaviour.toml"))
            base_summary["cable"] == problem.system.system_id || error("Base run folder belongs to a different cable")
            summary["base_folder"] = dc_base_folder
        end
        missing_base_indices = Int[]
        if run_dc_supplementary && dc_base_folder !== nothing
            for i in 1:35
                status_path = joinpath(dc_base_folder, "dc_pscad_$i", "status.toml")
                if !isfile(status_path) || get(TOML.parsefile(status_path), "status", "") ∉
                        ("completed", "nonfinite_response", "line_constants_error", "known_station_failure")
                    push!(missing_base_indices, i)
                end
            end
        end
        variant_indices = run_dc_selected ? dc_run_indices : run_dc_supplementary ?
            vcat(missing_base_indices, collect(36:length(variants))) : eachindex(variants)
        summary["missing_base_variant_indices"] = missing_base_indices
        summary["station_variant_indices"] = collect(variant_indices)
        summary["station_runs"] = length(variant_indices)
        summary["station_simulated_seconds"] = sum(variants[i].scenario === :short ? short_window.tmax : tmax for i in variant_indices)
        open(joinpath(folder, "dc_behaviour.toml"), "w") do io
            TOML.print(io, summary; sorted = true)
        end
        dielectric_groups = [
            [(nominal.inner_semicon_thickness, nominal.inner_semicon_rho, nominal.inner_semicon_eps_r),
                (nominal.insulation_thickness, nominal.insulation_rho, nominal.insulation_eps_r),
                (nominal.outer_semicon_thickness, nominal.outer_semicon_rho, nominal.outer_semicon_eps_r),
                (nominal.water_blocking_thickness, nominal.outer_semicon_rho, nominal.outer_semicon_eps_r)],
            [(nominal.inner_sheath_thickness, nominal.pe_rho, nominal.pe_eps_r),
                (nominal.jacket_thickness, nominal.pp_rho, nominal.pp_eps_r)]]
        expected_losses = zeros(4)
        radius = core_radius
        for (slot, layers) in enumerate(dielectric_groups)
            impedance = 0.0im
            for (thickness, rho, eps_r) in layers
                impedance += log1p(thickness / radius) / (inv(rho) + 2pi * 50im * 8.8541878128e-12 * eps_r)
                radius += thickness
            end
            admittance = 2pi / impedance
            expected_losses[slot] = min(10.0, real(admittance) / imag(admittance))
            if slot == 1
                radius += nominal.sheath_thickness
            elseif slot == 2 && length(dielectric_groups) == 3
                radius += nominal.armour_wire_diameter
            end
        end
        export_start = time_ns()
        for index in variant_indices
            variant = variants[index]
            directory = mkpath(joinpath(folder, "dc_pscad_$(index)"))
            selected = PSCADTransientFormulation(model = :fd_phase; line = variant.line, fit = variant.fit)
            selected_problem = variant.scenario === :short ? short_problem : long_problem
            path = export_data(:pscad, selected_problem, selected; file_name = joinpath(directory, "transient.pscx"))
            document = readxml(path)
            phase = only(findall("//User[@defn='master:Line_FrePhase_Options']", document))
            exported = Dict(p["name"] => parse(Float64, p["value"]) for p in findall("./paramlist/param", phase))
            for (key, value) in pairs(variant.fit)
                exported[string(key)] == round(value; sigdigits = 6) || error("Exported $key differs in $path")
            end
            all(key -> !haskey(exported, key), ("ET_PE", "MER_PE", "MIT_PE", "FSP", "FEP", "NFP", "FDIS")) ||
                error("Passivity defaults must come from the component definition in $path")
            losses = Vector{Float64}[]
            for cable in findall("//User[@defn='master:Cable_Coax']", document)
                values = [parse(Float64, only(findall("./paramlist/param[@name='LT$(slot)']", cable))["value"]) for slot in 1:4]
                if variant.line_name === :lossless
                    all(iszero, values) || error("Lossless LT differs in $path")
                else
                    all(isapprox.(values, round.(expected_losses; sigdigits = 6); rtol = 1e-12, atol = 0)) ||
                        error("Lossy LT differs from the radial-layer calculation in $path")
                end
                push!(losses, values)
            end
            length(losses) == 2 && all(==(first(losses)), losses) || error("Cable LT mismatch in $path")
            status = Dict{String, Any}("status" => "exported", "message" => "Station execution not requested",
                "scenario" => string(variant.scenario), "model" => "fd_phase",
                "fault_resistance_ohm" => R_fault, "fault_current_definition" => "receiving voltage divided by fault resistance; into fault",
                "measured_receiving_current_orientation" => "out of line into fault", "measured_current_to_fault_sign" => 1,
                "fault_current_consistency_assessed" => false,
                "line_formulation" => string(variant.line_name),
                "exported_line_formulation" => Dict(string(k) => string(formula_id(v)) for (k, v) in pairs(variant.line.formulas)),
                "fit" => Dict(string(k) => v for (k, v) in pairs(variant.fit)),
                "exported_fit" => exported,
                "dc_correction_method" => (variant.fit.DCCOR == 1 ? "functional_form" : "add_pole_residue"),
                "passivity_mode" => ("disabled", "scan_only", "enforcement")[variant.fit.CPASS + 1],
                "effective_fit" => merge(Dict{String, Any}(exported), Dict(
                    "ET_PE" => 1e-10, "MER_PE" => 2, "MIT_PE" => 5,
                    "FSP" => 0.001, "FEP" => 1e6, "NFP" => 1000, "FDIS" => 3)),
                "component_defaults_source" => "m7-master.pslx; Definition Line_FrePhase_Options; category Passivity Enforcement; parameter defaults",
                "LT" => first(losses), "project_path" => path,
                "export_check" => true, "line_constants_log_paths" => String[],
                "duration_seconds" => selected_problem.window.tmax, "step_seconds" => selected_problem.window.step)
            if variant.scenario === :long
                status["settled"] = false
                status["settling_assessed"] = false
            end
            open(joinpath(directory, "status.toml"), "w") do io
                TOML.print(io, status; sorted = true)
            end
        end
        summary["export_wall_seconds"] = (time_ns() - export_start) / 1e9
        println("DC export check: ", length(variant_indices), " projects, wall ", summary["export_wall_seconds"], " s")
        flush(stdout)

        reference_start = time_ns()
        N = 8192
        reference = try
            compute(long_problem, TransientFormulation(:nlt;
                options = (line = formulation, N, Tfactor = 2.0)); options = (timing = true,))
        catch exception
            exception isa InterruptException && rethrow()
            directory = mkpath(joinpath(folder, "dc_nlt"))
            open(joinpath(directory, "status.toml"), "w") do io
                TOML.print(io, Dict("status" => "julia_error", "message" => sprint(showerror, exception),
                    "settled" => false, "settling_assessed" => false, "N" => N, "wall_seconds" => (time_ns() - reference_start) / 1e9); sorted = true)
            end
            rethrow()
        end
        reference_wall = (time_ns() - reference_start) / 1e9
        tail = findall(>=(0.9tmax), reference.t)
        fault_current = reference.voltage[2, core_index, tail] ./ R_fault
        reference_final = [mean(fault_current), mean(reference.voltage[1, core_index, tail]),
            mean(reference.voltage[2, core_index, tail])]
        reference_variations = map((fault_current, reference.voltage[1, core_index, tail],
                reference.voltage[2, core_index, tail])) do values
            all(isfinite, values) || return Inf
            change = maximum(values) - minimum(values)
            scale = abs(mean(values))
            return iszero(scale) ? (iszero(change) ? 0.0 : Inf) : change / scale
        end
        variation = first(reference_variations)
        reference_status = Dict("status" => variation < 0.001 ? "completed" : "unsettled_reference",
            "message" => variation < 0.001 ? "" : "NLT final-window variation exceeds 0.1 percent",
            "fault_current_A" => reference_final[1], "sending_voltage_V" => reference_final[2],
            "receiving_voltage_V" => reference_final[3], "tail_variation" => variation,
            "fault_current_tail_variation" => reference_variations[1],
            "sending_voltage_tail_variation" => reference_variations[2],
            "receiving_voltage_tail_variation" => reference_variations[3],
            "settled" => variation < 0.001, "settling_assessed" => true,
            "fault_current_error_A" => 0.0, "fault_current_relative_error" => 0.0,
            "sending_voltage_error_V" => 0.0, "receiving_voltage_error_V" => 0.0,
            "sending_voltage_relative_error" => 0.0, "receiving_voltage_relative_error" => 0.0,
            "N" => N, "T_seconds" => details(reference).data.nlt.T, "wall_seconds" => reference_wall)
        reference_status["fault_resistance_ohm"] = R_fault
        reference_status["fault_current_definition"] = "receiving voltage divided by fault resistance; into fault"
        reference_status["measured_receiving_current_A"] = mean(reference.current[2, core_index, tail])
        reference_status["measured_receiving_current_orientation"] = "into line at receiving end"
        reference_status["measured_current_to_fault_sign"] = -1
        consistency = -1 .* reference.current[2, core_index, tail] .- reference.voltage[2, core_index, tail] ./ R_fault
        consistency_scale = maximum(abs, reference.voltage[2, core_index, tail] ./ R_fault)
        consistency_max = all(isfinite, consistency) ? maximum(abs, consistency) : Inf
        reference_status["fault_current_consistency_error_A"] = mean(consistency)
        reference_status["fault_current_consistency_max_abs_A"] = consistency_max
        reference_status["fault_current_consistency_relative_error"] = !isfinite(consistency_scale) ? Inf :
            iszero(consistency_scale) ? (iszero(consistency_max) ? 0.0 : Inf) : consistency_max / consistency_scale
        reference_status["fault_current_consistency_assessed"] = true
        directory = mkpath(joinpath(folder, "dc_nlt"))
        save(ObservedResult(reference; clip = false), joinpath(directory, "response.json"))
        open(joinpath(directory, "status.toml"), "w") do io
            TOML.print(io, reference_status; sorted = true)
        end
        println("DC NLT: ", reference_status)
        flush(stdout)
        variation < 0.001 || error("Extend the DC window: final NLT current variation is $variation")

        if !(run_dc_supplementary || run_dc_selected)
            # Compare all six fits with the same converged physical DC targets.
            dc_frequencies = [1e-18, 1e-15, 1e-12]
            dc_scan = compute(LineParametersProblem(problem.system, problem.temperature, problem.earth_props,
                dc_frequencies, 0.0), formulation)
            dc_differences = [maximum(abs, values[:, :, 1] - values[:, :, 2]) / maximum(abs, values[:, :, 1])
                for values in (dc_scan.Z.values, dc_scan.Y.values)]
            all(<=(1e-6), dc_differences) || error("Physical DC targets did not converge")
            Z0, Y0 = dc_scan.Z.values[:, :, 1], dc_scan.Y.values[:, :, 1]
            yc0 = Z0 \ sqrt(Z0 * Y0)
            targets = (admittance = real.((yc0 + transpose(yc0)) / 2),
                propagation = real.(exp(-sqrt(Y0 * Z0) * nominal.line_length)))
            # The NLT assembly is private: `src/engine/transientanalysis/formulas/nlt.jl`,
            # functions `sample_network!` and `stamp!`. For current waves (`docs/src/line_models.md`):
            # `Ys = (I - H^2) \ ((I + H^2) * Yc)`, `Ym = (I - H^2) \ (-2H * Yc)`.
            function dc_solution(yc, h)
                n = length(labels)
                blocks = (I - h * h) \ hcat((I + h * h) * yc, -2h * yc)
                ys, ym = blocks[:, 1:n], blocks[:, n+1:2n]
                line_admittance = [ys ym; ym ys]
                matrix = zeros(ComplexF64, 2n + 1, 2n + 1)
                matrix[1:2n, 1:2n] .= line_admittance
                rhs = zeros(ComplexF64, 2n + 1)
                for (offset, assignments) in ((0, sending), (n, receiving)), assignment in assignments
                    node = offset + only(findall(==(assignment.terminal), labels))
                    if assignment isa Step
                        matrix[node, end] = matrix[end, node] = 1
                        rhs[end] = assignment.amplitude
                    else
                        matrix[node, node] += inv(assignment.R)
                    end
                end
                volts = (matrix \ rhs)[1:2n]
                amps = line_admittance * volts
                # Circuit currents enter the line. The fault current follows Ohm's law.
                return real.([volts[n + core_index] / R_fault, volts[core_index], volts[n + core_index]]),
                    real(amps[n + core_index])
            end
            target_dc_final, target_dc_measured = dc_solution(targets.admittance, targets.propagation)
            reference_status["dc_solution_measured_receiving_current_A"] = target_dc_measured
            reference_status["dc_solution_current_consistency_error_A"] = -target_dc_measured - target_dc_final[1]
            for (i, field, unit) in ((1, "fault_current", "A"), (2, "sending_voltage", "V"), (3, "receiving_voltage", "V"))
                reference_status["dc_solution_$(field)_$unit"] = target_dc_final[i]
                reference_status["dc_solution_$(field)_error_$unit"] = target_dc_final[i] - reference_final[i]
                reference_status["dc_solution_$(field)_relative_error"] = (target_dc_final[i] - reference_final[i]) / reference_final[i]
            end
            open(joinpath(folder, "dc_nlt", "status.toml"), "w") do io
                TOML.print(io, reference_status; sorted = true)
            end
            println("DC target check: ", target_dc_final - reference_final)
            flush(stdout)
            fitting_problem = LineModelProblem(segment; line = problem)
            for (index, (method, enforce)) in enumerate(Iterators.product((:none, :add_pole, :functional_form), (false, true)))
                directory = mkpath(joinpath(folder, "dc_ulm_$(index)"))
                started = time_ns()
                status = Dict{String, Any}("method" => string(method), "passivity" => enforce,
                    "status" => "started", "message" => "", "settled" => false, "settling_assessed" => false,
                    "fault_resistance_ohm" => R_fault, "fault_current_definition" => "receiving voltage divided by fault resistance; into fault",
                    "measured_receiving_current_orientation" => "into line at receiving end", "measured_current_to_fault_sign" => -1,
                    "fault_current_consistency_assessed" => false)
                try
                    selected = LineModelFormulation(:ulm; options = (yc = (npoles = 16,),
                        h = (npoles = 16, td = true), passivity = enforce,
                        dc = method === :none ? (;) : (; method), line = formulation))
                    fit_wall = @elapsed fitted = compute(fitting_problem, selected; options = (timing = true,))
                    fitted.coordinates == labels || error("Fitted terminal order differs from the DC circuit")
                    pole_time_constants = [inv(minimum(p -> abs(real(p)), group; init = Inf))
                        for group in LineCableModels.poles(fitted)]
                    status["yc_time_constant_seconds"] = first(pole_time_constants)
                    status["h_time_constants_seconds"] = pole_time_constants[2:end]
                    status["slowest_time_constant_seconds"] = maximum(pole_time_constants)
                    dc_correction = if method === :add_pole
                        diagnostic = details(fitted).data.dc
                        for (quantity, correction) in (("yc", diagnostic.admittance), ("h", diagnostic.propagation))
                            status["$(quantity)_added_pole_a0_per_second"] = correction.a0
                            status["$(quantity)_added_pole_time_constant_seconds"] = inv(correction.a0)
                        end
                        (admittance = diagnostic.admittance, propagation = diagnostic.propagation)
                    else
                        nothing
                    end
                    achieved = (admittance = Yc(fitted, ComplexF64[0])[:, :, 1],
                        propagation = H(fitted, ComplexF64[0])[:, :, 1])
                    dc_final, dc_measured = dc_solution(achieved.admittance, achieved.propagation)
                    # Corrected fits keep the physical targets selected by the package.
                    circuit_targets = something(fitted.dc.targets, targets)
                    target_dc_final, target_dc_measured = dc_solution(circuit_targets.admittance, circuit_targets.propagation)
                    status["dc_solution_measured_receiving_current_A"] = dc_measured
                    status["dc_solution_current_consistency_error_A"] = -dc_measured - dc_final[1]
                    status["target_dc_measured_receiving_current_A"] = target_dc_measured
                    status["target_dc_current_consistency_error_A"] = -target_dc_measured - target_dc_final[1]
                    for (i, field, unit) in ((1, "fault_current", "A"), (2, "sending_voltage", "V"), (3, "receiving_voltage", "V"))
                        status["dc_solution_$(field)_$unit"] = dc_final[i]
                        status["dc_solution_$(field)_error_$unit"] = dc_final[i] - reference_final[i]
                        status["dc_solution_$(field)_relative_error"] = (dc_final[i] - reference_final[i]) / reference_final[i]
                        status["target_dc_$(field)_$unit"] = target_dc_final[i]
                        status["target_dc_$(field)_error_$unit"] = target_dc_final[i] - reference_final[i]
                        status["target_dc_$(field)_relative_error"] = (target_dc_final[i] - reference_final[i]) / reference_final[i]
                    end
                    dc_errors = map((actual, target) -> maximum(abs, actual - target) / maximum(abs, target), achieved, targets)
                    status["fit_wall_seconds"] = fit_wall
                    status["H0_relative_deviation"] = dc_errors.propagation
                    status["Yc0_relative_deviation"] = dc_errors.admittance
                    point = ObservedResult(fitted; clip = false)
                    evidence = (; targets, achieved, relative_deviation = dc_errors, dc_frequencies,
                        dc_differences, target_frequency = first(dc_frequencies), method, enforce,
                        pole_time_constants, dc_correction, dc_final, circuit_targets, target_dc_final, reference_final,
                        model = (admittance = fitted.model.admittance, propagation = fitted.model.propagation),
                        fitting_error = fitting_error(fitted), passivity = passivity(fitted))
                    save(ObservedResult(merge(point.gridpoint, (; dc_behaviour = evidence)),
                        NamedTuple[], point.errors, point.timings), joinpath(directory, "fit.json"))
                    response_wall = @elapsed response = compute(TransientAnalysisProblem(fitted;
                        sending, receiving, window), TransientFormulation(:td); options = (timing = true,))
                    response_tail = findall(>=(0.9tmax), response.t)
                    tail_values = (response.voltage[2, core_index, response_tail] ./ R_fault,
                        response.voltage[1, core_index, response_tail], response.voltage[2, core_index, response_tail])
                    final = mean.(tail_values)
                    status["fault_resistance_ohm"] = R_fault
                    status["fault_current_definition"] = "receiving voltage divided by fault resistance; into fault"
                    status["measured_receiving_current_A"] = mean(response.current[2, core_index, response_tail])
                    status["measured_receiving_current_orientation"] = "into line at receiving end"
                    status["measured_current_to_fault_sign"] = -1
                    consistency = -1 .* response.current[2, core_index, response_tail] .- response.voltage[2, core_index, response_tail] ./ R_fault
                    consistency_scale = maximum(abs, response.voltage[2, core_index, response_tail] ./ R_fault)
                    consistency_max = all(isfinite, consistency) ? maximum(abs, consistency) : Inf
                    status["fault_current_consistency_error_A"] = mean(consistency)
                    status["fault_current_consistency_max_abs_A"] = consistency_max
                    status["fault_current_consistency_relative_error"] = !isfinite(consistency_scale) ? Inf :
                        iszero(consistency_scale) ? (iszero(consistency_max) ? 0.0 : Inf) : consistency_max / consistency_scale
                    status["fault_current_consistency_assessed"] = true
                    for (i, field, unit) in ((1, "fault_current", "A"), (2, "sending_voltage", "V"), (3, "receiving_voltage", "V"))
                        status["$(field)_$unit"] = final[i]
                        status["$(field)_error_$unit"] = final[i] - reference_final[i]
                        values = tail_values[i]
                        change = maximum(values) - minimum(values)
                        scale = abs(final[i])
                        status["$(field)_tail_variation"] = !all(isfinite, values) ? Inf :
                            iszero(scale) ? (iszero(change) ? 0.0 : Inf) : change / scale
                        status["$(field)_relative_error"] = (final[i] - reference_final[i]) / reference_final[i]
                    end
                    status["tail_variation"] = status["fault_current_tail_variation"]
                    status["settled"] = status["tail_variation"] < 0.001
                    status["settling_assessed"] = true
                    status["response_wall_seconds"] = response_wall
                    status["status"] = all(isfinite, final) ? "completed" : "nonfinite_response"
                    status["message"] = all(isfinite, final) ? "" : "Nonfinite final response"
                    save(ObservedResult(response; clip = false), joinpath(directory, "response.json"))
                catch exception
                    exception isa InterruptException && rethrow()
                    status["status"] = "julia_error"
                    status["message"] = sprint(showerror, exception)
                end
                status["wall_seconds"] = (time_ns() - started) / 1e9
                open(joinpath(directory, "status.toml"), "w") do io
                    TOML.print(io, status; sorted = true)
                end
                println("DC ULM: ", status)
                flush(stdout)
            end
        end

        if run_pscad
            remote = RemoteConfig(ENV["LINECABLEMODELS_PSCAD_CONFIG"])
            station = identify(remote)
            open(joinpath(folder, "pscad_station.toml"), "w") do io
                TOML.print(io, station; sorted = true)
            end
            native_folder = mkpath(joinpath(remote.local_root, "$(problem.system.system_id)_$(basename(folder))"))
            if any(i -> variants[i].scenario === :short, variant_indices)
                short_reference = compute(short_problem, TransientFormulation(:nlt;
                    options = (line = formulation, N = 4096, Tfactor = 1.25)))
                save(ObservedResult(short_reference; clip = false), joinpath(folder, "dc_short_reference.json"))
            end
            for index in variant_indices
                variant = variants[index]
                directory = joinpath(folder, "dc_pscad_$(index)")
                status = TOML.parsefile(joinpath(directory, "status.toml"))
                before = Set(readdir(native_folder))
                started = time_ns()
                selected = PSCADTransientFormulation(model = :fd_phase; line = variant.line, fit = variant.fit)
                selected_problem = variant.scenario === :short ? short_problem : long_problem
                response = try
                    compute(selected_problem, selected; options = (; remote, work_root = native_folder))
                catch exception
                    exception isa InterruptException && rethrow()
                    if exception isa NonfiniteOutputError
                        status["status"] = "nonfinite_response"
                        status["message"] = sprint(showerror, exception)
                        status["nonfinite_output_path"] = exception.output_path
                        status["nonfinite_channel"] = exception.channel
                        status["first_nonfinite_time_seconds"] = exception.time
                        status["settled"] = false
                        status["settling_assessed"] = false
                    elseif exception isa LineConstantsError
                        status["status"] = "line_constants_error"
                        status["message"] = exception.message
                        status["log_path"] = exception.log_path
                    elseif exception isa ErrorException
                        status["status"] = "known_station_failure"
                        status["message"] = sprint(showerror, exception)
                    else
                        rethrow()
                    end
                    nothing
                end
                for name in setdiff(Set(readdir(native_folder)), before)
                    cp(joinpath(native_folder, name), joinpath(directory, name))
                end
                status["line_constants_log_paths"] = [joinpath(root, name) for (root, _, files) in walkdir(directory)
                    for name in files if endswith(lowercase(name), ".log")]
                if response !== nothing
                    exact = variant.scenario === :short ? short_reference : reference
                    deviation = compare(exact, response, voltage)
                    status["scaled_rms"] = deviation.absolute[2, core_index] / maximum(abs, exact.voltage[2, core_index, :])
                    status["peak_receiving_voltage_kV"] = maximum(abs, response.voltage[2, core_index, :]) / 1000
                    if variant.scenario === :long
                        response_tail = findall(>=(0.9tmax), response.t)
                        tail_values = (response.voltage[2, core_index, response_tail] ./ R_fault,
                            response.voltage[1, core_index, response_tail], response.voltage[2, core_index, response_tail])
                        final = mean.(tail_values)
                        status["fault_resistance_ohm"] = R_fault
                        status["fault_current_definition"] = "receiving voltage divided by fault resistance; into fault"
                        status["measured_receiving_current_A"] = mean(response.current[2, core_index, response_tail])
                        status["measured_receiving_current_orientation"] = "out of line into fault"
                        status["measured_current_to_fault_sign"] = 1
                        consistency = 1 .* response.current[2, core_index, response_tail] .- response.voltage[2, core_index, response_tail] ./ R_fault
                        consistency_scale = maximum(abs, response.voltage[2, core_index, response_tail] ./ R_fault)
                        consistency_max = all(isfinite, consistency) ? maximum(abs, consistency) : Inf
                        status["fault_current_consistency_error_A"] = mean(consistency)
                        status["fault_current_consistency_max_abs_A"] = consistency_max
                        status["fault_current_consistency_relative_error"] = !isfinite(consistency_scale) ? Inf :
                            iszero(consistency_scale) ? (iszero(consistency_max) ? 0.0 : Inf) : consistency_max / consistency_scale
                        status["fault_current_consistency_assessed"] = true
                        for (i, field, unit) in ((1, "fault_current", "A"), (2, "sending_voltage", "V"), (3, "receiving_voltage", "V"))
                            status["$(field)_$unit"] = final[i]
                            status["$(field)_error_$unit"] = final[i] - reference_final[i]
                            values = tail_values[i]
                            change = maximum(values) - minimum(values)
                            scale = abs(final[i])
                            status["$(field)_tail_variation"] = !all(isfinite, values) ? Inf :
                                iszero(scale) ? (iszero(change) ? 0.0 : Inf) : change / scale
                            status["$(field)_relative_error"] = (final[i] - reference_final[i]) / reference_final[i]
                        end
                        status["tail_variation"] = status["fault_current_tail_variation"]
                        status["settled"] = status["tail_variation"] < 0.001
                        status["settling_assessed"] = true
                    end
                    status["status"] = "completed"
                    status["message"] = ""
                    save(ObservedResult(response; clip = false), joinpath(directory, "response.json"))
                end
                status["station"] = station
                status["wall_seconds"] = (time_ns() - started) / 1e9
                open(joinpath(directory, "status.toml"), "w") do io
                    TOML.print(io, status; sorted = true)
                end
                println("DC PSCAD ", index, ": ", status["status"], ", wall ", status["wall_seconds"], " s")
                flush(stdout)
            end
        end
        summary["reference_wall_seconds"] = reference_wall
        summary["wall_seconds"] = (time_ns() - study_start) / 1e9
        open(joinpath(folder, "dc_behaviour.toml"), "w") do io
            TOML.print(io, summary; sorted = true)
        end
        println("DC behaviour: wall ", summary["wall_seconds"], " s, ", folder)
    end
end
