# Explicit stranded core, lead sheath, steel-wire armour.
# Two DC poles, 150 km long.

using LineCableModels
using LineCableModels.Engine: compare
using LineCableModels.ReportBuilder: BenchmarkTableDefinition
using LinearAlgebra: diag
using Dates
using Measurements, Distributions
using Statistics: mean, std

ARGS in (String[], ["time_domain"], ["pscad"], ["time_domain", "pscad"], ["model_settings"]) ||
    error("Usage: julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl [time_domain] [pscad] | model_settings")
run_time_domain = "time_domain" in ARGS
run_pscad = "pscad" in ARGS
run_model_settings = ARGS == ["model_settings"]

frequencies = sort!(vcat(exp10.(range(-3, 6; length = 226)), 50.0))

"""
Construct the cable line-parameter problem from physical data in SI units and temperatures in °C.
"""
function line_parameters_problem(
        strand_diameter, core_lay_ratio, inner_semicon_thickness,
        insulation_thickness, outer_semicon_thickness, water_blocking_thickness,
        sheath_thickness, inner_sheath_thickness, bedding_thickness, armour_wires,
        armour_wire_diameter, armour_lay_ratio, jacket_thickness, core_rho, core_mu_r,
        core_temperature_coefficient, inner_semicon_rho, inner_semicon_eps_r,
        outer_semicon_rho, outer_semicon_eps_r, insulation_rho, insulation_eps_r,
        insulation_mu_r, sheath_rho, sheath_mu_r, sheath_temperature_coefficient, pe_rho,
        pe_eps_r, pp_rho, pp_eps_r, pp_mu_r, armour_rho, armour_mu_r,
        armour_temperature_coefficient, reference_temperature, cable_x_left, cable_x_right,
        cable_y, line_length, voltage, conductor_cross_section, temperature, earth_rho,
        earth_eps_r, earth_mu_r, frequencies,
        stranded_core, core_diameter, armour_as_tube,
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
    steel = Material(:conductor, armour_rho, 1.0, armour_mu_r,
        reference_temperature, armour_temperature_coefficient)
    matrix = Material(kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0)

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
    armour_wire_radius = armour_wire_diameter / 2
    radius_under_bedding = sheath_outer + inner_sheath_thickness
    effective_bedding_thickness = armour_as_tube ?
        bedding_thickness : max(bedding_thickness,
            armour_wire_radius / sin(pi / armour_wires) - armour_wire_radius - radius_under_bedding)
    if !iszero(effective_bedding_thickness)
        push!(parts, Region(:sheath_bedding, Shell(effective_bedding_thickness), pp))
    end
    armour_inner_radius = radius_under_bedding + effective_bedding_thickness
    if armour_as_tube
        armour_outer = sqrt(armour_inner_radius^2 + armour_wires * armour_wire_radius^2)
        push!(parts, Group(:armor,
            Region(:armor_tube, Annulus(armour_inner_radius, armour_outer), steel)))
    else
        armour_outer = armour_inner_radius + 2armour_wire_radius
        armor = Group(:armor, Region(:armor_wires, Disk(armour_wire_radius), steel);
            pattern = Ring(armour_wires; r = armour_inner_radius + armour_wire_radius),
            path = Helix(LayRatio(armour_lay_ratio)))
        push!(parts, Enclosure(:armor_matrix, armor;
            primitive = Annulus(armour_inner_radius, armour_outer), fill = matrix))
    end
    if !iszero(jacket_thickness)
        push!(parts, Region(:armor_jacket, Shell(jacket_thickness), pp))
    end

    design = build(CableDesign, "cable_320kv_armoured", Stack(parts);
        nominal_data = (designation_code = "cable_320kv_armoured",
            U0 = voltage / 1000,
            conductor_cross_section = conductor_cross_section / 1e-6))
    earth = homogeneous(rho = earth_rho, eps_r = earth_eps_r, mu_r = earth_mu_r)
    cable_x = (cable_x_left, cable_x_right)
    connections = [Dict(terminal => length(design.terminal_order) * (index - 1) + offset
                        for (offset, terminal) in enumerate(design.terminal_order))
                   for index in 1:2]
    system = build(LineCableSystem, fill(design, 2),
        [Pose2(cable_x[index], cable_y) for index in 1:2];
        connections, environment = earth, system_id = "cable_320kv_armoured_dc_bipole",
        line_length = line_length)
    return LineParametersProblem(system; temperature, earth_props = earth, frequencies)
end

nominal = (;
    strand_diameter = 3.66e-3, core_lay_ratio = 15.0,
    inner_semicon_thickness = 2e-3, insulation_thickness = 20e-3,
    outer_semicon_thickness = 1.5e-3, water_blocking_thickness = 0.3e-3,
    sheath_thickness = 3.3e-3, inner_sheath_thickness = 3e-3, bedding_thickness = 6.210e-3,
    armour_wires = 68, armour_wire_diameter = 5.82e-3, armour_lay_ratio = 15.0,
    jacket_thickness = 10e-3, core_rho = 1.724e-8, core_mu_r = 1.0,
    core_temperature_coefficient = 0.00393, inner_semicon_rho = 1000.0,
    inner_semicon_eps_r = 1000.0, outer_semicon_rho = 500.0, outer_semicon_eps_r = 1000.0,
    insulation_rho = 1e14, insulation_eps_r = 2.3, insulation_mu_r = 1.0, sheath_rho = 2.14e-7,
    sheath_mu_r = 1.0, sheath_temperature_coefficient = 0.004, pe_rho = 1.97e14,
    pe_eps_r = 2.3, pp_rho = 1e15, pp_eps_r = 2.8, pp_mu_r = 1.0, armour_rho = 1.38e-7,
    armour_mu_r = 300.0, armour_temperature_coefficient = 0.0045, reference_temperature = 20.0,
    cable_x_left = -0.5, cable_x_right = 0.5, cable_y = -1.0, line_length = 150000.0,
    voltage = 320e3, conductor_cross_section = 1600e-6, temperature = 70.0, earth_rho = 1000.0,
    earth_eps_r = 10.0, earth_mu_r = 1.0, frequencies = frequencies,
    stranded_core = true, core_diameter = 46.3e-3, armour_as_tube = false,
)

problem = line_parameters_problem(nominal...)
formulation = Formulation(earth_impedance = :unified, earth_admittance = :unified,
    shunt_model = :equivalent, insulation_admittance = :lossy, semicon_admittance = :lossy;
    options = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
compute_time = @elapsed result = compute(problem, formulation)
folder = joinpath(pwd(), ".linecablemodels", "benchmarks",
    "cable_320kv_armoured_dc_bipole", Dates.format(now(), "yyyymmdd-HHMMSS"))
mkpath(folder)
export_data(:json, problem; file_name = joinpath(folder, "problem.json"))
save(ObservedResult(result, (R, X, G, B); length_unit = :base, clip = false),
    joinpath(folder, "baseline.json"))
println("cable_320kv_armoured_dc_bipole: ", length(problem.system.terminal_order),
    " terminals, compute ", compute_time, " s, ", folder)

if isempty(ARGS)
# One at a time, temperature, and earth resistivity.
swept = (
    :strand_diameter, :core_lay_ratio, :inner_semicon_thickness, :insulation_thickness,
    :outer_semicon_thickness, :water_blocking_thickness, :sheath_thickness,
    :inner_sheath_thickness, :bedding_thickness, :armour_wire_diameter, :armour_lay_ratio,
    :jacket_thickness, :core_rho, :core_mu_r, :core_temperature_coefficient,
    :inner_semicon_rho, :inner_semicon_eps_r, :outer_semicon_rho, :outer_semicon_eps_r,
    :insulation_rho, :insulation_eps_r, :insulation_mu_r, :sheath_rho, :sheath_mu_r,
    :sheath_temperature_coefficient, :pe_rho, :pe_eps_r, :pp_rho, :pp_eps_r, :pp_mu_r,
    :armour_rho, :armour_mu_r, :armour_temperature_coefficient,
)

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
    names = (:conductor, :semicons, :armour, :dielectric, :shunt)
    ranking = NamedTuple[]
    section_time = @elapsed begin
        run = [
            compute(line_parameters_problem(merge(nominal, (; stranded_core = false, core_rho = 2.3853e-8,
                core_mu_r = 1.01663, insulation_eps_r = 2.77781, insulation_mu_r = 1.59891,
                sheath_mu_r = 0.999983, pe_eps_r = 2.51862, pp_mu_r = 1.1263,
                armour_rho = 1.7475e-7, armour_mu_r = 36.6326))...), formulation),
            compute(line_parameters_problem(merge(nominal, (;
                inner_semicon_thickness = 0.0, outer_semicon_thickness = 0.0,
                water_blocking_thickness = 0.0,
                insulation_thickness = nominal.insulation_thickness + nominal.inner_semicon_thickness +
                    nominal.outer_semicon_thickness + nominal.water_blocking_thickness))...), formulation),
            compute(line_parameters_problem(merge(nominal, (; armour_as_tube = true))...), formulation),
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

if run_time_domain || run_pscad
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
        # Zip the sheath and armour grids so all four bonds have the same R.
        bonds = [TerminalImpedance(label; R = Grid((1.0, 10.0, 20.0)))
            for label in labels if endswith(label, ":sheath") || endswith(label, ":armor")]
        reference_problems = TransientAnalysisProblem[TransientAnalysisProblem(problem;
            sending = (Step(core1, nominal.voltage), assignments...), receiving = (), window)
            for assignments in zip(bonds...)]
        both_end_bonds = Tuple(TerminalImpedance(label; R = 10.0)
            for label in labels if endswith(label, ":sheath") || endswith(label, ":armor"))
        push!(reference_problems, TransientAnalysisProblem(problem;
            sending = (Step(core1, nominal.voltage), both_end_bonds...),
            receiving = both_end_bonds, window))
        alpha, beta = 1.47e4, 2.47e6
        pulse_peak_time = log(beta / alpha) / (beta - alpha)
        pulse_amplitude = nominal.voltage / (exp(-alpha * pulse_peak_time) - exp(-beta * pulse_peak_time))
        pulse_step = 0.25e-6
        pulse_tmax = ceil(max(3e-3, 3minimum(travel_times)) / pulse_step) * pulse_step
        pulse_window = (tmax = pulse_tmax, step = pulse_step)
        pulse_N = nextpow(2, ceil(Int, 2e6 * (1.25pulse_tmax) + 0.5))
        push!(reference_problems, TransientAnalysisProblem(problem;
            sending = (DoubleExponential(core1, pulse_amplitude, alpha, beta), both_end_bonds...),
            receiving = (), window = pulse_window))
        references = Vector{TransientResponse}(undef, length(reference_problems))
        reference_points = ObservedResult[]
        for (reference_index, reference_problem) in enumerate(reference_problems)
            run_time_domain || reference_index in (2, 5) || continue
            resistance = (1.0, 10.0, 20.0, 10.0, 10.0)[reference_index]
            scenario = reference_index <= 3 ? :step_single_point : reference_index == 4 ? :step_both_ends : :impulse
            N = scenario === :impulse ? pulse_N : 4096
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
                    bonds = [bond.terminal for bond in reference_problem.sending[2:end]]),
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
            native_folder = mkpath(joinpath(options.remote.local_root, "cable_320kv_armoured_dc_bipole_$(basename(folder))"))
            for (scenario, reference_index) in ((:step_single_point, 2), (:impulse, 5))
                reference = references[reference_index]
                reference_problem = reference_problems[reference_index]
                core_index = findfirst(==(core1), reference.coordinates)
                pscad_formulations = PSCADTransientFormulation[PSCADTransientFormulation(model = :fd_phase; line = formulation,
                    fit = (; YMaxP = order, AMaxP = order, DCenab = dc, CPASS = enforce))
                    for order in (10, 20, 30) for dc in (0, 1) for enforce in (0, 1)]
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
            for label in labels if endswith(label, ":sheath") || endswith(label, ":armor"))
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
