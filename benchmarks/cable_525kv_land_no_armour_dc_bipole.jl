# Solid equivalent core, aluminium sheath, unarmoured.
# Two DC poles, 150 km long.
# Values from the benchmark specification, Appendix A (DC525-L).

using LineCableModels
using Dates

frequencies = sort!(vcat(exp10.(range(-3, 6; length = 226)), 50.0))

"""
Construct the cable line-parameter problem from physical data in SI units and temperatures in °C.
"""
function line_parameters_problem(
        core_diameter, core_dc_resistance, inner_semicon_thickness, insulation_thickness,
        outer_semicon_thickness, water_blocking_thickness, pre_sheath_filler_thickness,
        sheath_thickness, post_sheath_filler_thickness, jacket_thickness,
        outer_semicon_skin_thickness, core_mu_r, core_temperature_coefficient,
        inner_semicon_rho, inner_semicon_eps_r, outer_semicon_rho, outer_semicon_eps_r,
        insulation_rho, insulation_eps_r, insulation_mu_r, sheath_rho, sheath_mu_r,
        sheath_temperature_coefficient, pe_rho, pe_eps_r, reference_temperature, cable_x_left,
        cable_x_right, cable_y, line_length, voltage, conductor_cross_section, temperature,
        earth_rho, earth_eps_r, earth_mu_r, frequencies,
)
    core_rho = core_dc_resistance * (pi * (core_diameter / 2)^2)
    core = Material(:conductor, core_rho, 1.0, core_mu_r,
        reference_temperature, core_temperature_coefficient)
    semicon_inner = Material(:semicon, inner_semicon_rho, inner_semicon_eps_r)
    semicon_outer = Material(:semicon, outer_semicon_rho, outer_semicon_eps_r)
    xlpe = Material(:insulator, insulation_rho, insulation_eps_r, insulation_mu_r)
    sheath = Material(:conductor, sheath_rho, 1.0, sheath_mu_r,
        reference_temperature, sheath_temperature_coefficient)
    pe = Material(:insulator, pe_rho, pe_eps_r)

    radius = core_diameter / 2
    parts = AbstractCablePart[
        Group(:core, Region(:core_conductor, Disk(radius), core)),
    ]
    for (tag, layer_thickness, material) in (
        (:core_semicon_inner, inner_semicon_thickness, semicon_inner),
        (:core_insulation, insulation_thickness, xlpe),
        (:core_semicon_outer, outer_semicon_thickness, semicon_outer),
        (:core_water_blocking, water_blocking_thickness, semicon_outer),
        (:core_pre_sheath_filler, pre_sheath_filler_thickness, pe),
    )
        push!(parts, Region(tag, Shell(layer_thickness), material))
        radius += layer_thickness
    end
    sheath_outer = radius + sheath_thickness
    push!(parts, Group(:sheath, Region(:sheath_aluminum, Annulus(radius, sheath_outer), sheath)))
    for (tag, layer_thickness, material) in (
        (:sheath_post_filler, post_sheath_filler_thickness, pe),
        (:sheath_outer_jacket, jacket_thickness, pe),
        (:sheath_outer_semicon, outer_semicon_skin_thickness, semicon_outer),
    )
        push!(parts, Region(tag, Shell(layer_thickness), material))
    end

    design = build(CableDesign, "cable_525kv_land_no_armour", Stack(parts);
        nominal_data = (designation_code = "cable_525kv_land_no_armour",
            U0 = voltage / 1000,
            conductor_cross_section = conductor_cross_section / 1e-6))
    earth = homogeneous(rho = earth_rho, eps_r = earth_eps_r, mu_r = earth_mu_r)
    cable_x = (cable_x_left, cable_x_right)
    connections = [Dict(:core => 2index - 1, :sheath => 2index)
                   for index in 1:2]
    system = build(LineCableSystem, fill(design, 2),
        [Pose2(cable_x[index], cable_y) for index in 1:2];
        connections, environment = earth, system_id = "cable_525kv_land_no_armour_dc_bipole",
        line_length = line_length)
    return LineParametersProblem(system; temperature, earth_props = earth, frequencies)
end

nominal = (;
    core_diameter = 68e-3, core_dc_resistance = 0.006e-3, inner_semicon_thickness = 1.8e-3,
    insulation_thickness = 26.5e-3, outer_semicon_thickness = 1.5e-3,
    water_blocking_thickness = 0.3e-3, pre_sheath_filler_thickness = 3.7e-3,
    sheath_thickness = 1.2e-3, post_sheath_filler_thickness = 2.2e-3, jacket_thickness = 5e-3,
    outer_semicon_skin_thickness = 0.3e-3, core_mu_r = 1.0,
    core_temperature_coefficient = 0.00393, inner_semicon_rho = 1000.0,
    inner_semicon_eps_r = 1000.0, outer_semicon_rho = 500.0, outer_semicon_eps_r = 1000.0,
    insulation_rho = 1e14, insulation_eps_r = 2.4, insulation_mu_r = 1.0,
    sheath_rho = 2.8264e-8, sheath_mu_r = 1.000022, sheath_temperature_coefficient = 0.00429,
    pe_rho = 1.97e14, pe_eps_r = 2.5, reference_temperature = 20.0, cable_x_left = -0.5,
    cable_x_right = 0.5, cable_y = -1.0, line_length = 150000.0, voltage = 525e3,
    conductor_cross_section = 3000e-6, temperature = 70.0, earth_rho = 1000.0,
    earth_eps_r = 10.0, earth_mu_r = 1.0, frequencies = frequencies,
)

problem = line_parameters_problem(nominal...)
formulation = Formulation(earth_impedance = :unified, earth_admittance = :unified,
    shunt_model = :equivalent, insulation_admittance = :lossy, semicon_admittance = :lossy;
    options = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
compute_time = @elapsed result = compute(problem, formulation)
folder = joinpath(pwd(), ".linecablemodels", "benchmarks",
    "cable_525kv_land_no_armour_dc_bipole", Dates.format(now(), "yyyymmdd-HHMMSS"))
mkpath(folder)
export_data(:json, problem; file_name = joinpath(folder, "problem.json"))
save(ObservedResult(result, (R, X, G, B); length_unit = :base, clip = false),
    joinpath(folder, "baseline.json"))
println("cable_525kv_land_no_armour_dc_bipole: ", length(problem.system.terminal_order),
    " terminals, compute ", compute_time, " s, ", folder)
