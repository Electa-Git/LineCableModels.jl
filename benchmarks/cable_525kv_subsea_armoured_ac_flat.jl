# Solid equivalent core, lead sheath, steel-wire armour.
# Three AC phases in a flat arrangement, 20 km long.
# Values from the benchmark specification, Appendix A (AC525-S).

using LineCableModels
using Dates

frequencies = sort!(vcat(exp10.(range(-3, 6; length = 226)), 50.0))

"""
Construct the cable line-parameter problem from physical data in SI units and temperatures in °C.
"""
function line_parameters_problem(
        core_diameter, core_dc_resistance, inner_semicon_thickness, insulation_thickness,
        outer_semicon_thickness, water_blocking_thickness, sheath_thickness,
        inner_sheath_thickness, bedding_thickness, armour_wires, armour_wire_diameter,
        armour_length_factor, jacket_thickness, core_mu_r, core_temperature_coefficient,
        inner_semicon_rho, inner_semicon_eps_r, outer_semicon_rho, outer_semicon_eps_r,
        insulation_rho, insulation_eps_r, insulation_mu_r, sheath_rho, sheath_mu_r,
        sheath_temperature_coefficient, pe_rho, pe_eps_r, pp_rho, pp_eps_r, pp_mu_r,
        armour_rho, armour_mu_r, armour_temperature_coefficient, reference_temperature,
        cable_x_left, cable_x_middle, cable_x_right, cable_y, line_length, voltage,
        conductor_cross_section, temperature, earth_rho, earth_eps_r, earth_mu_r, frequencies,
)
    core_rho = core_dc_resistance * (pi * (core_diameter / 2)^2)
    armour_lay_ratio = pi / sqrt(armour_length_factor^2 - 1)
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

    radius = core_diameter / 2
    parts = AbstractCablePart[
        Group(:core, Region(:core_conductor, Disk(radius), core)),
    ]
    for (tag, layer_thickness, material) in (
        (:core_semicon_inner, inner_semicon_thickness, semicon_inner),
        (:core_insulation, insulation_thickness, xlpe),
        (:core_semicon_outer, outer_semicon_thickness, semicon_outer),
        (:core_water_blocking, water_blocking_thickness, semicon_outer),
    )
        push!(parts, Region(tag, Shell(layer_thickness), material))
        radius += layer_thickness
    end
    sheath_outer = radius + sheath_thickness
    push!(parts, Group(:sheath, Region(:sheath_lead, Annulus(radius, sheath_outer), sheath)))
    push!(parts, Region(:sheath_inner_pe, Shell(inner_sheath_thickness), pe))
    push!(parts, Region(:sheath_bedding, Shell(bedding_thickness), pp))
    armour_wire_radius = armour_wire_diameter / 2
    armour_inner_radius = sheath_outer + inner_sheath_thickness + bedding_thickness
    armour_outer = armour_inner_radius + 2armour_wire_radius
    armor = Group(:armor, Region(:armor_wires, Disk(armour_wire_radius), steel);
        pattern = Ring(armour_wires; r = armour_inner_radius + armour_wire_radius),
        path = Helix(LayRatio(armour_lay_ratio)))
    push!(parts, Enclosure(:armor_matrix, armor;
        primitive = Annulus(armour_inner_radius, armour_outer), fill = matrix))
    push!(parts, Region(:armor_jacket, Shell(jacket_thickness), pp))

    design = build(CableDesign, "cable_525kv_subsea_armoured", Stack(parts);
        nominal_data = (designation_code = "cable_525kv_subsea_armoured",
            U0 = voltage / (1000sqrt(3)), U = voltage / 1000,
            conductor_cross_section = conductor_cross_section / 1e-6))
    earth = homogeneous(rho = earth_rho, eps_r = earth_eps_r, mu_r = earth_mu_r)
    cable_x = (cable_x_left, cable_x_middle, cable_x_right)
    connections = [Dict(:core => 3index - 2, :sheath => 3index - 1, :armor => 3index)
                   for index in 1:3]
    system = build(LineCableSystem, fill(design, 3),
        [Pose2(cable_x[index], cable_y) for index in 1:3];
        connections, environment = earth, system_id = "cable_525kv_subsea_armoured_ac_flat",
        line_length = line_length)
    return LineParametersProblem(system; temperature, earth_props = earth, frequencies)
end

nominal = (;
    core_diameter = 60e-3, core_dc_resistance = 0.0072e-3, inner_semicon_thickness = 2e-3,
    insulation_thickness = 26e-3, outer_semicon_thickness = 1.8e-3,
    water_blocking_thickness = 0.6e-3, sheath_thickness = 3.2e-3,
    inner_sheath_thickness = 3.3e-3, bedding_thickness = 1e-3, armour_wires = 68,
    armour_wire_diameter = 6e-3, armour_length_factor = 1.048, jacket_thickness = 6.6e-3,
    core_mu_r = 1.0, core_temperature_coefficient = 0.00393, inner_semicon_rho = 1000.0,
    inner_semicon_eps_r = 1000.0, outer_semicon_rho = 500.0, outer_semicon_eps_r = 1000.0,
    insulation_rho = 1e14, insulation_eps_r = 2.4, insulation_mu_r = 1.0, sheath_rho = 2.14e-7,
    sheath_mu_r = 0.999983, sheath_temperature_coefficient = 0.004, pe_rho = 1.97e14,
    pe_eps_r = 2.5, pp_rho = 1e15, pp_eps_r = 2.8, pp_mu_r = 1.0, armour_rho = 1.38e-7,
    armour_mu_r = 10.0, armour_temperature_coefficient = 0.0045, reference_temperature = 20.0,
    cable_x_left = -1.0, cable_x_middle = 0.0, cable_x_right = 1.0, cable_y = -1.0,
    line_length = 20000.0, voltage = 525e3, conductor_cross_section = 2500e-6,
    temperature = 70.0, earth_rho = 1000.0, earth_eps_r = 10.0, earth_mu_r = 1.0,
    frequencies = frequencies,
)

problem = line_parameters_problem(nominal...)
formulation = Formulation(earth_impedance = :unified, earth_admittance = :unified,
    shunt_model = :equivalent, insulation_admittance = :lossy, semicon_admittance = :lossy;
    options = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
compute_time = @elapsed result = compute(problem, formulation)
folder = joinpath(pwd(), ".linecablemodels", "benchmarks",
    "cable_525kv_subsea_armoured_ac_flat", Dates.format(now(), "yyyymmdd-HHMMSS"))
mkpath(folder)
export_data(:json, problem; file_name = joinpath(folder, "problem.json"))
save(ObservedResult(result, (R, X, G, B); length_unit = :base, clip = false),
    joinpath(folder, "baseline.json"))
println("cable_525kv_subsea_armoured_ac_flat: ", length(problem.system.terminal_order),
    " terminals, compute ", compute_time, " s, ", folder)
