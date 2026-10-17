function known_case_problem(::Val{:cable_525kv_subsea_armoured_ac_flat};
        frequencies=10.0 .^ range(-1, 7; length=101))
    library = MaterialsLibrary(add_defaults=true)
    core_diameter = 0.06
    core_r20 = 0.0072
    inner_semicon_thickness = 0.002
    insulation_thickness = 0.026
    outer_semicon_thickness = 0.0018
    water_blocking_thickness = 0.0006
    lead_sheath_thickness = 0.0032
    inner_pe_thickness = 0.0033
    bedding_thickness = 0.001
    armor_wires = 68
    armor_wire_diameter = 0.006
    armor_lay_ratio = 10.019920548139554
    jacket_thickness = 0.0066
    xlpe_rho = 1.0e14
    xlpe_eps_r = 2.4
    lead_rho = 2.14e-7
    lead_mu_r = 0.999983
    pe_rho = 1.97e14
    pe_eps_r = 2.5
    pp_rho = 1.0e15
    pp_eps_r = 2.8
    steel_rho = 1.38e-7
    steel_mu_r = 10.0
    cable_x = (-1.0, 0.0, 1.0)
    cable_y = -1.0
    line_length = 20000.0
    voltage = 525.0
    temperature = 70.0
    earth_rho = 1000.0
    earth_eps_r = 10.0

    core_rho = core_r20 * (pi * (core_diameter / 2)^2) / 1_000
    core = Material(
        :conductor, core_rho, 1.0, 1.0, 20.0, 0.00393
    )
    semicon_inner = Material(:semicon, 1_000.0, 1_000.0)
    semicon_outer = Material(:semicon, 500.0, 1_000.0)
    xlpe = Material(:insulator, xlpe_rho, xlpe_eps_r)
    lead = Material(
        :conductor, lead_rho, 1.0, lead_mu_r, 20.0, 0.004
    )
    pe = Material(:insulator, pe_rho, pe_eps_r)
    pp = Material(:insulator, pp_rho, pp_eps_r)
    steel = Material(
        :conductor, steel_rho, 1.0, steel_mu_r, 20.0, 0.0045
    )
    matrix = Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )

    radius = core_diameter / 2
    parts = AbstractCablePart[
    Group(
        :core,
        Region(:core_conductor, Disk(radius), core)
    ),
]
    for (tag, layer_thickness, material) in (
        (:core_semicon_inner, inner_semicon_thickness, semicon_inner),
        (:core_insulation, insulation_thickness, xlpe),
        (:core_semicon_outer, outer_semicon_thickness, semicon_outer),
        (:core_water_blocking, water_blocking_thickness, semicon_outer)
    )
        push!(parts,
            Region(tag, Shell(layer_thickness), material))
        radius += layer_thickness
    end
    lead_outer = radius + lead_sheath_thickness
    push!(parts,
        Group(
            :sheath,
            Region(
                :sheath_lead,
                Annulus(radius, lead_outer),
                lead
            )
        ))
    push!(parts,
        Region(
            :sheath_inner_pe,
            Shell(inner_pe_thickness),
            pe
        ))
    push!(parts,
        Region(
            :sheath_bedding,
            Shell(bedding_thickness),
            pp
        ))
    armor_wire_radius = armor_wire_diameter / 2
    armor_inner_radius = lead_outer + inner_pe_thickness + bedding_thickness
    armor_outer = armor_inner_radius + 2armor_wire_radius
    armor = Group(
        :armor,
        Region(
            :armor_wires,
            Disk(armor_wire_radius),
            steel
        );
        pattern = Ring(
            armor_wires; r = armor_inner_radius + armor_wire_radius
        ),
        path = Helix(LayRatio(armor_lay_ratio))
    )
    push!(parts,
        Enclosure(
            :armor_matrix,
            armor;
            primitive = Annulus(armor_inner_radius, armor_outer),
            fill = matrix
        ))
    push!(parts,
        Region(
            :armor_jacket,
            Shell(jacket_thickness),
            pp
        ))

    design = build(
        CableDesign,
        "cable_525kv_subsea_armoured",
        Stack(parts);
        nominal_data = (
            designation_code = "cable_525kv_subsea_armoured",
            U0 = voltage / sqrt(3),
            U = voltage,
            conductor_cross_section = 2_500.0
        )
    )
    earth = homogeneous(
        rho = earth_rho, eps_r = earth_eps_r, mu_r = 1.0
    )
    connections = [Dict(:core => 3index - 2, :sheath => 3index - 1, :armor => 3index)
                   for index in 1:3]
    system = build(
        LineCableSystem,
        fill(design, 3),
        [Pose2(cable_x[index], cable_y) for index in 1:3];
        connections,
        environment = earth,
        system_id = "cable_525kv_subsea_armoured_ac_flat",
        line_length = line_length
    )
    return LineParametersProblem(system; temperature, earth_props=earth, frequencies)
end
