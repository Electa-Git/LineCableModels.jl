function known_case_problem(::Val{:cable_525kv_land_no_armour_ac_flat};
        frequencies=10.0 .^ range(-1, 7; length=101))
    library = MaterialsLibrary(add_defaults=true)
    core_diameter = 68.0e-3
    core_r20 = 0.006
    inner_semicon_thickness = 1.8e-3
    insulation_thickness = 26.5e-3
    outer_semicon_thickness = 1.5e-3
    water_blocking_thickness = 0.3e-3
    pre_sheath_filler = 3.7e-3
    aluminum_sheath_thickness = 1.2e-3
    post_sheath_filler = 2.2e-3
    jacket_thickness = 5.0e-3
    outer_semicon_skin = 0.3e-3
    xlpe_rho = 1.0e14
    xlpe_eps_r = 2.4
    aluminum_rho = 2.8264e-8
    pe_rho = 1.97e14
    pe_eps_r = 2.5
    cable_x = (-1.0, 0.0, 1.0)
    cable_y = -1.0
    line_length = 20_000.0
    voltage = 525.0
    temperature = 70.0
    earth_rho = 1_000.0
    earth_eps_r = 10.0

    core_rho = core_r20 * (pi * (core_diameter / 2)^2) / 1_000
    core = Material(
        :conductor, core_rho, 1.0, 1.0, 20.0, 0.00393
    )
    semicon_inner = Material(:semicon, 1_000.0, 1_000.0)
    semicon_outer = Material(:semicon, 500.0, 1_000.0)
    xlpe = Material(:insulator, xlpe_rho, xlpe_eps_r)
    aluminum = Material(
        :conductor, aluminum_rho, 1.0, 1.000022, 20.0, 0.00429
    )
    pe = Material(:insulator, pe_rho, pe_eps_r)
    
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
        (:core_water_blocking, water_blocking_thickness, semicon_outer),
        (:core_pre_sheath_filler, pre_sheath_filler, pe)
    )
        push!(parts,
            Region(tag, Shell(layer_thickness), material))
        radius += layer_thickness
    end
    aluminum_outer = radius + aluminum_sheath_thickness
    push!(parts,
        Group(
            :sheath,
            Region(
                :sheath_aluminum,
                Annulus(radius, aluminum_outer),
                aluminum
            )
        ))
    for (tag, layer_thickness, material) in (
        (:sheath_post_filler, post_sheath_filler, pe),
        (:sheath_outer_jacket, jacket_thickness, pe),
        (:sheath_outer_semicon, outer_semicon_skin, semicon_outer)
    )
        push!(parts,
            Region(tag, Shell(layer_thickness), material))
    end
    
    design = build(
        CableDesign,
        "cable_525kv_land_no_armour",
        Stack(parts);
        nominal_data = (
            designation_code = "cable_525kv_land_no_armour",
            U0 = voltage / sqrt(3),
            U = voltage,
            conductor_cross_section = 3_000.0
        )
    )
    earth = homogeneous(
        rho = earth_rho, eps_r = earth_eps_r, mu_r = 1.0
    )
    connections = [Dict(:core => 2index - 1, :sheath => 2index) for index in 1:3]
    system = build(
        LineCableSystem,
        fill(design, 3),
        [Pose2(cable_x[index], cable_y) for index in 1:3];
        connections,
        environment = earth,
        system_id = "cable_525kv_land_no_armour_ac_flat",
        line_length = line_length
    )
    return LineParametersProblem(system; temperature, earth_props=earth, frequencies)
end
