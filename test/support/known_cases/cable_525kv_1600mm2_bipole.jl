function known_case_problem(::Val{:cable_525kv_1600mm2_bipole};
        frequencies=10.0 .^ range(-1, 7; length=101))
    library = MaterialsLibrary(add_defaults=true)
    strand_layers = 6
    strands_per_layer = 6
    strand_diameter = 0.0036649
    core_lay_ratio = 11.0
    inner_semicon_thickness = 0.002
    insulation_thickness = 0.026
    outer_semicon_thickness = 0.0018
    water_blocking_thickness = 0.0003
    lead_screen_thickness = 0.0033
    inner_sheath_thickness = 0.003
    bedding_thickness = 0.003
    armor_wires = 68
    armor_wire_diameter = 0.005827
    armor_lay_ratio = 10.0
    armor_packing_clearance_ratio = 0.0
    jacket_thickness = 0.01
    cable_x = (-0.5, 0.5)
    cable_y = -1.0
    line_length = 1000.0
    temperature = 20.0
    earth_rho = 100.0
    earth_eps_r = 10.0

    copper = Material(library, :copper)
    semicon1 = Material(library, :semicon1)
    semicon2 = Material(library, :semicon2)
    pe = Material(library, :pe)
    polyacrylate = Material(library, :polyacrylate)
    lead = Material(library, :lead)
    pp = Material(library, :pp)
    steel = Material(library, :steel)
    matrix = Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )
    radius_before_bedding = (strand_layers + 0.5) * strand_diameter +
                            inner_semicon_thickness +
                            insulation_thickness +
                            outer_semicon_thickness +
                            water_blocking_thickness +
                            lead_screen_thickness +
                            inner_sheath_thickness
    armor_wire_radius = armor_wire_diameter / 2
    unbuffered_armor_radius = radius_before_bedding + bedding_thickness
    unbuffered_outer_radius = unbuffered_armor_radius + armor_wire_diameter +
                              jacket_thickness
    packing_buffer = armor_packing_clearance_ratio * unbuffered_outer_radius
    effective_bedding_thickness = bedding_thickness + packing_buffer
    strand_radius = strand_diameter / 2
    radius = (2strand_layers + 1) * strand_radius
    stranded_core = stranded(
        copper;
        shape = Disk(strand_radius),
        lay = LayRatio(core_lay_ratio),
        boundary = Disk(radius),
        fill = matrix
    )
    core = stranded_core
    parts = AbstractCablePart[core]
    for (tag, thickness, material) in (
        (:core_semicon_inner, inner_semicon_thickness, semicon1),
        (:core_insulation, insulation_thickness, pe),
        (:core_semicon_outer, outer_semicon_thickness, semicon2),
        (:core_water_blocking, water_blocking_thickness, polyacrylate)
    )
        push!(parts, Region(
            tag, Shell(thickness), material
        ))
        radius += thickness
    end

    lead_outer = radius + lead_screen_thickness
    push!(parts,
        Group(
            :sheath,
            Region(
                :sheath_lead_screen,
                Annulus(radius, lead_outer),
                lead
            )
        ))
    radius = lead_outer
    push!(parts,
        Region(
            :sheath_inner,
            Shell(inner_sheath_thickness),
            pe
        ))
    radius += inner_sheath_thickness
    push!(parts,
        Region(
            :sheath_bedding,
            Shell(effective_bedding_thickness),
            pp
        ))
    radius += effective_bedding_thickness

    armor_outer = radius + 2armor_wire_radius
    armor_centre = radius + armor_wire_radius
    armor = Group(
        :armor,
        Region(
            :armor_wires, Disk(armor_wire_radius), steel
        );
        pattern = Ring(armor_wires; r = armor_centre),
        path = Helix(LayRatio(armor_lay_ratio))
    )
    push!(parts,
        Enclosure(
            :armor_matrix,
            armor;
            primitive = Annulus(radius, armor_outer),
            fill = matrix
        ))
    push!(parts, Region(
        :armor_jacket,
        Shell(jacket_thickness),
        pp
    ))
    design = build(
        CableDesign,
        "525kV_1600mm2",
        Stack(parts)
    )
    earth = homogeneous(
        rho = earth_rho, eps_r = earth_eps_r, mu_r = 1.0
    )
    connections = [Dict(:core => 3index - 2, :sheath => 3index - 1, :armor => 3index)
                   for index in 1:2]
    system = build(
        LineCableSystem,
        fill(design, 2),
        [Pose2(cable_x[index], cable_y) for index in 1:2];
        connections,
        system_id = "cable_525kv_1600mm2_bipole",
        line_length = line_length
    )
    return LineParametersProblem(system; temperature, earth_props=earth, frequencies)
end
