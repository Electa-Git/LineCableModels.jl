function known_case_problem(::Val{:cable_640kv_2000mm2_bipole};
        frequencies=10.0 .^ range(-1, 7; length=101))
    library = MaterialsLibrary(add_defaults=true)
    core_strand_radius = 0.0032718537481013924
    ring_counts = (6, 12, 18, 24)
    core_lay_ratio = 13.0
    semicon_tape_thickness = 0.0003
    inner_semicon_thickness = 0.000768
    insulation_thickness = 0.032
    outer_semicon_thickness = 0.000472
    lead_screen_thickness = 0.0033
    inner_sheath_thickness = 0.003
    aluminum_tape_thickness = 0.00015
    jacket_thickness = 0.00605
    cable_x = (-0.5, 0.5)
    cable_y = -1.0
    line_length = 1000.0
    temperature = 20.0
    earth_rho = 100.0
    earth_eps_r = 10.0

    aluminum = Material(library, :aluminum)
    copper = Material(library, :copper)
    pe = Material(library, :pe)
    xlpe = Material(library, :xlpe)
    semicon1 = Material(library, :semicon1)
    semicon2 = Material(library, :semicon2)
    polyacrylate = Material(library, :polyacrylate)
    lead = Material(library, :lead)
    matrix = Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )
    strand_count = 1 + sum(ring_counts)
    radius = sqrt(strand_count) * core_strand_radius
    core = stranded(
        copper;
        shape = Disk(core_strand_radius),
        lay = LayRatio(core_lay_ratio),
        compact = true,
        boundary = Disk(radius),
        fill = matrix
    )
    parts = AbstractCablePart[core]
    for (tag, thickness, material) in (
        (:core_semicon_tape_inner, semicon_tape_thickness, polyacrylate),
        (:core_semicon_inner, inner_semicon_thickness, semicon1),
        (:core_insulation, insulation_thickness, xlpe),
        (:core_semicon_outer, outer_semicon_thickness, semicon2),
        (:core_semicon_tape_outer, semicon_tape_thickness, polyacrylate)
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

    aluminum_outer = radius + aluminum_tape_thickness
    push!(parts,
        Group(
            :jacket,
            Region(
                :jacket_aluminum_tape,
                Annulus(radius, aluminum_outer),
                aluminum
            )
        ))
    push!(parts,
        Region(
            :jacket_insulation,
            Shell(jacket_thickness),
            pe
        ))
    design = build(
        CableDesign,
        "640kV_2000mm2",
        Stack(parts)
    )
    earth = homogeneous(
        rho = earth_rho, eps_r = earth_eps_r, mu_r = 1.0
    )
    connections = [Dict(:core => 3index - 2, :sheath => 3index - 1, :jacket => 3index)
                   for index in 1:2]
    system = build(
        LineCableSystem,
        fill(design, 2),
        [Pose2(cable_x[index], cable_y) for index in 1:2];
        connections,
        system_id = "cable_640kv_2000mm2_bipole",
        line_length = line_length
    )
    return LineParametersProblem(system; temperature, earth_props=earth, frequencies)
end
