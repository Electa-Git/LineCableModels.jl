function known_case_problem(::Val{:cable_18kv_1000mm2_trefoil_homogenized};
        frequencies=10.0 .^ range(-1, 7; length=101))
    library = MaterialsLibrary(add_defaults=true)
    core_strand_diameter = 0.00419
    core_ring_lay_ratio_1 = 15.0
    core_ring_lay_ratio_2 = 13.5
    core_ring_lay_ratio_3 = 12.5
    core_ring_lay_ratio_4 = 11.0
    core_ring_wire_counts = (1, 6, 12, 18, 24)
    screen_wires = 49
    screen_wire_diameter = 0.00094
    screen_wire_lay_ratio = 10.0
    semicon_tape_thickness = 0.0003
    inner_semicon_thickness = 0.000768
    insulation_thickness = 0.0083
    outer_semicon_thickness = 0.000472
    copper_tape_thickness = 0.0001
    copper_tape_width = 0.01
    copper_tape_lay_ratio = 10.0
    water_blocking_thickness = 0.00094
    aluminum_tape_thickness = 0.00015
    pe_face_thickness = 5.0e-5
    jacket_thickness = 0.0034
    formation_x = 0.0
    formation_y = -1.0
    formation_clearance_ratio = 0.2
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
    matrix = Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )

    counts = core_ring_wire_counts
    lay_ratios = LayRatio(
        core_ring_lay_ratio_1,
        core_ring_lay_ratio_2,
        core_ring_lay_ratio_3,
        core_ring_lay_ratio_4
    )
    wire_radius = core_strand_diameter / 2
    radius = (2length(counts) - 1) * wire_radius
    stranded_core = stranded(
        aluminum;
        shape = Disk(wire_radius),
        lay = lay_ratios,
        boundary = Disk(radius),
        fill = matrix
    )
    core = stranded_core
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

    screen_radius = screen_wire_diameter / 2
    screen_inner = radius
    screen_centre = radius + screen_radius
    screen_part = Group(
        :sheath,
        Region(
            :sheath_wires, Disk(screen_radius), copper
        );
        pattern = Ring(screen_wires; r = screen_centre),
        path = Helix(LayRatio(screen_wire_lay_ratio))
    )
    radius += 2screen_radius
    tape_outer = radius + copper_tape_thickness
    tape = Group(
        :sheath,
        Region(
            :sheath_copper_tape,
            Rectangle(
                copper_tape_width,
                copper_tape_thickness
            ),
            copper
        );
        pattern = Ring(1; r = (radius + tape_outer) / 2),
        path = Helix(LayRatio(copper_tape_lay_ratio))
    )
    push!(parts,
        Enclosure(
            :screen_matrix,
            Stack(screen_part, tape);
            primitive = Annulus(screen_inner, tape_outer),
            fill = matrix
        ))
    radius = tape_outer
    push!(parts,
        Region(
            :sheath_water_blocking,
            Shell(water_blocking_thickness),
            polyacrylate
        ))
    radius += water_blocking_thickness

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
            :jacket_pe_face, Shell(pe_face_thickness), pe
        ))
    push!(parts,
        Region(
            :jacket_insulation, Shell(jacket_thickness), pe
        ))
    design = build(
        CableDesign,
        "18kV_1000mm2",
        Stack(parts)
    )
    formation_spacing = (2 + formation_clearance_ratio) *
                        outer_radius(design)
    earth = homogeneous(
        rho = earth_rho, eps_r = earth_eps_r, mu_r = 1.0
    )
    placements = trefoil(
        design;
        center = at(formation_x, formation_y),
        spacing = formation_spacing,
        connections = (
            core = (1, 4, 7),
            sheath = (2, 5, 8),
            jacket = (3, 6, 9)
        )
    )
    system = build(
        LineCableSystem,
        placements;
        system_id = "cable_18kv_1000mm2_trefoil",
        line_length = line_length
    )
    source = system
    designs = map(homogenize, source.designs)
    system = build(
        LineCableSystem,
        designs,
        source.positions;
        connections = source.connections,
        environment = source.environment,
        system_id = "cable_18kv_1000mm2_trefoil_homogenized",
        line_length = source.line_length
    )
    return LineParametersProblem(system; temperature, earth_props=earth, frequencies)
end
