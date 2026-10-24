function known_case_problem(::Val{:cable_132kv_630mm2_flathor};
        frequencies=10.0 .^ range(-1, 7; length=101))
    library = MaterialsLibrary(add_defaults=true)
    core_layers = 5
    core_strand_diameter = 0.0036648500338805653
    core_lay_ratio = 13.0
    screen_wires = 19
    screen_wire_diameter = 0.0025881867280128637
    screen_lay_ratio = 10.0
    semicon_tape_thickness = 0.0003
    inner_semicon_thickness = 0.000768
    insulation_thickness = 0.016
    outer_semicon_thickness = 0.000472
    copper_tape_thickness = 0.0001
    copper_tape_width = 0.01
    copper_tape_lay_ratio = 10.0
    water_blocking_thickness = 0.00094
    aluminum_tape_thickness = 0.00015
    jacket_thickness = 0.00505
    cable_x = (-0.5, 0.0, 0.5)
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
    matrix = Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )
    wire_radius = core_strand_diameter / 2
    radius = (2core_layers - 1) * wire_radius
    stranded_core = stranded(
        copper;
        shape = Disk(wire_radius),
        lay = LayRatio(core_lay_ratio),
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
        path = Helix(LayRatio(screen_lay_ratio))
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
    radius = aluminum_outer
    push!(parts,
        Region(
            :jacket_insulation,
            Shell(jacket_thickness),
            pe
        ))
    design = build(
        CableDesign,
        "132kV_630mm2",
        Stack(parts)
    )
    earth = homogeneous(
        rho = earth_rho, eps_r = earth_eps_r, mu_r = 1.0
    )
    designs = fill(design, 3)
    positions = [Pose2(cable_x[index], cable_y) for index in 1:3]
    connections = [Dict(:core => 3index - 2, :sheath => 3index - 1, :jacket => 3index)
                   for index in 1:3]
    system = build(
        LineCableSystem,
        designs,
        positions;
        connections,
        system_id = "cable_132kv_630mm2_flathor",
        line_length = line_length
    )
    return LineParametersProblem(system; temperature, earth_props=earth, frequencies)
end
