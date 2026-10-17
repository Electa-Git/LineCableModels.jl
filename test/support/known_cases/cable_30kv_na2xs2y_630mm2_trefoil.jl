function known_case_problem(::Val{:cable_30kv_na2xs2y_630mm2_trefoil};
        frequencies=10.0 .^ range(-1, 7; length=101))
    library = MaterialsLibrary(add_defaults=true)
    core_courses = 4
    core_strand_diameter = 0.00365
    core_lay_ratio = 15.0
    inner_semicon_thickness = 0.0005
    insulation_thickness = 0.008
    outer_semicon_thickness = 0.0004
    screen_wires = 60
    screen_wire_diameter = 0.00085
    screen_lay_ratio = 15.0
    copper_tape_thickness = 0.0001
    copper_tape_width = 0.015
    jacket_thickness = 0.0027
    formation_x = 0.0
    formation_y = -1.0
    formation_clearance = 0.0
    line_length = 1000.0
    temperature = 20.0
    earth_rho = 100.0
    earth_eps_r = 10.0

    aluminum = Material(library, :aluminum)
    copper = Material(library, :copper)
    xlpe = Material(library, :xlpe)
    semicon1 = Material(library, :semicon1)
    semicon2 = Material(library, :semicon2)
    pe = Material(library, :pe)
    air = Material(library, :air)

    wire_radius = core_strand_diameter / 2
    core_radius = (2core_courses + 1) * wire_radius
    core = terminal(
        :core,
        stranded(
            aluminum;
            shape = Disk(wire_radius),
            boundary = Disk(core_radius),
            lay = LayRatio(core_lay_ratio),
            fill = air
        ),
        screen(
            semicon1; t = inner_semicon_thickness, tag = :inner_semiconductor
        ),
        LineCableModels.insulation(xlpe; t = insulation_thickness),
        screen(
            semicon2; t = outer_semicon_thickness, tag = :outer_semiconductor
        )
    )

    screen_inner = core_radius + inner_semicon_thickness +
                   insulation_thickness + outer_semicon_thickness
    screen_radius = screen_wire_diameter / 2
    screen_outer = screen_inner + screen_wire_diameter
    screen_part = Group(
        :sheath,
        Region(
            :sheath_wires, Disk(screen_radius), copper
        );
        pattern = Ring(screen_wires; r = screen_inner + screen_radius),
        path = Helix(LayRatio(screen_lay_ratio))
    )
    tape_outer = screen_outer + copper_tape_thickness
    tape = Group(
        :sheath,
        Region(
            :sheath_copper_tape,
            Rectangle(copper_tape_width, copper_tape_thickness),
            copper
        );
        pattern = Ring(1; r = (screen_outer + tape_outer) / 2),
        path = Helix(LayRatio(screen_lay_ratio))
    )
    design = build(
        CableDesign,
        "na2xs2y_1x630_35_18_30kv",
        Stack(
            core,
            Enclosure(
                :screen_matrix,
                Stack(screen_part, tape);
                primitive = Annulus(screen_inner, tape_outer),
                fill = air
            ),
            jacket(pe; t = jacket_thickness)
        );
        nominal_data = (
            designation_code = "NA2XS2Y 1x630/35 18/30kV",
            U0 = 18.0,
            U = 30.0,
            conductor_cross_section = 630.0,
            screen_cross_section = 35.0,
            resistance = 0.0469,
            capacitance = 0.29,
            inductance = 0.4775
        )
    )
    placements = trefoil(
        design;
        center = at(formation_x, formation_y),
        spacing = 2 * outer_radius(design) + formation_clearance,
        connections = (core = (1, 3, 5), sheath = (2, 4, 6))
    )
    system = build(
        LineCableSystem, placements;
        system_id = "cable_30kv_na2xs2y_630mm2_trefoil",
        line_length = line_length
    )
    earth = homogeneous(
            rho = earth_rho, eps_r = earth_eps_r, mu_r = 1.0
        )
    return LineParametersProblem(system; temperature, earth_props=earth, frequencies)
end
