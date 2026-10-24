function known_case_problem(::Val{:cable_220kv_milliken_1x2500_252_trefoil};
        frequencies=10.0 .^ range(-1, 7; length=101))
    library = MaterialsLibrary(add_defaults=true)
    core_outer_radius = 0.03115
    core_wire_radius = 0.001475
    core_sectors = 6
    core_fillet_factor = 0.04
    binder_thickness = 0.0005
    inner_semicon_thickness = 0.0016
    insulation_thickness = 0.0203
    outer_semicon_thickness = 0.0013
    screen_inner_bedding_thickness = 0.0003
    screen_wires = 65
    screen_wire_radius = 0.00111
    screen_lay_ratio = 5.5012442232492
    screen_outer_bedding_thickness = 0.00073
    aluminum_foil_thickness = 0.0002
    jacket_thickness = 0.0048
    semicon_rho = 0.06
    semicon_eps_r = 1000.0
    xlpe_tan_delta = 0.001
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
    hdpe = pe
    xlpe = Material(
        kind = :insulator, rho = 1.97e14, eps_r = 2.5,
        tan_delta = xlpe_tan_delta
    )
    carbon_pe = Material(
        kind = :semicon, rho = semicon_rho, eps_r = semicon_eps_r
    )

    span = 2pi / core_sectors
    fillet = core_fillet_factor * core_outer_radius
    segment = Sector(
        span = span,
        r_base = fillet * (inv(sin(span / 2)) - 1),
        r_back = core_outer_radius,
        fillet = fillet
    )
    core = terminal(
        :core,
        milliken(
            aluminum;
            shape = Disk(core_wire_radius),
            segment,
            segments = core_sectors,
            fill = pe
        ),
        bedding(pe; t = binder_thickness, tag = :conductor_binder),
        screen(carbon_pe; t = inner_semicon_thickness, tag = :inner_semicon),
        LineCableModels.insulation(xlpe; t = insulation_thickness, tag = :xlpe),
        screen(carbon_pe; t = outer_semicon_thickness, tag = :outer_semicon),
        bedding(pe; t = screen_inner_bedding_thickness, tag = :screen_bedding_inner)
    )
    screen_inner = core_outer_radius + binder_thickness +
                   inner_semicon_thickness + insulation_thickness +
                   outer_semicon_thickness + screen_inner_bedding_thickness
    screen_outer = screen_inner + 2screen_wire_radius
    wires = Group(
        :sheath,
        Region(:copper_screen_wires, Disk(screen_wire_radius), copper);
        pattern = Ring(screen_wires; r = screen_inner +
                                                           screen_wire_radius),
        path = Helix(LayRatio(screen_lay_ratio))
    )
    metallic_screen = Stack(
        pipe(
            wires; shape = Annulus(screen_inner, screen_outer), fill = pe
        ),
        bedding(pe; t = screen_outer_bedding_thickness, tag = :screen_bedding_outer),
        terminal(
            :foil,
            sheath(aluminum; t = aluminum_foil_thickness, tag = :aluminium_foil),
            jacket(hdpe; t = jacket_thickness, tag = :hdpe_oversheath)
        )
    )
    design = build(
        CableDesign,
        "220kv_eaxecew_1x2500_252",
        Stack(core, metallic_screen)
    )
    spacing = (2 + formation_clearance_ratio) * outer_radius(design)
    placements = trefoil(
        design;
        center = at(formation_x, formation_y),
        spacing,
        connections = (core = (1, 4, 7), sheath = (2, 5, 8), foil = (3, 6, 9))
    )
    system = build(
        LineCableSystem, placements;
        system_id = "cable_220kv_milliken_1x2500_252_trefoil",
        line_length = line_length
    )
    earth = homogeneous(rho = earth_rho, eps_r = earth_eps_r)
    return LineParametersProblem(system; temperature, earth_props=earth, frequencies)
end
