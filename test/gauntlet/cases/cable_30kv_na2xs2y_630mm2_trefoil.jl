# Translated from the legacy dataset/trefoil_NA2XS2Y_630.jl declaration.
# Construction dimensions are in metres. The legacy five-layer core means
# one centre strand plus four 6k courses: 61 aluminium wires in total.
# Preserve strand dimensions rather than imposing the nominal 630 mm² area
# or the tabulated conductor/insulation/overall diameters of 30.2/47.7/56 mm.
# Legacy trefoil d=0 selected the minimum non-overlapping spacing: touching
# jackets. Gauntlet retains all six terminals for primitive-matrix comparison;
# the source instead connected cores to phases 1/2/3 and grounded the screens.
# Normalize the source's 50:50:2500 sweep to the standard trefoil gauntlet
# range: 101 logarithmic samples from 0.1 Hz to 1 MHz.
# Unspecified interstitial space is explicit air, as in the other wire cases.
case_definition(
    :cable_30kv_na2xs2y_630mm2_trefoil,
    (
        core_courses = case_parameter(
            :core_courses, 4; tags = (:topology, :cable_layer)
        ),
        core_strand_diameter = case_parameter(
            :core_strand_diameter, 3.65e-3; tags = (:geometry, :cable_layer)
        ),
        core_lay_ratio = case_parameter(
            :core_lay_ratio, 15.0; tags = (:geometry, :cable_layer)
        ),
        inner_semicon_thickness = case_parameter(
            :inner_semicon_thickness, 0.50e-3; tags = (:geometry, :cable_layer)
        ),
        insulation_thickness = case_parameter(
            :insulation_thickness, 8.00e-3; tags = (:geometry, :cable_layer)
        ),
        outer_semicon_thickness = case_parameter(
            :outer_semicon_thickness, 0.40e-3; tags = (:geometry, :cable_layer)
        ),
        screen_wires = case_parameter(
            :screen_wires, 60; tags = (:topology, :cable_layer)
        ),
        screen_wire_diameter = case_parameter(
            :screen_wire_diameter, 0.85e-3; tags = (:geometry, :cable_layer)
        ),
        screen_lay_ratio = case_parameter(
            :screen_lay_ratio, 15.0; tags = (:geometry, :cable_layer)
        ),
        copper_tape_thickness = case_parameter(
            :copper_tape_thickness, 0.10e-3; tags = (:geometry, :cable_layer)
        ),
        copper_tape_width = case_parameter(
            :copper_tape_width, 15.0e-3; tags = (:geometry, :cable_layer)
        ),
        jacket_thickness = case_parameter(
            :jacket_thickness, 2.70e-3; tags = (:geometry, :cable_layer)
        ),
        formation_x = case_parameter(:formation_x, 0.0; tags = (:geometry, :system)),
        formation_y = case_parameter(:formation_y, -1.0; tags = (:geometry, :system)),
        formation_clearance = case_parameter(
            :formation_clearance, 0.0; tags = (:geometry, :system)
        ),
        line_length = case_parameter(:line_length, 1000.0; tags = (:operation, :system)),
        temperature = case_parameter(:temperature, 20.0; tags = (:operation, :system)),
        earth_rho = case_parameter(:earth_rho, 100.0; tags = (:material, :earth)),
        earth_eps_r = case_parameter(:earth_eps_r, 10.0; tags = (:material, :earth)),
        frequencies = case_parameter(
            :frequencies, collect(10.0 .^ range(-1, stop = 6, length = 101));
            tags = (:operation, :frequency)
        )
    ),
    [
        "cable:1:core", "cable:1:sheath",
        "cable:2:core", "cable:2:sheath",
        "cable:3:core", "cable:3:sheath"
    ];
    description = "18/30 kV NA2XS2Y 630 mm² Al cables in touching trefoil"
) do p
    materials = LineCableModels.MaterialsLibrary(add_defaults = true)
    aluminum = LineCableModels.Material(materials, :aluminum)
    copper = LineCableModels.Material(materials, :copper)
    xlpe = LineCableModels.Material(materials, :xlpe)
    semicon1 = LineCableModels.Material(materials, :semicon1)
    semicon2 = LineCableModels.Material(materials, :semicon2)
    pe = LineCableModels.Material(materials, :pe)
    air = LineCableModels.Material(materials, :air)

    wire_radius = p.core_strand_diameter / 2
    core_radius = (2p.core_courses + 1) * wire_radius
    core = LineCableModels.terminal(
        :core,
        LineCableModels.stranded(
            aluminum;
            shape = LineCableModels.Disk(wire_radius),
            boundary = LineCableModels.Disk(core_radius),
            lay = LineCableModels.LayRatio(p.core_lay_ratio),
            fill = air
        ),
        LineCableModels.screen(
            semicon1; t = p.inner_semicon_thickness, tag = :inner_semiconductor
        ),
        LineCableModels.insulation(xlpe; t = p.insulation_thickness),
        LineCableModels.screen(
            semicon2; t = p.outer_semicon_thickness, tag = :outer_semiconductor
        )
    )

    screen_inner = core_radius + p.inner_semicon_thickness +
                   p.insulation_thickness + p.outer_semicon_thickness
    screen_radius = p.screen_wire_diameter / 2
    screen_outer = screen_inner + p.screen_wire_diameter
    screen = LineCableModels.Group(
        :sheath,
        LineCableModels.Region(
            :sheath_wires, LineCableModels.Disk(screen_radius), copper
        );
        pattern = LineCableModels.Ring(p.screen_wires; r = screen_inner + screen_radius),
        path = LineCableModels.Helix(LineCableModels.LayRatio(p.screen_lay_ratio))
    )
    tape_outer = screen_outer + p.copper_tape_thickness
    tape = LineCableModels.Group(
        :sheath,
        LineCableModels.Region(
            :sheath_copper_tape,
            LineCableModels.Rectangle(p.copper_tape_width, p.copper_tape_thickness),
            copper
        );
        pattern = LineCableModels.Ring(1; r = (screen_outer + tape_outer) / 2),
        path = LineCableModels.Helix(LineCableModels.LayRatio(p.screen_lay_ratio))
    )
    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "na2xs2y_1x630_35_18_30kv",
        LineCableModels.Stack(
            core,
            LineCableModels.Enclosure(
                :screen_matrix,
                LineCableModels.Stack(screen, tape);
                primitive = LineCableModels.Annulus(screen_inner, tape_outer),
                fill = air
            ),
            LineCableModels.jacket(pe; t = p.jacket_thickness)
        );
        # Datasheet metadata retains the source units: kV, mm², Ω/km,
        # μF/km and mH/km. These values do not override the geometry.
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
    placements = LineCableModels.trefoil(
        design;
        center = LineCableModels.at(p.formation_x, p.formation_y),
        spacing = 2 * LineCableModels.outer_radius(design) + p.formation_clearance,
        connections = (core = (1, 3, 5), sheath = (2, 4, 6))
    )
    system = LineCableModels.build(
        LineCableModels.LineCableSystem, placements;
        system_id = "cable_30kv_na2xs2y_630mm2_trefoil",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = LineCableModels.homogeneous(
            rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
        ),
        frequencies = p.frequencies
    )
end
