# Adapted from 220kv_eaxecew_1x2500_252.json, CableBuilderSpec v1,
# and its neighbouring SeraingThermal/specs/materials.json.
# Dimensions are in metres. The legacy sector count of 60 is replaced by the
# current inferred 61-strand inventory per segment and a tangent shared centre.
# The old dimensionless contour-rounding factor is interpreted as fillet/R;
# r_base is the minimum value keeping the filleted sectors inside their pitches.
# This is a current-API reconstruction, not a bit-exact legacy geometry import.
# The legacy screen group contains two metals separated by PE. Following the
# gauntlet convention, the copper screen and aluminium foil have separate ports;
# no bonding or reduction is baked into the case definition.
# Nominal 2500/252 mm² labels do not override the explicit source dimensions.
# Operating conditions follow the existing trefoil gauntlet, not the thermal job.
case_definition(
    :cable_220kv_eaxecew_1x2500_252_trefoil,
    (
        core_outer_radius = case_parameter(
            :core_outer_radius, 0.03115; tags = (:geometry, :cable_layer)
        ),
        core_wire_radius = case_parameter(
            :core_wire_radius, 0.001475; tags = (:geometry, :cable_layer)
        ),
        core_sectors = case_parameter(
            :core_sectors, 6; tags = (:topology, :cable_layer)
        ),
        core_fillet_factor = case_parameter(
            :core_fillet_factor, 0.04; tags = (:geometry, :cable_layer)
        ),
        binder_thickness = case_parameter(
            :binder_thickness, 0.0005; tags = (:geometry, :cable_layer)
        ),
        inner_semicon_thickness = case_parameter(
            :inner_semicon_thickness, 0.0016; tags = (:geometry, :cable_layer)
        ),
        insulation_thickness = case_parameter(
            :insulation_thickness, 0.0203; tags = (:geometry, :cable_layer)
        ),
        outer_semicon_thickness = case_parameter(
            :outer_semicon_thickness, 0.0013; tags = (:geometry, :cable_layer)
        ),
        screen_inner_bedding_thickness = case_parameter(
            :screen_inner_bedding_thickness, 0.0003; tags = (:geometry, :cable_layer)
        ),
        screen_wires = case_parameter(
            :screen_wires, 65; tags = (:topology, :cable_layer)
        ),
        screen_wire_radius = case_parameter(
            :screen_wire_radius, 0.00111; tags = (:geometry, :cable_layer)
        ),
        screen_lay_ratio = case_parameter(
            :screen_lay_ratio, 5.5012442232492; tags = (:geometry, :cable_layer)
        ),
        screen_outer_bedding_thickness = case_parameter(
            :screen_outer_bedding_thickness, 0.00073; tags = (:geometry, :cable_layer)
        ),
        aluminum_foil_thickness = case_parameter(
            :aluminum_foil_thickness, 0.0002; tags = (:geometry, :cable_layer)
        ),
        jacket_thickness = case_parameter(
            :jacket_thickness, 0.0048; tags = (:geometry, :cable_layer)
        ),
        semicon_rho = case_parameter(
            :semicon_rho, 0.06; tags = (:material, :cable_layer)
        ),
        semicon_eps_r = case_parameter(
            :semicon_eps_r, 1000.0; tags = (:material, :cable_layer)
        ),
        xlpe_tan_delta = case_parameter(
            :xlpe_tan_delta, 0.001; tags = (:material, :cable_layer)
        ),
        formation_x = case_parameter(:formation_x, 0.0; tags = (:geometry, :system)),
        formation_y = case_parameter(:formation_y, -1.0; tags = (:geometry, :system)),
        formation_clearance_ratio = case_parameter(
            :formation_clearance_ratio, 0.2; tags = (:geometry, :system)
        ),
        line_length = case_parameter(:line_length, 1000.0; tags = (:operation, :system)),
        temperature = case_parameter(:temperature, 20.0; tags = (:operation, :system)),
        earth_rho = case_parameter(:earth_rho, 100.0; tags = (:material, :earth)),
        earth_eps_r = case_parameter(:earth_eps_r, 10.0; tags = (:material, :earth)),
        frequencies = case_parameter(
            :frequencies, collect(10.0 .^ range(0, stop = 6, length = 101));
            tags = (:operation, :frequency)
        )
    ),
    [
        "cable:1:core", "cable:1:sheath", "cable:1:foil",
        "cable:2:core", "cable:2:sheath", "cable:2:foil",
        "cable:3:core", "cable:3:sheath", "cable:3:foil"
    ];
    description = "220 kV EAXeCeW 2500 mm² Al / 252 mm² Cu cables in trefoil"
) do p
    materials = LineCableModels.MaterialsLibrary(add_defaults = true)
    aluminum = LineCableModels.Material(materials, :aluminum)
    copper = LineCableModels.Material(materials, :copper)
    pe = LineCableModels.Material(materials, :pe)
    # HDPE and PE have identical electrical properties in the source library.
    hdpe = pe
    xlpe = LineCableModels.Material(
        kind = :insulator, rho = 1.97e14, eps_r = 2.5,
        tan_delta = p.xlpe_tan_delta
    )
    carbon_pe = LineCableModels.Material(
        kind = :semicon, rho = p.semicon_rho, eps_r = p.semicon_eps_r
    )

    span = 2pi / p.core_sectors
    fillet = p.core_fillet_factor * p.core_outer_radius
    segment = LineCableModels.Sector(
        span = span,
        r_base = fillet * (inv(sin(span / 2)) - 1),
        r_back = p.core_outer_radius,
        fillet = fillet
    )
    core = LineCableModels.terminal(
        :core,
        LineCableModels.milliken(
            aluminum;
            shape = LineCableModels.Disk(p.core_wire_radius),
            segment,
            segments = p.core_sectors,
            fill = pe
        ),
        LineCableModels.bedding(pe; t = p.binder_thickness, tag = :conductor_binder),
        LineCableModels.screen(carbon_pe; t = p.inner_semicon_thickness, tag = :inner_semicon),
        LineCableModels.insulation(xlpe; t = p.insulation_thickness, tag = :xlpe),
        LineCableModels.screen(carbon_pe; t = p.outer_semicon_thickness, tag = :outer_semicon),
        LineCableModels.bedding(pe; t = p.screen_inner_bedding_thickness, tag = :screen_bedding_inner)
    )
    screen_inner = p.core_outer_radius + p.binder_thickness +
                   p.inner_semicon_thickness + p.insulation_thickness +
                   p.outer_semicon_thickness + p.screen_inner_bedding_thickness
    screen_outer = screen_inner + 2p.screen_wire_radius
    wires = LineCableModels.Group(
        :sheath,
        LineCableModels.Region(:copper_screen_wires, LineCableModels.Disk(p.screen_wire_radius), copper);
        pattern = LineCableModels.Ring(p.screen_wires; r = screen_inner + p.screen_wire_radius),
        path = LineCableModels.Helix(LineCableModels.LayRatio(p.screen_lay_ratio))
    )
    metallic_screen = LineCableModels.Stack(
        LineCableModels.pipe(
            wires; shape = LineCableModels.Annulus(screen_inner, screen_outer), fill = pe
        ),
        LineCableModels.bedding(pe; t = p.screen_outer_bedding_thickness, tag = :screen_bedding_outer),
        LineCableModels.terminal(
            :foil,
            LineCableModels.sheath(aluminum; t = p.aluminum_foil_thickness, tag = :aluminium_foil),
            LineCableModels.jacket(hdpe; t = p.jacket_thickness, tag = :hdpe_oversheath)
        )
    )
    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "220kv_eaxecew_1x2500_252",
        LineCableModels.Stack(core, metallic_screen)
    )
    spacing = (2 + p.formation_clearance_ratio) * LineCableModels.outer_radius(design)
    placements = LineCableModels.trefoil(
        design;
        center = LineCableModels.at(p.formation_x, p.formation_y),
        spacing,
        connections = (core = (1, 4, 7), sheath = (2, 5, 8), foil = (3, 6, 9))
    )
    system = LineCableModels.build(
        LineCableModels.LineCableSystem, placements;
        system_id = "cable_220kv_eaxecew_1x2500_252_trefoil",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = LineCableModels.homogeneous(rho = p.earth_rho, eps_r = p.earth_eps_r),
        frequencies = p.frequencies
    )
end
