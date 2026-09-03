case_definition(
    :cable_525kv_1600mm2_bipole,
    (
        strand_layers = case_parameter(
            :strand_layers, 6; tags = (:topology, :cable_layer)
        ),
        strands_per_layer = case_parameter(
            :strands_per_layer, 6; tags = (:topology, :cable_layer)
        ),
        strand_diameter = case_parameter(
            :strand_diameter, 3.6649e-3; tags = (:geometry, :cable_layer)
        ),
        core_lay_ratio = case_parameter(
            :core_lay_ratio, 11.0; tags = (:geometry, :cable_layer)
        ),
        inner_semicon_thickness = case_parameter(
            :inner_semicon_thickness, 2.0e-3; tags = (:geometry, :cable_layer)
        ),
        insulation_thickness = case_parameter(
            :insulation_thickness, 26.0e-3; tags = (:geometry, :cable_layer)
        ),
        outer_semicon_thickness = case_parameter(
            :outer_semicon_thickness, 1.8e-3; tags = (:geometry, :cable_layer)
        ),
        water_blocking_thickness = case_parameter(
            :water_blocking_thickness, 0.3e-3; tags = (:geometry, :cable_layer)
        ),
        lead_screen_thickness = case_parameter(
            :lead_screen_thickness, 3.3e-3; tags = (:geometry, :cable_layer)
        ),
        inner_sheath_thickness = case_parameter(
            :inner_sheath_thickness, 3.0e-3; tags = (:geometry, :cable_layer)
        ),
        bedding_thickness = case_parameter(
            :bedding_thickness, 3.0e-3; tags = (:geometry, :cable_layer)
        ),
        armor_wires = case_parameter(
            :armor_wires, 68; tags = (:topology, :cable_layer)
        ),
        armor_wire_diameter = case_parameter(
            :armor_wire_diameter, 5.827e-3; tags = (:geometry, :cable_layer)
        ),
        armor_lay_ratio = case_parameter(
            :armor_lay_ratio, 10.0; tags = (:geometry, :cable_layer)
        ),
        armor_packing_clearance_ratio = case_parameter(
            :armor_packing_clearance_ratio, 0.2;
            tags = (:constraint, :cable_layer)
        ),
        jacket_thickness = case_parameter(
            :jacket_thickness, 10.0e-3; tags = (:geometry, :cable_layer)
        ),
        cable_x = case_parameter(
            :cable_x, (-0.5, 0.5); tags = (:geometry, :system)
        ),
        cable_y = case_parameter(:cable_y, -1.0; tags = (:geometry, :system)),
        line_length = case_parameter(
            :line_length, 1000.0; tags = (:operation, :system)
        ),
        temperature = case_parameter(
            :temperature, 20.0; tags = (:operation, :system)
        ),
        earth_rho = case_parameter(:earth_rho, 100.0; tags = (:material, :earth)),
        earth_eps_r = case_parameter(
            :earth_eps_r, 10.0; tags = (:material, :earth)
        ),
        frequencies = case_parameter(
            :frequencies, collect(10.0 .^ range(0, stop = 6, length = 101));
            tags = (:operation, :frequency)
        )
    ),
    [
        "cable:1:core", "cable:1:sheath", "cable:1:armor",
        "cable:2:core", "cable:2:sheath", "cable:2:armor"
    ];
    description = "525 kV 1600 mm² armoured cable DC bipole"
) do p
    materials = LineCableModels.MaterialsLibrary(add_defaults = true)
    copper = LineCableModels.Material(materials, :copper)
    semicon1 = LineCableModels.Material(materials, :semicon1)
    semicon2 = LineCableModels.Material(materials, :semicon2)
    pe = LineCableModels.Material(materials, :pe)
    polyacrylate = LineCableModels.Material(materials, :polyacrylate)
    lead = LineCableModels.Material(materials, :lead)
    pp = LineCableModels.Material(materials, :pp)
    steel = LineCableModels.Material(materials, :steel)
    matrix = LineCableModels.Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )
    radius_before_bedding = (p.strand_layers + 0.5) * p.strand_diameter +
                            p.inner_semicon_thickness +
                            p.insulation_thickness +
                            p.outer_semicon_thickness +
                            p.water_blocking_thickness +
                            p.lead_screen_thickness +
                            p.inner_sheath_thickness
    armor_wire_radius = p.armor_wire_diameter / 2
    minimum_armor_radius = armor_wire_radius / sinpi(1 / p.armor_wires) - armor_wire_radius
    unbuffered_armor_radius = radius_before_bedding + p.bedding_thickness
    unbuffered_outer_radius = unbuffered_armor_radius + p.armor_wire_diameter +
                              p.jacket_thickness
    packing_buffer = p.armor_packing_clearance_ratio * unbuffered_outer_radius
    buffered_armor_radius = unbuffered_armor_radius + packing_buffer
    # Fixed-count armor wires cannot occupy a ring smaller than their packing
    # radius. A fixed design clearance keeps ordinary 10% manufacturing draws
    # away from that boundary; the compliant bedding still absorbs any extreme
    # positive shortfall instead of constructing an unsupported wire packing.
    packing_shortfall = max(
        minimum_armor_radius - buffered_armor_radius,
        zero(minimum_armor_radius - buffered_armor_radius)
    )
    bedding_thickness = p.bedding_thickness + packing_buffer + packing_shortfall
    strand_radius = p.strand_diameter / 2
    radius = (2p.strand_layers + 1) * strand_radius
    stranded_core = LineCableModels.stranded(
        copper;
        shape = LineCableModels.Disk(strand_radius),
        layers = p.strand_layers,
        n = p.strands_per_layer,
        lay = LineCableModels.LayRatio(p.core_lay_ratio),
        boundary = LineCableModels.Disk(radius)
    )
    core = LineCableModels.Enclosure(
        :core_matrix,
        stranded_core;
        primitive = LineCableModels.Disk(radius),
        fill = matrix
    )
    parts = LineCableModels.AbstractCablePart[core]
    for (tag, thickness, material) in (
        (:core_semicon_inner, p.inner_semicon_thickness, semicon1),
        (:core_insulation, p.insulation_thickness, pe),
        (:core_semicon_outer, p.outer_semicon_thickness, semicon2),
        (:core_water_blocking, p.water_blocking_thickness, polyacrylate)
    )
        push!(parts, LineCableModels.Region(
            tag, LineCableModels.Shell(thickness), material
        ))
        radius += thickness
    end

    lead_outer = radius + p.lead_screen_thickness
    push!(parts,
        LineCableModels.Group(
            :sheath,
            LineCableModels.Region(
                :sheath_lead_screen,
                LineCableModels.Annulus(radius, lead_outer),
                lead
            )
        ))
    radius = lead_outer
    push!(parts,
        LineCableModels.Region(
            :sheath_inner,
            LineCableModels.Shell(p.inner_sheath_thickness),
            pe
        ))
    radius += p.inner_sheath_thickness
    push!(parts,
        LineCableModels.Region(
            :sheath_bedding,
            LineCableModels.Shell(bedding_thickness),
            pp
        ))
    radius += bedding_thickness

    armor_outer = radius + 2armor_wire_radius
    armor_centre = radius + armor_wire_radius
    armor = LineCableModels.Group(
        :armor,
        LineCableModels.Region(
            :armor_wires, LineCableModels.Disk(armor_wire_radius), steel
        );
        pattern = LineCableModels.Ring(p.armor_wires; r = armor_centre),
        path = LineCableModels.Helix(LineCableModels.LayRatio(p.armor_lay_ratio))
    )
    push!(parts,
        LineCableModels.Enclosure(
            :armor_matrix,
            armor;
            primitive = LineCableModels.Annulus(radius, armor_outer),
            fill = matrix
        ))
    push!(parts, LineCableModels.Region(
        :armor_jacket,
        LineCableModels.Shell(p.jacket_thickness),
        pp
    ))
    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "525kV_1600mm2",
        LineCableModels.Stack(parts)
    )
    earth = LineCableModels.homogeneous(
        rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
    )
    connections = [Dict(:core => 3index - 2, :sheath => 3index - 1, :armor => 3index)
                   for index in 1:2]
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        fill(design, 2),
        [LineCableModels.Pose2(p.cable_x[index], p.cable_y) for index in 1:2];
        connections,
        system_id = "cable_525kv_1600mm2_bipole",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = earth,
        frequencies = p.frequencies
    )
end
