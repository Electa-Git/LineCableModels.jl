case_definition(
    :cable_380kv_armoured_ac_flat,
    (
        strand_layers = case_parameter(
            :strand_layers, 6; tags = (:topology, :cable_layer)
        ),
        strand_diameter = case_parameter(
            :strand_diameter, 3.66e-3; tags = (:geometry, :cable_layer)
        ),
        core_lay_ratio = case_parameter(
            :core_lay_ratio, 15.0; tags = (:geometry, :cable_layer)
        ),
        inner_semicon_thickness = case_parameter(
            :inner_semicon_thickness, 2.0e-3; tags = (:geometry, :cable_layer)
        ),
        insulation_thickness = case_parameter(
            :insulation_thickness, 20.0e-3; tags = (:geometry, :cable_layer)
        ),
        outer_semicon_thickness = case_parameter(
            :outer_semicon_thickness, 1.5e-3; tags = (:geometry, :cable_layer)
        ),
        water_blocking_thickness = case_parameter(
            :water_blocking_thickness, 0.3e-3; tags = (:geometry, :cable_layer)
        ),
        lead_sheath_thickness = case_parameter(
            :lead_sheath_thickness, 3.3e-3; tags = (:geometry, :cable_layer)
        ),
        inner_pe_thickness = case_parameter(
            :inner_pe_thickness, 3.0e-3; tags = (:geometry, :cable_layer)
        ),
        bedding_thickness = case_parameter(
            :bedding_thickness, 3.0e-3; tags = (:geometry, :cable_layer)
        ),
        armor_wires = case_parameter(
            :armor_wires, 68; tags = (:topology, :cable_layer)
        ),
        armor_wire_diameter = case_parameter(
            :armor_wire_diameter, 5.82e-3; tags = (:geometry, :cable_layer)
        ),
        armor_lay_ratio = case_parameter(
            :armor_lay_ratio, 15.0; tags = (:geometry, :cable_layer)
        ),
        jacket_thickness = case_parameter(
            :jacket_thickness, 10.0e-3; tags = (:geometry, :cable_layer)
        ),
        core_rho = case_parameter(
            :core_rho, 2.3853e-8; tags = (:material, :cable_layer)
        ),
        core_mu_r = case_parameter(
            :core_mu_r, 1.01663; tags = (:material, :cable_layer)
        ),
        xlpe_rho = case_parameter(
            :xlpe_rho, 1.0e14; tags = (:material, :cable_layer)
        ),
        xlpe_eps_r = case_parameter(
            :xlpe_eps_r, 2.77781; tags = (:material, :cable_layer)
        ),
        xlpe_mu_r = case_parameter(
            :xlpe_mu_r, 1.59891; tags = (:material, :cable_layer)
        ),
        lead_rho = case_parameter(
            :lead_rho, 2.14e-7; tags = (:material, :cable_layer)
        ),
        lead_mu_r = case_parameter(
            :lead_mu_r, 0.999983; tags = (:material, :cable_layer)
        ),
        pe_rho = case_parameter(
            :pe_rho, 1.97e14; tags = (:material, :cable_layer)
        ),
        pe_eps_r = case_parameter(
            :pe_eps_r, 2.51862; tags = (:material, :cable_layer)
        ),
        pp_rho = case_parameter(
            :pp_rho, 1.0e15; tags = (:material, :cable_layer)
        ),
        pp_eps_r = case_parameter(
            :pp_eps_r, 2.8; tags = (:material, :cable_layer)
        ),
        pp_mu_r = case_parameter(
            :pp_mu_r, 1.1263; tags = (:material, :cable_layer)
        ),
        steel_rho = case_parameter(
            :steel_rho, 1.7475e-7; tags = (:material, :cable_layer)
        ),
        steel_mu_r = case_parameter(
            :steel_mu_r, 36.6326; tags = (:material, :cable_layer)
        ),
        cable_x = case_parameter(
            :cable_x, (-1.0, 0.0, 1.0); tags = (:geometry, :system)
        ),
        cable_y = case_parameter(:cable_y, -1.0; tags = (:geometry, :system)),
        line_length = case_parameter(
            :line_length, 20_000.0; tags = (:operation, :system)
        ),
        voltage = case_parameter(:voltage, 380.0; tags = (:operation, :system)),
        temperature = case_parameter(
            :temperature, 70.0; tags = (:operation, :system)
        ),
        earth_rho = case_parameter(
            :earth_rho, 1_000.0; tags = (:material, :earth)
        ),
        earth_eps_r = case_parameter(
            :earth_eps_r, 10.0; tags = (:material, :earth)
        ),
        frequencies = case_parameter(
            :frequencies,
            sort!(unique!(vcat(
                collect(10.0 .^ range(-3, stop = 6, length = 226)),
                50.0
            )));
            tags = (:operation, :frequency)
        )
    ),
    [
        "cable:1:core", "cable:1:sheath", "cable:1:armor",
        "cable:2:core", "cable:2:sheath", "cable:2:armor",
        "cable:3:core", "cable:3:sheath", "cable:3:armor"
    ];
    description = "380 kV armoured cables in AC flat formation"
) do p
    core = LineCableModels.Material(
        :conductor, p.core_rho, 1.0, p.core_mu_r, 20.0, 0.00393
    )
    semicon_inner = LineCableModels.Material(:semicon, 1_000.0, 1_000.0)
    semicon_outer = LineCableModels.Material(:semicon, 500.0, 1_000.0)
    xlpe = LineCableModels.Material(
        :insulator, p.xlpe_rho, p.xlpe_eps_r, p.xlpe_mu_r
    )
    lead = LineCableModels.Material(
        :conductor, p.lead_rho, 1.0, p.lead_mu_r, 20.0, 0.004
    )
    pe = LineCableModels.Material(:insulator, p.pe_rho, p.pe_eps_r)
    pp = LineCableModels.Material(
        :insulator, p.pp_rho, p.pp_eps_r, p.pp_mu_r
    )
    steel = LineCableModels.Material(
        :conductor, p.steel_rho, 1.0, p.steel_mu_r, 20.0, 0.0045
    )
    matrix = LineCableModels.Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )

    strand_radius = p.strand_diameter / 2
    strand_count = 1 + 3p.strand_layers * (p.strand_layers + 1)
    radius = sqrt(strand_count) * strand_radius
    packed_core = LineCableModels.stranded(
        core;
        shape = LineCableModels.Disk(strand_radius),
        layers = p.strand_layers,
        n = 6,
        lay = LineCableModels.LayRatio(p.core_lay_ratio),
        compact = LineCableModels.FillFactor(1),
        boundary = LineCableModels.Disk(radius)
    )
    parts = LineCableModels.AbstractCablePart[packed_core]
    for (tag, layer_thickness, material) in (
        (:core_semicon_inner, p.inner_semicon_thickness, semicon_inner),
        (:core_insulation, p.insulation_thickness, xlpe),
        (:core_semicon_outer, p.outer_semicon_thickness, semicon_outer),
        (:core_water_blocking, p.water_blocking_thickness, semicon_outer)
    )
        push!(parts,
            LineCableModels.Region(tag, LineCableModels.Shell(layer_thickness), material))
        radius += layer_thickness
    end
    lead_outer = radius + p.lead_sheath_thickness
    push!(parts,
        LineCableModels.Group(
            :sheath,
            LineCableModels.Region(
                :sheath_lead,
                LineCableModels.Annulus(radius, lead_outer),
                lead
            )
        ))
    radius = lead_outer
    push!(parts,
        LineCableModels.Region(
            :sheath_inner_pe,
            LineCableModels.Shell(p.inner_pe_thickness),
            pe
        ))
    radius += p.inner_pe_thickness

    armor_wire_radius = p.armor_wire_diameter / 2
    minimum_armor_inner_radius = armor_wire_radius / sinpi(1 / p.armor_wires) -
                                 armor_wire_radius
    nominal_armor_inner_radius = radius + p.bedding_thickness
    packing_shortfall = max(
        minimum_armor_inner_radius - nominal_armor_inner_radius,
        zero(minimum_armor_inner_radius)
    )
    effective_bedding_thickness = p.bedding_thickness + packing_shortfall
    push!(parts,
        LineCableModels.Region(
            :sheath_bedding,
            LineCableModels.Shell(effective_bedding_thickness),
            pp
        ))
    radius += effective_bedding_thickness
    armor_outer = radius + 2armor_wire_radius
    armor = LineCableModels.Group(
        :armor,
        LineCableModels.Region(
            :armor_wires,
            LineCableModels.Disk(armor_wire_radius),
            steel
        );
        pattern = LineCableModels.Ring(
            p.armor_wires; r = radius + armor_wire_radius
        ),
        path = LineCableModels.Helix(LineCableModels.LayRatio(p.armor_lay_ratio))
    )
    push!(parts,
        LineCableModels.Enclosure(
            :armor_matrix,
            armor;
            primitive = LineCableModels.Annulus(radius, armor_outer),
            fill = matrix
        ))
    radius = armor_outer
    push!(parts,
        LineCableModels.Region(
            :armor_jacket,
            LineCableModels.Shell(p.jacket_thickness),
            pp
        ))

    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "cable_380kv_armoured",
        LineCableModels.Stack(parts);
        nominal_data = (
            designation_code = "cable_380kv_armoured",
            U0 = p.voltage / sqrt(3),
            U = p.voltage,
            conductor_cross_section = 1_600.0
        )
    )
    earth = LineCableModels.homogeneous(
        rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
    )
    connections = [Dict(:core => 3index - 2, :sheath => 3index - 1, :armor => 3index)
                   for index in 1:3]
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        fill(design, 3),
        [LineCableModels.Pose2(p.cable_x[index], p.cable_y) for index in 1:3];
        connections,
        environment = earth,
        system_id = "cable_380kv_armoured_ac_flat",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = earth,
        frequencies = p.frequencies
    )
end
