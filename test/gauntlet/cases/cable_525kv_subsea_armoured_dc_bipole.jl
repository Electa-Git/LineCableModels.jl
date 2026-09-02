case_definition(
    :cable_525kv_subsea_armoured_dc_bipole,
    (
        core_diameter = case_parameter(
            :core_diameter, 60.0e-3; tags = (:geometry, :cable_layer)
        ),
        core_r20 = case_parameter(
            :core_r20, 0.0072; tags = (:material, :cable_layer)
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
            :water_blocking_thickness, 0.6e-3; tags = (:geometry, :cable_layer)
        ),
        lead_sheath_thickness = case_parameter(
            :lead_sheath_thickness, 3.2e-3; tags = (:geometry, :cable_layer)
        ),
        inner_pe_thickness = case_parameter(
            :inner_pe_thickness, 3.3e-3; tags = (:geometry, :cable_layer)
        ),
        bedding_thickness = case_parameter(
            :bedding_thickness, 1.0e-3; tags = (:geometry, :cable_layer)
        ),
        armor_wires = case_parameter(
            :armor_wires, 68; tags = (:topology, :cable_layer)
        ),
        armor_wire_diameter = case_parameter(
            :armor_wire_diameter, 6.0e-3; tags = (:geometry, :cable_layer)
        ),
        armor_lay_ratio = case_parameter(
            :armor_lay_ratio, pi / sqrt(1.048^2 - 1);
            tags = (:geometry, :cable_layer)
        ),
        jacket_thickness = case_parameter(
            :jacket_thickness, 6.6e-3; tags = (:geometry, :cable_layer)
        ),
        xlpe_rho = case_parameter(
            :xlpe_rho, 1.0e14; tags = (:material, :cable_layer)
        ),
        xlpe_eps_r = case_parameter(
            :xlpe_eps_r, 2.4; tags = (:material, :cable_layer)
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
            :pe_eps_r, 2.5; tags = (:material, :cable_layer)
        ),
        pp_rho = case_parameter(
            :pp_rho, 1.0e15; tags = (:material, :cable_layer)
        ),
        pp_eps_r = case_parameter(
            :pp_eps_r, 2.8; tags = (:material, :cable_layer)
        ),
        steel_rho = case_parameter(
            :steel_rho, 1.38e-7; tags = (:material, :cable_layer)
        ),
        steel_mu_r = case_parameter(
            :steel_mu_r, 10.0; tags = (:material, :cable_layer)
        ),
        cable_x = case_parameter(
            :cable_x, (-0.5, 0.5); tags = (:geometry, :system)
        ),
        cable_y = case_parameter(:cable_y, -1.0; tags = (:geometry, :system)),
        line_length = case_parameter(
            :line_length, 150_000.0; tags = (:operation, :system)
        ),
        voltage = case_parameter(:voltage, 525.0; tags = (:operation, :system)),
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
        "cable:2:core", "cable:2:sheath", "cable:2:armor"
    ]
) do p
    core_rho = p.core_r20 * (pi * (p.core_diameter / 2)^2) / 1_000
    core = LineCableModels.Material(
        :conductor, core_rho, 1.0, 1.0, 20.0, 0.00393
    )
    semicon_inner = LineCableModels.Material(:semicon, 1_000.0, 1_000.0)
    semicon_outer = LineCableModels.Material(:semicon, 500.0, 1_000.0)
    xlpe = LineCableModels.Material(:insulator, p.xlpe_rho, p.xlpe_eps_r)
    lead = LineCableModels.Material(
        :conductor, p.lead_rho, 1.0, p.lead_mu_r, 20.0, 0.004
    )
    pe = LineCableModels.Material(:insulator, p.pe_rho, p.pe_eps_r)
    pp = LineCableModels.Material(:insulator, p.pp_rho, p.pp_eps_r)
    steel = LineCableModels.Material(
        :conductor, p.steel_rho, 1.0, p.steel_mu_r, 20.0, 0.0045
    )
    matrix = LineCableModels.Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )

    radius = p.core_diameter / 2
    parts = LineCableModels.AbstractCablePart[
        LineCableModels.Group(
        :core,
        LineCableModels.Region(:core_conductor, LineCableModels.Disk(radius), core)
    ),
    ]
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
    push!(parts,
        LineCableModels.Region(
            :sheath_inner_pe,
            LineCableModels.Shell(p.inner_pe_thickness),
            pe
        ))
    push!(parts,
        LineCableModels.Region(
            :sheath_bedding,
            LineCableModels.Shell(p.bedding_thickness),
            pp
        ))
    armor_wire_radius = p.armor_wire_diameter / 2
    armor_inner_radius = lead_outer + p.inner_pe_thickness + p.bedding_thickness
    armor_outer = armor_inner_radius + 2armor_wire_radius
    armor = LineCableModels.Group(
        :armor,
        LineCableModels.Region(
            :armor_wires,
            LineCableModels.Disk(armor_wire_radius),
            steel
        );
        pattern = LineCableModels.Ring(
            p.armor_wires; r = armor_inner_radius + armor_wire_radius
        ),
        path = LineCableModels.Helix(LineCableModels.LayRatio(p.armor_lay_ratio))
    )
    push!(parts,
        LineCableModels.Enclosure(
            :armor_matrix,
            armor;
            primitive = LineCableModels.Annulus(armor_inner_radius, armor_outer),
            fill = matrix
        ))
    push!(parts,
        LineCableModels.Region(
            :armor_jacket,
            LineCableModels.Shell(p.jacket_thickness),
            pp
        ))

    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "cable_525kv_subsea_armoured",
        LineCableModels.Stack(parts);
        nominal_data = (
            designation_code = "cable_525kv_subsea_armoured",
            U0 = p.voltage,
            conductor_cross_section = 2_500.0
        )
    )
    earth = LineCableModels.Earth(
        rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
    )
    connections = [Dict(:core => 3index - 2, :sheath => 3index - 1, :armor => 3index)
                   for index in 1:2]
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        fill(design, 2),
        [LineCableModels.Pose2(p.cable_x[index], p.cable_y) for index in 1:2];
        connections,
        environment = earth,
        system_id = "cable_525kv_subsea_armoured_dc_bipole",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = earth,
        frequencies = p.frequencies
    )
end
