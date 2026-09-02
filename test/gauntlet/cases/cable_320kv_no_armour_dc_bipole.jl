case_definition(
    :cable_320kv_no_armour_dc_bipole,
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
        cable_x = case_parameter(
            :cable_x, (-0.5, 0.5); tags = (:geometry, :system)
        ),
        cable_y = case_parameter(:cable_y, -1.0; tags = (:geometry, :system)),
        line_length = case_parameter(
            :line_length, 150_000.0; tags = (:operation, :system)
        ),
        voltage = case_parameter(:voltage, 320.0; tags = (:operation, :system)),
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
        "cable:1:core", "cable:1:sheath",
        "cable:2:core", "cable:2:sheath"
    ]
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
    matrix = LineCableModels.Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )

    parts = LineCableModels.AbstractCablePart[]
    strand_radius = p.strand_diameter / 2
    radius = zero(strand_radius)
    for layer in 0:p.strand_layers
        outer = layer == 0 ? strand_radius : radius + 2strand_radius
        push!(parts,
            LineCableModels.Group(
                :core,
                LineCableModels.Region(
                    Symbol(:core_strands_, layer + 1),
                    LineCableModels.Disk(strand_radius),
                    core
                );
                pattern = layer == 0 ?
                          LineCableModels.Ring(1; r = zero(radius)) :
                          LineCableModels.Hexa(layer),
                path = layer == 0 ? nothing :
                       LineCableModels.Helix(LineCableModels.LayRatio(p.core_lay_ratio))
            ))
        radius = outer
    end
    core_enclosure = LineCableModels.Enclosure(
        :core_matrix,
        LineCableModels.Stack(parts);
        primitive = LineCableModels.Disk(radius),
        fill = matrix
    )
    parts = LineCableModels.AbstractCablePart[core_enclosure]
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
            :sheath_outer_jacket,
            LineCableModels.Shell(p.jacket_thickness),
            pp
        ))

    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "cable_320kv_no_armour",
        LineCableModels.Stack(parts);
        nominal_data = (
            designation_code = "cable_320kv_no_armour",
            U0 = p.voltage,
            conductor_cross_section = 1_600.0
        )
    )
    earth = LineCableModels.homogeneous(
        rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
    )
    connections = [Dict(:core => 2index - 1, :sheath => 2index) for index in 1:2]
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        fill(design, 2),
        [LineCableModels.Pose2(p.cable_x[index], p.cable_y) for index in 1:2];
        connections,
        environment = earth,
        system_id = "cable_320kv_no_armour_dc_bipole",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = earth,
        frequencies = p.frequencies
    )
end
