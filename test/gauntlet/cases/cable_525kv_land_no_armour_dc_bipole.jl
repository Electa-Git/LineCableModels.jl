case_definition(
    :cable_525kv_land_no_armour_dc_bipole,
    (
        core_diameter = case_parameter(
            :core_diameter, 68.0e-3; tags = (:geometry, :cable_layer)
        ),
        core_r20 = case_parameter(
            :core_r20, 0.006; tags = (:material, :cable_layer)
        ),
        inner_semicon_thickness = case_parameter(
            :inner_semicon_thickness, 1.8e-3; tags = (:geometry, :cable_layer)
        ),
        insulation_thickness = case_parameter(
            :insulation_thickness, 26.5e-3; tags = (:geometry, :cable_layer)
        ),
        outer_semicon_thickness = case_parameter(
            :outer_semicon_thickness, 1.5e-3; tags = (:geometry, :cable_layer)
        ),
        water_blocking_thickness = case_parameter(
            :water_blocking_thickness, 0.3e-3; tags = (:geometry, :cable_layer)
        ),
        pre_sheath_filler = case_parameter(
            :pre_sheath_filler, 3.7e-3; tags = (:geometry, :cable_layer)
        ),
        aluminum_sheath_thickness = case_parameter(
            :aluminum_sheath_thickness, 1.2e-3; tags = (:geometry, :cable_layer)
        ),
        post_sheath_filler = case_parameter(
            :post_sheath_filler, 2.2e-3; tags = (:geometry, :cable_layer)
        ),
        jacket_thickness = case_parameter(
            :jacket_thickness, 5.0e-3; tags = (:geometry, :cable_layer)
        ),
        outer_semicon_skin = case_parameter(
            :outer_semicon_skin, 0.3e-3; tags = (:geometry, :cable_layer)
        ),
        xlpe_rho = case_parameter(
            :xlpe_rho, 1.0e14; tags = (:material, :cable_layer)
        ),
        xlpe_eps_r = case_parameter(
            :xlpe_eps_r, 2.4; tags = (:material, :cable_layer)
        ),
        aluminum_rho = case_parameter(
            :aluminum_rho, 2.8264e-8; tags = (:material, :cable_layer)
        ),
        pe_rho = case_parameter(
            :pe_rho, 1.97e14; tags = (:material, :cable_layer)
        ),
        pe_eps_r = case_parameter(
            :pe_eps_r, 2.5; tags = (:material, :cable_layer)
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
        "cable:1:core", "cable:1:sheath",
        "cable:2:core", "cable:2:sheath"
    ]
) do p
    core_rho = p.core_r20 * (pi * (p.core_diameter / 2)^2) / 1_000
    core = LineCableModels.Material(
        :conductor, core_rho, 1.0, 1.0, 20.0, 0.00393
    )
    semicon_inner = LineCableModels.Material(:semicon, 1_000.0, 1_000.0)
    semicon_outer = LineCableModels.Material(:semicon, 500.0, 1_000.0)
    xlpe = LineCableModels.Material(:insulator, p.xlpe_rho, p.xlpe_eps_r)
    aluminum = LineCableModels.Material(
        :conductor, p.aluminum_rho, 1.0, 1.000022, 20.0, 0.00429
    )
    pe = LineCableModels.Material(:insulator, p.pe_rho, p.pe_eps_r)

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
        (:core_water_blocking, p.water_blocking_thickness, semicon_outer),
        (:core_pre_sheath_filler, p.pre_sheath_filler, pe)
    )
        push!(parts,
            LineCableModels.Region(tag, LineCableModels.Shell(layer_thickness), material))
        radius += layer_thickness
    end
    aluminum_outer = radius + p.aluminum_sheath_thickness
    push!(parts,
        LineCableModels.Group(
            :sheath,
            LineCableModels.Region(
                :sheath_aluminum,
                LineCableModels.Annulus(radius, aluminum_outer),
                aluminum
            )
        ))
    for (tag, layer_thickness, material) in (
        (:sheath_post_filler, p.post_sheath_filler, pe),
        (:sheath_outer_jacket, p.jacket_thickness, pe),
        (:sheath_outer_semicon, p.outer_semicon_skin, semicon_outer)
    )
        push!(parts,
            LineCableModels.Region(tag, LineCableModels.Shell(layer_thickness), material))
    end

    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "cable_525kv_land_no_armour",
        LineCableModels.Stack(parts);
        nominal_data = (
            designation_code = "cable_525kv_land_no_armour",
            U0 = p.voltage,
            conductor_cross_section = 3_000.0
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
        system_id = "cable_525kv_land_no_armour_dc_bipole",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = earth,
        frequencies = p.frequencies
    )
end
