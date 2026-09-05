case_definition(
    :solid_1000mm2_single,
    (
        core_cross_section = case_parameter(
            :core_cross_section, 1000.0e-6; tags = (:geometry, :cable_layer)
        ),
        insulation_thickness = case_parameter(
            :insulation_thickness, 8.3e-3; tags = (:geometry, :cable_layer)
        ),
        cable_x = case_parameter(:cable_x, 0.0; tags = (:geometry, :system)),
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
    ["cable:1:core"];
    description = "Single 1000 mm² solid conductor"
) do p
    materials = LineCableModels.MaterialsLibrary(add_defaults = true)
    aluminum = LineCableModels.Material(materials, :aluminum)
    xlpe = LineCableModels.Material(materials, :xlpe)
    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "solid_1000mm2",
        LineCableModels.Stack(
            LineCableModels.Group(
                :core,
                LineCableModels.Region(
                    :core_metal,
                    LineCableModels.Disk(sqrt(p.core_cross_section / π)),
                    aluminum
                )
            ),
            LineCableModels.Region(
                :core_insulation,
                LineCableModels.Shell(p.insulation_thickness),
                xlpe
            )
        )
    )
    earth = LineCableModels.homogeneous(
        rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
    )
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        design,
        LineCableModels.Pose2(p.cable_x, p.cable_y);
        connections = Dict(:core => 1),
        system_id = "solid_1000mm2_single",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = earth,
        frequencies = p.frequencies
    )
end
