case_definition(
    :two_bare_wires,
    (
        core_radius = case_parameter(
            :core_radius, 0.0425; tags = (:geometry, :cable_layer)
        ),
        insulation_thickness = case_parameter(
            :insulation_thickness, 1.0e-3; tags = (:geometry, :cable_layer)
        ),
        first_x = case_parameter(:first_x, 0.0; tags = (:geometry, :system)),
        second_x = case_parameter(:second_x, 1.0; tags = (:geometry, :system)),
        cable_y = case_parameter(:cable_y, -1.0; tags = (:geometry, :system)),
        line_length = case_parameter(:line_length, 1.0; tags = (:operation, :system)),
        temperature = case_parameter(
            :temperature, 20.0; tags = (:operation, :system)
        ),
        earth_rho = case_parameter(:earth_rho, 0.1; tags = (:material, :earth)),
        frequencies = case_parameter(
            :frequencies, collect(10.0 .^ range(0, stop = 6, length = 101));
            tags = (:operation, :frequency)
        )
    ),
    ["cable:1:core", "cable:2:core"]
) do p
    materials = LineCableModels.MaterialsLibrary(add_defaults = true)
    copper = LineCableModels.Material(materials, :copper)
    artificial_insulation = LineCableModels.Material(
        kind = :insulator,
        rho = 1.97e14, eps_r = 2.3, mu_r = 1.0, T0 = 20.0, alpha = 0.0
    )
    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "two_bare_wires",
        LineCableModels.Stack(
            LineCableModels.Group(
                :core,
                LineCableModels.Region(
                    :core_metal,
                    LineCableModels.Disk(p.core_radius),
                    copper
                )
            ),
            LineCableModels.Region(
                :core_insulation,
                LineCableModels.Shell(p.insulation_thickness),
                artificial_insulation
            )
        )
    )
    earth = LineCableModels.Earth(rho = p.earth_rho, eps_r = 1.0, mu_r = 1.0)
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        [design, design],
        [
            LineCableModels.Pose2(p.first_x, p.cable_y),
            LineCableModels.Pose2(p.second_x, p.cable_y)
        ];
        connections = [Dict(:core => 1), Dict(:core => 2)],
        system_id = "two_bare_wires",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = earth,
        frequencies = p.frequencies
    )
end
