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
    artificial_conductor = LineCableModels.Material(
        kind = :conductor,
        rho = eps(Float64), eps_r = 1.0, mu_r = 1.0, T0 = 20.0, alpha = 0.0
    )
    artificial_insulation = LineCableModels.Material(
        kind = :insulator,
        rho = 1.97e14, eps_r = 2.3, mu_r = 1.0, T0 = 20.0, alpha = 0.0
    )
    design = LineCableModels.CableBuilder(
        "two_bare_wires",
        LineCableModels.Conductor.Solid(
            :core; radius = p.core_radius, material = artificial_conductor
        ),
        LineCableModels.Insulator.Tubular(
            :core;
            thickness = p.insulation_thickness,
            material = artificial_insulation
        );
        nominal = LineCableModels.DataModel.NominalData()
    )
    earth = LineCableModels.Earth(rho = p.earth_rho, eps_r = 1.0, mu_r = 1.0)
    positions = (
        LineCableModels.at(
            x = p.first_x, y = p.cable_y, phases = (:core => 1,)
        ),
        LineCableModels.at(
            x = p.second_x, y = p.cable_y, phases = (:core => 2,)
        )
    )
    return LineCableModels.SystemBuilder(
        "two_bare_wires",
        design,
        positions;
        length = p.line_length,
        temperature = p.temperature,
        earth,
        frequencies = p.frequencies
    )
end
