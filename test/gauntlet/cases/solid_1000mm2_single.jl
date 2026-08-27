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
    ["cable:1:core"]
) do p
    materials = LineCableModels.MaterialsLibrary(add_defaults = true)
    aluminum = LineCableModels.Material(materials, :aluminum)
    xlpe = LineCableModels.Material(materials, :xlpe)
    design = LineCableModels.CableBuilder(
        "solid_1000mm2",
        LineCableModels.Conductor.Solid(
            :core; radius = sqrt(p.core_cross_section / π), material = aluminum
        ),
        LineCableModels.Insulator.Tubular(
            :core; thickness = p.insulation_thickness, material = xlpe
        );
        nominal = LineCableModels.DataModel.NominalData()
    )
    earth = LineCableModels.Earth(
        rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
    )
    positions = (
        LineCableModels.at(
            x = p.cable_x, y = p.cable_y, phases = (:core => 1,)
        ),
    )
    return LineCableModels.SystemBuilder(
        "solid_1000mm2_single",
        design,
        positions;
        length = p.line_length,
        temperature = p.temperature,
        earth,
        frequencies = p.frequencies
    )
end
