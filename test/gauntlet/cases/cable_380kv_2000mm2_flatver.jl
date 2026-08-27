case_definition(
    :cable_380kv_2000mm2_flatver,
    (
        core_layers = case_parameter(:core_layers, 5; tags = (:topology, :cable_layer)),
        core_strand_diameter = case_parameter(
            :core_strand_diameter, 6.543707496202785e-3;
            tags = (:geometry, :cable_layer)
        ),
        core_lay_ratio = case_parameter(
            :core_lay_ratio, 13.0; tags = (:geometry, :cable_layer)
        ),
        screen_wires = case_parameter(
            :screen_wires, 14; tags = (:topology, :cable_layer)
        ),
        screen_wire_diameter = case_parameter(
            :screen_wire_diameter, 6.543707496202785e-3;
            tags = (:geometry, :cable_layer)
        ),
        screen_lay_ratio = case_parameter(
            :screen_lay_ratio, 10.0; tags = (:geometry, :cable_layer)
        ),
        semicon_tape_thickness = case_parameter(
            :semicon_tape_thickness, 0.3e-3; tags = (:geometry, :cable_layer)
        ),
        inner_semicon_thickness = case_parameter(
            :inner_semicon_thickness, 0.768e-3; tags = (:geometry, :cable_layer)
        ),
        insulation_thickness = case_parameter(
            :insulation_thickness, 28.0e-3; tags = (:geometry, :cable_layer)
        ),
        outer_semicon_thickness = case_parameter(
            :outer_semicon_thickness, 0.472e-3; tags = (:geometry, :cable_layer)
        ),
        copper_tape_thickness = case_parameter(
            :copper_tape_thickness, 0.1e-3; tags = (:geometry, :cable_layer)
        ),
        copper_tape_width = case_parameter(
            :copper_tape_width, 10.0e-3; tags = (:geometry, :cable_layer)
        ),
        copper_tape_lay_ratio = case_parameter(
            :copper_tape_lay_ratio, 10.0; tags = (:geometry, :cable_layer)
        ),
        water_blocking_thickness = case_parameter(
            :water_blocking_thickness, 0.94e-3; tags = (:geometry, :cable_layer)
        ),
        aluminum_tape_thickness = case_parameter(
            :aluminum_tape_thickness, 0.15e-3; tags = (:geometry, :cable_layer)
        ),
        jacket_thickness = case_parameter(
            :jacket_thickness, 5.05e-3; tags = (:geometry, :cable_layer)
        ),
        cable_x = case_parameter(:cable_x, 0.0; tags = (:geometry, :system)),
        cable_y = case_parameter(
            :cable_y, (-1.0, -1.5, -2.0); tags = (:geometry, :system)
        ),
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
        "cable:1:core", "cable:1:sheath", "cable:1:jacket",
        "cable:2:core", "cable:2:sheath", "cable:2:jacket",
        "cable:3:core", "cable:3:sheath", "cable:3:jacket"
    ]
) do p
    materials = LineCableModels.MaterialsLibrary(add_defaults = true)
    aluminum = LineCableModels.Material(materials, :aluminum)
    copper = LineCableModels.Material(materials, :copper)
    pe = LineCableModels.Material(materials, :pe)
    xlpe = LineCableModels.Material(materials, :xlpe)
    semicon1 = LineCableModels.Material(materials, :semicon1)
    semicon2 = LineCableModels.Material(materials, :semicon2)
    polyacrylate = LineCableModels.Material(materials, :polyacrylate)
    core_parts = (
        LineCableModels.Conductor.Stranded(
            :core;
            layers = p.core_layers,
            wire_radius = p.core_strand_diameter / 2,
            num_wires = 6,
            lay_ratio = p.core_lay_ratio,
            material = copper
        ),
        LineCableModels.Insulator.Semicon(
            :core; thickness = p.semicon_tape_thickness, material = polyacrylate
        ),
        LineCableModels.Insulator.Semicon(
            :core; thickness = p.inner_semicon_thickness, material = semicon1
        ),
        LineCableModels.Insulator.Tubular(
            :core; thickness = p.insulation_thickness, material = xlpe
        ),
        LineCableModels.Insulator.Semicon(
            :core; thickness = p.outer_semicon_thickness, material = semicon2
        ),
        LineCableModels.Insulator.Semicon(
            :core; thickness = p.semicon_tape_thickness, material = polyacrylate
        )
    )
    sheath_parts = (
        LineCableModels.Conductor.Wires(
            :sheath;
            wire_radius = p.screen_wire_diameter / 2,
            num_wires = p.screen_wires,
            lay_ratio = p.screen_lay_ratio,
            material = copper
        ),
        LineCableModels.Conductor.Strip(
            :sheath;
            thickness = p.copper_tape_thickness,
            width = p.copper_tape_width,
            lay_ratio = p.copper_tape_lay_ratio,
            material = copper
        ),
        LineCableModels.Insulator.Semicon(
            :sheath; thickness = p.water_blocking_thickness, material = polyacrylate
        )
    )
    jacket_parts = (
        LineCableModels.Conductor.Tubular(
            :jacket; thickness = p.aluminum_tape_thickness, material = aluminum
        ),
        LineCableModels.Insulator.Tubular(
            :jacket; thickness = p.jacket_thickness, material = pe
        )
    )
    design = LineCableModels.CableBuilder(
        "380kV_2000mm2", core_parts, sheath_parts, jacket_parts;
        nominal = LineCableModels.DataModel.NominalData()
    )
    earth = LineCableModels.Earth(
        rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
    )
    positions = Tuple(
        LineCableModels.at(
            x = p.cable_x,
            y = p.cable_y[index],
            phases = (
                :core => 3index - 2,
                :sheath => 3index - 1,
                :jacket => 3index
            )
        ) for index in 1:3
    )
    return LineCableModels.SystemBuilder(
        "cable_380kv_2000mm2_flatver",
        design,
        positions;
        length = p.line_length,
        temperature = p.temperature,
        earth,
        frequencies = p.frequencies
    )
end
