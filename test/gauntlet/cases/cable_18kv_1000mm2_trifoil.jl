case_definition(
    :cable_18kv_1000mm2_trifoil,
    (
        core_strand_diameter = case_parameter(
            :core_strand_diameter, 4.19e-3; tags = (:geometry, :cable_layer)
        ),
        core_ring_lay_ratio_1 = case_parameter(
            :core_ring_lay_ratio_1, 15.0; tags = (:geometry, :cable_layer)
        ),
        core_ring_lay_ratio_2 = case_parameter(
            :core_ring_lay_ratio_2, 13.5; tags = (:geometry, :cable_layer)
        ),
        core_ring_lay_ratio_3 = case_parameter(
            :core_ring_lay_ratio_3, 12.5; tags = (:geometry, :cable_layer)
        ),
        core_ring_lay_ratio_4 = case_parameter(
            :core_ring_lay_ratio_4, 11.0; tags = (:geometry, :cable_layer)
        ),
        core_ring_wire_counts = case_parameter(
            :core_ring_wire_counts, (1, 6, 12, 18, 24);
            tags = (:topology, :cable_layer)
        ),
        screen_wires = case_parameter(
            :screen_wires, 49; tags = (:topology, :cable_layer)
        ),
        screen_wire_diameter = case_parameter(
            :screen_wire_diameter, 0.94e-3; tags = (:geometry, :cable_layer)
        ),
        screen_wire_lay_ratio = case_parameter(
            :screen_wire_lay_ratio, 10.0; tags = (:geometry, :cable_layer)
        ),
        semicon_tape_thickness = case_parameter(
            :semicon_tape_thickness, 0.3e-3; tags = (:geometry, :cable_layer)
        ),
        inner_semicon_thickness = case_parameter(
            :inner_semicon_thickness, 0.768e-3; tags = (:geometry, :cable_layer)
        ),
        insulation_thickness = case_parameter(
            :insulation_thickness, 8.3e-3; tags = (:geometry, :cable_layer)
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
        pe_face_thickness = case_parameter(
            :pe_face_thickness, 0.05e-3; tags = (:geometry, :cable_layer)
        ),
        jacket_thickness = case_parameter(
            :jacket_thickness, 3.4e-3; tags = (:geometry, :cable_layer)
        ),
        formation_x = case_parameter(
            :formation_x, 0.0; tags = (:geometry, :system)
        ),
        formation_y = case_parameter(
            :formation_y, -1.0; tags = (:geometry, :system)
        ),
        formation_clearance_ratio = case_parameter(
            :formation_clearance_ratio, 0.2; tags = (:geometry, :system)
        ),
        line_length = case_parameter(
            :line_length, 1000.0; tags = (:operation, :system)
        ),
        temperature = case_parameter(
            :temperature, 20.0; tags = (:operation, :system)
        ),
        earth_rho = case_parameter(
            :earth_rho, 100.0; tags = (:material, :earth)
        ),
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

    counts = p.core_ring_wire_counts
    core_parts = (
        LineCableModels.Conductor.Wires(
            :core;
            wire_radius = p.core_strand_diameter / 2,
            num_wires = counts[1],
            lay_ratio = 0.0,
            material = aluminum
        ),
        LineCableModels.Conductor.Wires(
            :core;
            wire_radius = p.core_strand_diameter / 2,
            num_wires = counts[2],
            lay_ratio = p.core_ring_lay_ratio_1,
            material = aluminum
        ),
        LineCableModels.Conductor.Wires(
            :core;
            wire_radius = p.core_strand_diameter / 2,
            num_wires = counts[3],
            lay_ratio = p.core_ring_lay_ratio_2,
            material = aluminum
        ),
        LineCableModels.Conductor.Wires(
            :core;
            wire_radius = p.core_strand_diameter / 2,
            num_wires = counts[4],
            lay_ratio = p.core_ring_lay_ratio_3,
            material = aluminum
        ),
        LineCableModels.Conductor.Wires(
            :core;
            wire_radius = p.core_strand_diameter / 2,
            num_wires = counts[5],
            lay_ratio = p.core_ring_lay_ratio_4,
            material = aluminum
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
            lay_ratio = p.screen_wire_lay_ratio,
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
            :jacket; thickness = p.pe_face_thickness, material = pe
        ),
        LineCableModels.Insulator.Tubular(
            :jacket; thickness = p.jacket_thickness, material = pe
        )
    )
    nominal_data = (
        designation_code = "NA2XS(FL)2Y",
        U0 = 18.0,
        U = 30.0,
        conductor_cross_section = 1000.0,
        screen_cross_section = 35.0,
        resistance = 0.0291,
        capacitance = 0.39,
        inductance = 0.3
    )
    design = LineCableModels.CableBuilder(
        "18kV_1000mm2", core_parts, sheath_parts, jacket_parts;
        nominal = nominal_data
    )
    outer_radius = last(design.components).insulator_group.r_ex
    formation_spacing = (2 + p.formation_clearance_ratio) * outer_radius
    earth = LineCableModels.Earth(
        rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
    )
    positions = (
        LineCableModels.trifoil(
            x = p.formation_x,
            y = p.formation_y,
            spacing = formation_spacing,
            phases = (
                :core => (1, 4, 7),
                :sheath => (2, 5, 8),
                :jacket => (3, 6, 9)
            )
        ),
    )
    return LineCableModels.SystemBuilder(
        "cable_18kv_1000mm2_trifoil",
        design,
        positions;
        length = p.line_length,
        temperature = p.temperature,
        earth,
        frequencies = p.frequencies
    )
end
