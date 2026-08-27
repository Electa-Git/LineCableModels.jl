case_definition(
    :cable_525kv_1600mm2_bipole,
    (
        strand_layers = case_parameter(
            :strand_layers, 6; tags = (:topology, :cable_layer)
        ),
        strands_per_layer = case_parameter(
            :strands_per_layer, 6; tags = (:topology, :cable_layer)
        ),
        strand_diameter = case_parameter(
            :strand_diameter, 3.6649e-3; tags = (:geometry, :cable_layer)
        ),
        core_lay_ratio = case_parameter(
            :core_lay_ratio, 11.0; tags = (:geometry, :cable_layer)
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
            :water_blocking_thickness, 0.3e-3; tags = (:geometry, :cable_layer)
        ),
        lead_screen_thickness = case_parameter(
            :lead_screen_thickness, 3.3e-3; tags = (:geometry, :cable_layer)
        ),
        inner_sheath_thickness = case_parameter(
            :inner_sheath_thickness, 3.0e-3; tags = (:geometry, :cable_layer)
        ),
        bedding_thickness = case_parameter(
            :bedding_thickness, 3.0e-3; tags = (:geometry, :cable_layer)
        ),
        armor_wires = case_parameter(
            :armor_wires, 68; tags = (:topology, :cable_layer)
        ),
        armor_wire_diameter = case_parameter(
            :armor_wire_diameter, 5.827e-3; tags = (:geometry, :cable_layer)
        ),
        armor_lay_ratio = case_parameter(
            :armor_lay_ratio, 10.0; tags = (:geometry, :cable_layer)
        ),
        armor_packing_clearance_ratio = case_parameter(
            :armor_packing_clearance_ratio, 0.2;
            tags = (:constraint, :cable_layer)
        ),
        jacket_thickness = case_parameter(
            :jacket_thickness, 10.0e-3; tags = (:geometry, :cable_layer)
        ),
        cable_x = case_parameter(
            :cable_x, (-0.5, 0.5); tags = (:geometry, :system)
        ),
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
    [
        "cable:1:core", "cable:1:sheath", "cable:1:armor",
        "cable:2:core", "cable:2:sheath", "cable:2:armor"
    ]
) do p
    materials = LineCableModels.MaterialsLibrary(add_defaults = true)
    copper = LineCableModels.Material(materials, :copper)
    semicon1 = LineCableModels.Material(materials, :semicon1)
    semicon2 = LineCableModels.Material(materials, :semicon2)
    pe = LineCableModels.Material(materials, :pe)
    polyacrylate = LineCableModels.Material(materials, :polyacrylate)
    lead = LineCableModels.Material(materials, :lead)
    pp = LineCableModels.Material(materials, :pp)
    steel = LineCableModels.Material(materials, :steel)
    radius_before_bedding =
        (p.strand_layers + 0.5) * p.strand_diameter +
        p.inner_semicon_thickness +
        p.insulation_thickness +
        p.outer_semicon_thickness +
        p.water_blocking_thickness +
        p.lead_screen_thickness +
        p.inner_sheath_thickness
    armor_wire_radius = p.armor_wire_diameter / 2
    minimum_armor_radius =
        armor_wire_radius / sinpi(1 / p.armor_wires) - armor_wire_radius
    unbuffered_armor_radius = radius_before_bedding + p.bedding_thickness
    unbuffered_outer_radius =
        unbuffered_armor_radius + p.armor_wire_diameter + p.jacket_thickness
    packing_buffer =
        p.armor_packing_clearance_ratio * unbuffered_outer_radius
    buffered_armor_radius = unbuffered_armor_radius + packing_buffer
    # Fixed-count armor wires cannot occupy a ring smaller than their packing
    # radius. A fixed design clearance keeps ordinary 10% manufacturing draws
    # away from that boundary; the compliant bedding still absorbs any extreme
    # positive shortfall instead of constructing an unsupported wire packing.
    packing_shortfall = max(
        minimum_armor_radius - buffered_armor_radius,
        zero(minimum_armor_radius - buffered_armor_radius)
    )
    bedding_thickness = p.bedding_thickness + packing_buffer + packing_shortfall
    core_parts = (
        LineCableModels.Conductor.Stranded(
            :core;
            layers = p.strand_layers + 1,
            wire_radius = p.strand_diameter / 2,
            num_wires = p.strands_per_layer,
            lay_ratio = p.core_lay_ratio,
            material = copper
        ),
        LineCableModels.Insulator.Semicon(
            :core; thickness = p.inner_semicon_thickness, material = semicon1
        ),
        LineCableModels.Insulator.Tubular(
            :core; thickness = p.insulation_thickness, material = pe
        ),
        LineCableModels.Insulator.Semicon(
            :core; thickness = p.outer_semicon_thickness, material = semicon2
        ),
        LineCableModels.Insulator.Semicon(
            :core; thickness = p.water_blocking_thickness, material = polyacrylate
        )
    )
    sheath_parts = (
        LineCableModels.Conductor.Tubular(
            :sheath; thickness = p.lead_screen_thickness, material = lead
        ),
        LineCableModels.Insulator.Tubular(
            :sheath; thickness = p.inner_sheath_thickness, material = pe
        ),
        LineCableModels.Insulator.Tubular(
            :sheath; thickness = bedding_thickness, material = pp
        )
    )
    armor_parts = (
        LineCableModels.Conductor.Wires(
            :armor;
            wire_radius = armor_wire_radius,
            num_wires = p.armor_wires,
            lay_ratio = p.armor_lay_ratio,
            material = steel
        ),
        LineCableModels.Insulator.Tubular(
            :armor; thickness = p.jacket_thickness, material = pp
        )
    )
    nominal_data = (
        designation_code = "(N)2XH(F)RK2Y",
        U0 = 500.0,
        U = 525.0,
        conductor_cross_section = 1600.0,
        screen_cross_section = 1000.0,
        resistance = nothing,
        capacitance = nothing,
        inductance = nothing
    )
    design = LineCableModels.CableBuilder(
        "525kV_1600mm2", core_parts, sheath_parts, armor_parts;
        nominal = nominal_data
    )
    earth = LineCableModels.Earth(
        rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
    )
    positions = Tuple(
        LineCableModels.at(
            x = p.cable_x[index],
            y = p.cable_y,
            phases = (
                :core => 3index - 2,
                :sheath => 3index - 1,
                :armor => 3index
            )
        ) for index in 1:2
    )
    return LineCableModels.SystemBuilder(
        "cable_525kv_1600mm2_bipole",
        design,
        positions;
        length = p.line_length,
        temperature = p.temperature,
        earth,
        frequencies = p.frequencies
    )
end
