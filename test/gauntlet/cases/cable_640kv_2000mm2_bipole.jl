case_definition(
    :cable_640kv_2000mm2_bipole,
    (
        core_strand_radius = case_parameter(
            :core_strand_radius, 6.543707496202785e-3 / 2;
            tags = (:geometry, :cable_layer)
        ),
        ring_counts = case_parameter(
            :ring_counts, (6, 12, 18, 24); tags = (:topology, :cable_layer)
        ),
        core_lay_ratio = case_parameter(
            :core_lay_ratio, 13.0; tags = (:geometry, :cable_layer)
        ),
        semicon_tape_thickness = case_parameter(
            :semicon_tape_thickness, 0.3e-3; tags = (:geometry, :cable_layer)
        ),
        inner_semicon_thickness = case_parameter(
            :inner_semicon_thickness, 0.768e-3; tags = (:geometry, :cable_layer)
        ),
        insulation_thickness = case_parameter(
            :insulation_thickness, 32.0e-3; tags = (:geometry, :cable_layer)
        ),
        outer_semicon_thickness = case_parameter(
            :outer_semicon_thickness, 0.472e-3; tags = (:geometry, :cable_layer)
        ),
        lead_screen_thickness = case_parameter(
            :lead_screen_thickness, 3.3e-3; tags = (:geometry, :cable_layer)
        ),
        inner_sheath_thickness = case_parameter(
            :inner_sheath_thickness, 3.0e-3; tags = (:geometry, :cable_layer)
        ),
        aluminum_tape_thickness = case_parameter(
            :aluminum_tape_thickness, 0.15e-3; tags = (:geometry, :cable_layer)
        ),
        jacket_thickness = case_parameter(
            :jacket_thickness, 6.05e-3; tags = (:geometry, :cable_layer)
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
        "cable:1:core", "cable:1:sheath", "cable:1:jacket",
        "cable:2:core", "cable:2:sheath", "cable:2:jacket"
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
    lead = LineCableModels.Material(materials, :lead)
    DM = LineCableModels.DataModel
    core_strand_area = π * p.core_strand_radius^2
    core_conductor = let group = DM.ConductorGroup(DM.Tubular(
            0.0, p.core_strand_radius, copper
        ))
        for count in p.ring_counts
            radius_in = group.r_ex
            strand_width = prevfloat(2π * radius_in / count)
            strand_thickness = core_strand_area / strand_width
            radius_ext = sqrt(
                radius_in^2 + count * strand_width * strand_thickness / π
            )
            LineCableModels.add!(
                group,
                DM.RectStrands(
                    radius_in,
                    radius_ext,
                    strand_thickness,
                    strand_width,
                    count,
                    p.core_lay_ratio,
                    copper
                )
            )
        end
        group
    end
    radius_in = core_conductor.r_ex
    core_insulation = DM.InsulatorGroup(DM.Semicon(
        radius_in, radius_in + p.semicon_tape_thickness, polyacrylate
    ))
    radius_in = core_insulation.r_ex
    LineCableModels.add!(core_insulation, DM.Semicon(
        radius_in, radius_in + p.inner_semicon_thickness, semicon1
    ))
    radius_in = core_insulation.r_ex
    LineCableModels.add!(core_insulation, DM.Insulator(
        radius_in, radius_in + p.insulation_thickness, xlpe
    ))
    radius_in = core_insulation.r_ex
    LineCableModels.add!(core_insulation, DM.Semicon(
        radius_in, radius_in + p.outer_semicon_thickness, semicon2
    ))
    radius_in = core_insulation.r_ex
    LineCableModels.add!(core_insulation, DM.Semicon(
        radius_in, radius_in + p.semicon_tape_thickness, polyacrylate
    ))
    core_component = DM.CableComponent("core", core_conductor, core_insulation)

    radius_in = core_insulation.r_ex
    sheath_conductor = DM.ConductorGroup(DM.Tubular(
        radius_in, radius_in + p.lead_screen_thickness, lead
    ))
    radius_in = sheath_conductor.r_ex
    sheath_insulation = DM.InsulatorGroup(DM.Insulator(
        radius_in, radius_in + p.inner_sheath_thickness, pe
    ))
    sheath_component = DM.CableComponent(
        "sheath", sheath_conductor, sheath_insulation
    )

    radius_in = sheath_insulation.r_ex
    jacket_conductor = DM.ConductorGroup(DM.Tubular(
        radius_in, radius_in + p.aluminum_tape_thickness, aluminum
    ))
    radius_in = jacket_conductor.r_ex
    jacket_insulation = DM.InsulatorGroup(DM.Insulator(
        radius_in, radius_in + p.jacket_thickness, pe
    ))
    jacket_component = DM.CableComponent(
        "jacket", jacket_conductor, jacket_insulation
    )
    design = DM.CableDesign(
        "640kV_2000mm2",
        [core_component, sheath_component, jacket_component];
        nominal_data = DM.NominalData()
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
                :jacket => 3index
            )
        ) for index in 1:2
    )
    return LineCableModels.SystemBuilder(
        "cable_640kv_2000mm2_bipole",
        design,
        positions;
        length = p.line_length,
        temperature = p.temperature,
        earth,
        frequencies = p.frequencies
    )
end
