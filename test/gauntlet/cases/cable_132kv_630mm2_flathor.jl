case_definition(
    :cable_132kv_630mm2_flathor,
    (
        core_layers = case_parameter(:core_layers, 5; tags = (:topology, :cable_layer)),
        core_strand_diameter = case_parameter(
            :core_strand_diameter, 3.6648500338805653e-3;
            tags = (:geometry, :cable_layer)
        ),
        core_lay_ratio = case_parameter(
            :core_lay_ratio, 13.0; tags = (:geometry, :cable_layer)
        ),
        screen_wires = case_parameter(
            :screen_wires, 19; tags = (:topology, :cable_layer)
        ),
        screen_wire_diameter = case_parameter(
            :screen_wire_diameter, 2.5881867280128637e-3;
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
            :insulation_thickness, 16.0e-3; tags = (:geometry, :cable_layer)
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
        cable_x = case_parameter(
            :cable_x, (-0.5, 0.0, 0.5); tags = (:geometry, :system)
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
    parts = LineCableModels.AbstractCablePart[]
    wire_radius = p.core_strand_diameter / 2
    radius = zero(wire_radius)
    for layer in 0:(p.core_layers - 1)
        count = layer == 0 ? 1 : 6layer
        centre_radius = count == 1 ? zero(radius) : radius + wire_radius
        push!(parts,
            LineCableModels.Group(
                :core,
                LineCableModels.Region(
                    Symbol(:core_strands_, layer + 1),
                    LineCableModels.Disk(wire_radius),
                    copper
                );
                pattern = LineCableModels.Ring(count; r = centre_radius),
                path = count == 1 ? nothing :
                       LineCableModels.Helix(LineCableModels.LayRatio(p.core_lay_ratio))
            ))
        radius = count == 1 ? wire_radius : radius + 2wire_radius
    end
    for (tag, thickness, material) in (
        (:core_semicon_tape_inner, p.semicon_tape_thickness, polyacrylate),
        (:core_semicon_inner, p.inner_semicon_thickness, semicon1),
        (:core_insulation, p.insulation_thickness, xlpe),
        (:core_semicon_outer, p.outer_semicon_thickness, semicon2),
        (:core_semicon_tape_outer, p.semicon_tape_thickness, polyacrylate)
    )
        push!(parts, LineCableModels.Region(
            tag, LineCableModels.Shell(thickness), material
        ))
        radius += thickness
    end

    screen_radius = p.screen_wire_diameter / 2
    screen_centre = radius + screen_radius
    push!(parts,
        LineCableModels.Group(
            :sheath,
            LineCableModels.Region(
                :sheath_wires, LineCableModels.Disk(screen_radius), copper
            );
            pattern = LineCableModels.Ring(p.screen_wires; r = screen_centre),
            path = LineCableModels.Helix(LineCableModels.LayRatio(p.screen_lay_ratio))
        ))
    radius += 2screen_radius
    tape_outer = radius + p.copper_tape_thickness
    tape_span = p.copper_tape_width / ((radius + tape_outer) / 2)
    push!(parts,
        LineCableModels.Group(
            :sheath,
            LineCableModels.Region(
                :sheath_copper_tape,
                LineCableModels.Sector(radius, tape_outer, -tape_span / 2, tape_span),
                copper
            );
            path = LineCableModels.Helix(LineCableModels.LayRatio(p.copper_tape_lay_ratio))
        ))
    radius = tape_outer
    push!(parts,
        LineCableModels.Region(
            :sheath_water_blocking,
            LineCableModels.Shell(p.water_blocking_thickness),
            polyacrylate
        ))
    radius += p.water_blocking_thickness

    aluminum_outer = radius + p.aluminum_tape_thickness
    push!(parts,
        LineCableModels.Group(
            :jacket,
            LineCableModels.Region(
                :jacket_aluminum_tape,
                LineCableModels.Annulus(radius, aluminum_outer),
                aluminum
            )
        ))
    radius = aluminum_outer
    push!(parts, LineCableModels.Region(
        :jacket_insulation,
        LineCableModels.Shell(p.jacket_thickness),
        pe
    ))
    nominal_data = (
        designation_code = "A2XS(FL)2Y",
        U0 = 76.0,
        U = 132.0,
        conductor_cross_section = 630.0,
        screen_cross_section = 95.0,
        resistance = 0.0395,
        capacitance = 0.2,
        inductance = 0.12
    )
    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "132kV_630mm2",
        LineCableModels.Stack(parts);
        nominal_data = LineCableModels.NominalData(; nominal_data...)
    )
    earth = LineCableModels.Earth(
        rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
    )
    designs = fill(design, 3)
    positions = [LineCableModels.Pose2(p.cable_x[index], p.cable_y) for index in 1:3]
    connections = [Dict(:core => 3index - 2, :sheath => 3index - 1, :jacket => 3index)
                   for index in 1:3]
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        designs,
        positions;
        connections,
        system_id = "cable_132kv_630mm2_flathor",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = earth,
        frequencies = p.frequencies
    )
end
