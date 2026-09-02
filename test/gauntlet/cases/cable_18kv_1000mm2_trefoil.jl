case_definition(
    :cable_18kv_1000mm2_trefoil,
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
    matrix = LineCableModels.Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )

    counts = p.core_ring_wire_counts
    lay_ratios = (
        0.0,
        p.core_ring_lay_ratio_1,
        p.core_ring_lay_ratio_2,
        p.core_ring_lay_ratio_3,
        p.core_ring_lay_ratio_4
    )
    parts = LineCableModels.AbstractCablePart[]
    wire_radius = p.core_strand_diameter / 2
    radius = zero(wire_radius)
    for (layer, (count, lay_ratio)) in enumerate(zip(counts, lay_ratios))
        centre_radius = count == 1 ? zero(radius) : radius + wire_radius
        push!(parts,
            LineCableModels.Group(
                :core,
                LineCableModels.Region(
                    Symbol(:core_strands_, layer),
                    LineCableModels.Disk(wire_radius),
                    aluminum
                );
                pattern = LineCableModels.Ring(count; r = centre_radius),
                path = iszero(lay_ratio) ? nothing :
                       LineCableModels.Helix(LineCableModels.LayRatio(lay_ratio))
            ))
        radius = count == 1 ? wire_radius : radius + 2wire_radius
    end
    core = LineCableModels.Enclosure(
        :core_matrix,
        LineCableModels.Stack(parts);
        primitive = LineCableModels.Disk(radius),
        fill = matrix
    )
    parts = LineCableModels.AbstractCablePart[core]
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
    screen_inner = radius
    screen_centre = radius + screen_radius
    screen = LineCableModels.Group(
        :sheath,
        LineCableModels.Region(
            :sheath_wires, LineCableModels.Disk(screen_radius), copper
        );
        pattern = LineCableModels.Ring(p.screen_wires; r = screen_centre),
        path = LineCableModels.Helix(LineCableModels.LayRatio(p.screen_wire_lay_ratio))
    )
    radius += 2screen_radius
    tape_outer = radius + p.copper_tape_thickness
    tape_span = p.copper_tape_width / ((radius + tape_outer) / 2)
    tape = LineCableModels.Group(
        :sheath,
        LineCableModels.Region(
            :sheath_copper_tape,
            LineCableModels.Sector(radius, tape_outer, -tape_span / 2, tape_span),
            copper
        );
        path = LineCableModels.Helix(LineCableModels.LayRatio(p.copper_tape_lay_ratio))
    )
    push!(parts,
        LineCableModels.Enclosure(
            :screen_matrix,
            LineCableModels.Stack(screen, tape);
            primitive = LineCableModels.Annulus(screen_inner, tape_outer),
            fill = matrix
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
    push!(parts,
        LineCableModels.Region(
            :jacket_pe_face, LineCableModels.Shell(p.pe_face_thickness), pe
        ))
    push!(parts,
        LineCableModels.Region(
            :jacket_insulation, LineCableModels.Shell(p.jacket_thickness), pe
        ))
    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "18kV_1000mm2",
        LineCableModels.Stack(parts)
    )
    formation_spacing = (2 + p.formation_clearance_ratio) *
                        LineCableModels.outer_radius(design)
    earth = LineCableModels.homogeneous(
        rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
    )
    placements = LineCableModels.trefoil(
        design;
        center = LineCableModels.at(p.formation_x, p.formation_y),
        spacing = formation_spacing,
        connections = (
            core = (1, 4, 7),
            sheath = (2, 5, 8),
            jacket = (3, 6, 9)
        )
    )
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        placements;
        system_id = "cable_18kv_1000mm2_trefoil",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = earth,
        frequencies = p.frequencies
    )
end
