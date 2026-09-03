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
    ];
    description = "640 kV 2000 mm² cable DC bipole"
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
    matrix = LineCableModels.Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )
    core_strand_area = π * p.core_strand_radius^2
    radius = p.core_strand_radius
    course_shapes = LineCableModels.Rectangle[]
    for count in p.ring_counts
        strand_width = prevfloat(2π * radius / count)
        strand_thickness = core_strand_area / strand_width
        outer = sqrt(
            radius^2 + count * strand_width * strand_thickness / π
        )
        push!(course_shapes,
            LineCableModels.Rectangle(strand_width, strand_thickness))
        radius = outer
    end
    core = LineCableModels.stranded(
        copper;
        center = LineCableModels.Disk(p.core_strand_radius),
        shape = Tuple(course_shapes),
        layers = length(p.ring_counts),
        n = p.ring_counts,
        lay = LineCableModels.LayRatio(p.core_lay_ratio),
        compact = LineCableModels.FillFactor(1.0),
        boundary = LineCableModels.Disk(radius)
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

    lead_outer = radius + p.lead_screen_thickness
    push!(parts,
        LineCableModels.Group(
            :sheath,
            LineCableModels.Region(
                :sheath_lead_screen,
                LineCableModels.Annulus(radius, lead_outer),
                lead
            )
        ))
    radius = lead_outer
    push!(parts,
        LineCableModels.Region(
            :sheath_inner,
            LineCableModels.Shell(p.inner_sheath_thickness),
            pe
        ))
    radius += p.inner_sheath_thickness

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
            :jacket_insulation,
            LineCableModels.Shell(p.jacket_thickness),
            pe
        ))
    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "640kV_2000mm2",
        LineCableModels.Stack(parts)
    )
    earth = LineCableModels.homogeneous(
        rho = p.earth_rho, eps_r = p.earth_eps_r, mu_r = 1.0
    )
    connections = [Dict(:core => 3index - 2, :sheath => 3index - 1, :jacket => 3index)
                   for index in 1:2]
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        fill(design, 2),
        [LineCableModels.Pose2(p.cable_x[index], p.cable_y) for index in 1:2];
        connections,
        system_id = "cable_640kv_2000mm2_bipole",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = earth,
        frequencies = p.frequencies
    )
end
