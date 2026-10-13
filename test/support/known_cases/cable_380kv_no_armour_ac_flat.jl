function known_case_problem(::Val{:cable_380kv_no_armour_ac_flat};
        frequencies=10.0 .^ range(-1, 7; length=101))
    library = MaterialsLibrary(add_defaults=true)
    strand_layers = 6
    strand_diameter = 3.66e-3
    core_lay_ratio = 15.0
    inner_semicon_thickness = 2.0e-3
    insulation_thickness = 20.0e-3
    outer_semicon_thickness = 1.5e-3
    water_blocking_thickness = 0.3e-3
    lead_sheath_thickness = 3.3e-3
    inner_pe_thickness = 3.0e-3
    jacket_thickness = 10.0e-3
    core_rho = 2.3853e-8
    core_mu_r = 1.01663
    xlpe_rho = 1.0e14
    xlpe_eps_r = 2.77781
    xlpe_mu_r = 1.59891
    lead_rho = 2.14e-7
    lead_mu_r = 0.999983
    pe_rho = 1.97e14
    pe_eps_r = 2.51862
    pp_rho = 1.0e15
    pp_eps_r = 2.8
    pp_mu_r = 1.1263
    cable_x = (-1.0, 0.0, 1.0)
    cable_y = -1.0
    line_length = 20_000.0
    voltage = 380.0
    temperature = 70.0
    earth_rho = 1_000.0
    earth_eps_r = 10.0

    core = Material(:conductor, core_rho, 1.0, core_mu_r, 20.0, 0.00393)
    semicon_inner = Material(:semicon, 1_000.0, 1_000.0)
    semicon_outer = Material(:semicon, 500.0, 1_000.0)
    xlpe = Material(:insulator, xlpe_rho, xlpe_eps_r, xlpe_mu_r)
    lead = Material(:conductor, lead_rho, 1.0, lead_mu_r, 20.0, 0.004)
    pe = Material(:insulator, pe_rho, pe_eps_r)
    pp = Material(:insulator, pp_rho, pp_eps_r, pp_mu_r)
    strand_radius = strand_diameter / 2
    strand_count = 1 + 3strand_layers * (strand_layers + 1)
    radius = sqrt(strand_count) * strand_radius
    packed_core = stranded(
        core;
        shape=Disk(strand_radius),
        lay=LayRatio(core_lay_ratio),
        compact=true,
        boundary=Disk(radius),
        fill=Material(kind=:insulator, rho=Inf)
    )
    parts = AbstractCablePart[packed_core]
    for (tag, layer_thickness, material) in (
        (:core_semicon_inner, inner_semicon_thickness, semicon_inner),
        (:core_insulation, insulation_thickness, xlpe),
        (:core_semicon_outer, outer_semicon_thickness, semicon_outer),
        (:core_water_blocking, water_blocking_thickness, semicon_outer)
    )
        push!(parts,
            Region(tag, Shell(layer_thickness), material))
        radius += layer_thickness
    end
    lead_outer = radius + lead_sheath_thickness
    push!(parts,
        Group(
            :sheath,
            Region(
                :sheath_lead,
                Annulus(radius, lead_outer),
                lead
            )
        ))
    push!(parts,
        Region(
            :sheath_inner_pe,
            Shell(inner_pe_thickness),
            pe
        ))
    push!(parts,
        Region(
            :sheath_outer_jacket,
            Shell(jacket_thickness),
            pp
        ))

    design = build(
        CableDesign,
        "cable_380kv_no_armour",
        Stack(parts);
        nominal_data=(
            designation_code="cable_380kv_no_armour",
            U0=voltage / sqrt(3),
            U=voltage,
            conductor_cross_section=1_600.0
        )
    )
    earth = homogeneous(rho=earth_rho, eps_r=earth_eps_r, mu_r=1.0)
    connections = [Dict(:core => 2index - 1, :sheath => 2index) for index in 1:3]
    system = build(
        LineCableSystem,
        fill(design, 3),
        [Pose2(cable_x[index], cable_y) for index in 1:3];
        connections,
        environment=earth,
        system_id="cable_380kv_no_armour_ac_flat",
        line_length=line_length
    )
    return LineParametersProblem(system; temperature, earth_props=earth, frequencies)
end
