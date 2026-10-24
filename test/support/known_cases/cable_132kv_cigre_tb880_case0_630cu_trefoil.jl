function known_case_problem(::Val{:cable_132kv_cigre_tb880_case0_630cu_trefoil};
        frequencies=10.0 .^ range(-1, 7; length=101))
    library = MaterialsLibrary(add_defaults=true)
    core_radius = 0.01515
    inner_semicon_thickness = 0.0015
    insulation_thickness = 0.0155
    outer_semicon_thickness = 0.0013
    sheath_thickness = 0.0008
    jacket_thickness = 0.0035
    semicon_rho = 0.06
    semicon_eps_r = 1000.0
    xlpe_tan_delta = 0.001
    formation_x = 0.0
    formation_y = -1.0
    formation_clearance_ratio = 0.2
    line_length = 1000.0
    temperature = 20.0
    earth_rho = 100.0
    earth_eps_r = 10.0

    copper = Material(library, :copper)
    aluminum = Material(library, :aluminum)
    hdpe = Material(library, :pe)
    xlpe = Material(
        kind=:insulator, rho=1.97e14, eps_r=2.5,
        tan_delta=xlpe_tan_delta
    )
    carbon_pe = Material(kind=:semicon, rho=semicon_rho, eps_r=semicon_eps_r)
    design = build(
        CableDesign,
        "cigre_tb880_case0_132kv_630cu",
        Stack(
            terminal(
                :core,
                LineCableModels.core(copper; r=core_radius, tag=:compacted_copper_conductor),
                screen(carbon_pe; t=inner_semicon_thickness, tag=:inner_semiconductor),
                LineCableModels.insulation(xlpe; t=insulation_thickness, tag=:xlpe_insulation),
                screen(carbon_pe; t=outer_semicon_thickness, tag=:outer_semiconductor)
            ),
            terminal(
                :sheath,
                sheath(aluminum; t=sheath_thickness, tag=:aluminum_sheath),
                jacket(hdpe; t=jacket_thickness, tag=:hdpe_oversheath)
            )
        )
    )
    spacing = (2 + formation_clearance_ratio) * outer_radius(design)
    placements = trefoil(
        design;
        center=at(formation_x, formation_y),
        spacing,
        connections=(core=(1, 3, 5), sheath=(2, 4, 6))
    )
    system = build(
        LineCableSystem, placements;
        system_id="cable_132kv_cigre_tb880_case0_630cu_trefoil",
        line_length=line_length
    )
    earth = homogeneous(rho=earth_rho, eps_r=earth_eps_r)
    return LineParametersProblem(system; temperature, earth_props=earth, frequencies)
end
