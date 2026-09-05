# Adapted from cigre_tb880_case0_132kv_630cu.json, CableBuilderSpec v1,
# and its neighbouring SeraingThermal/specs/materials.json.
# Preserve the declared solid-core radius, 0.01515 m, rather than rescaling it
# to the nominal 630 mm² name: its represented metal area is about 721.06 mm².
# Layer dimensions are literal source values, in metres. The trefoil, earth,
# temperature and frequency grid follow existing gauntlet operating defaults.
# This is a case derived from the supplied design, not a CIGRE reference result.
case_definition(
    :cable_132kv_cigre_tb880_case0_630cu_trefoil,
    (
        core_radius = case_parameter(
            :core_radius, 0.01515; tags = (:geometry, :cable_layer)
        ),
        inner_semicon_thickness = case_parameter(
            :inner_semicon_thickness, 0.0015; tags = (:geometry, :cable_layer)
        ),
        insulation_thickness = case_parameter(
            :insulation_thickness, 0.0155; tags = (:geometry, :cable_layer)
        ),
        outer_semicon_thickness = case_parameter(
            :outer_semicon_thickness, 0.0013; tags = (:geometry, :cable_layer)
        ),
        sheath_thickness = case_parameter(
            :sheath_thickness, 0.0008; tags = (:geometry, :cable_layer)
        ),
        jacket_thickness = case_parameter(
            :jacket_thickness, 0.0035; tags = (:geometry, :cable_layer)
        ),
        semicon_rho = case_parameter(
            :semicon_rho, 0.06; tags = (:material, :cable_layer)
        ),
        semicon_eps_r = case_parameter(
            :semicon_eps_r, 1000.0; tags = (:material, :cable_layer)
        ),
        xlpe_tan_delta = case_parameter(
            :xlpe_tan_delta, 0.001; tags = (:material, :cable_layer)
        ),
        formation_x = case_parameter(:formation_x, 0.0; tags = (:geometry, :system)),
        formation_y = case_parameter(:formation_y, -1.0; tags = (:geometry, :system)),
        formation_clearance_ratio = case_parameter(
            :formation_clearance_ratio, 0.2; tags = (:geometry, :system)
        ),
        line_length = case_parameter(:line_length, 1000.0; tags = (:operation, :system)),
        temperature = case_parameter(:temperature, 20.0; tags = (:operation, :system)),
        earth_rho = case_parameter(:earth_rho, 100.0; tags = (:material, :earth)),
        earth_eps_r = case_parameter(:earth_eps_r, 10.0; tags = (:material, :earth)),
        frequencies = case_parameter(
            :frequencies, collect(10.0 .^ range(0, stop = 6, length = 101));
            tags = (:operation, :frequency)
        )
    ),
    [
        "cable:1:core", "cable:1:sheath",
        "cable:2:core", "cable:2:sheath",
        "cable:3:core", "cable:3:sheath"
    ];
    description = "CIGRE TB 880 Case 0 — 132 kV 630 mm² Cu cables in trefoil"
) do p
    materials = LineCableModels.MaterialsLibrary(add_defaults = true)
    copper = LineCableModels.Material(materials, :copper)
    aluminum = LineCableModels.Material(materials, :aluminum)
    hdpe = LineCableModels.Material(materials, :pe)
    xlpe = LineCableModels.Material(
        kind = :insulator, rho = 1.97e14, eps_r = 2.5,
        tan_delta = p.xlpe_tan_delta
    )
    carbon_pe = LineCableModels.Material(
        kind = :semicon, rho = p.semicon_rho, eps_r = p.semicon_eps_r
    )
    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "cigre_tb880_case0_132kv_630cu",
        LineCableModels.Stack(
            LineCableModels.terminal(
                :core,
                LineCableModels.core(copper; r = p.core_radius, tag = :compacted_copper_conductor),
                LineCableModels.screen(carbon_pe; t = p.inner_semicon_thickness, tag = :inner_semiconductor),
                LineCableModels.insulation(xlpe; t = p.insulation_thickness, tag = :xlpe_insulation),
                LineCableModels.screen(carbon_pe; t = p.outer_semicon_thickness, tag = :outer_semiconductor)
            ),
            LineCableModels.terminal(
                :sheath,
                LineCableModels.sheath(aluminum; t = p.sheath_thickness, tag = :aluminum_sheath),
                LineCableModels.jacket(hdpe; t = p.jacket_thickness, tag = :hdpe_oversheath)
            )
        )
    )
    spacing = (2 + p.formation_clearance_ratio) * LineCableModels.outer_radius(design)
    placements = LineCableModels.trefoil(
        design;
        center = LineCableModels.at(p.formation_x, p.formation_y),
        spacing,
        connections = (core = (1, 3, 5), sheath = (2, 4, 6))
    )
    system = LineCableModels.build(
        LineCableModels.LineCableSystem, placements;
        system_id = "cable_132kv_cigre_tb880_case0_630cu_trefoil",
        line_length = p.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = p.temperature,
        earth_props = LineCableModels.homogeneous(rho = p.earth_rho, eps_r = p.earth_eps_r),
        frequencies = p.frequencies
    )
end
