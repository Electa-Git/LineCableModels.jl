function build_case(::Val{:cable_18kv_1000mm2_trefoil_homogenized}, p)
    problem = build_case(Val(:cable_18kv_1000mm2_trefoil), p)
    source = problem.system
    designs = map(LineCableModels.homogenize, source.designs)
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        designs,
        source.positions;
        connections = source.connections,
        environment = source.environment,
        system_id = "cable_18kv_1000mm2_trefoil_homogenized",
        line_length = source.line_length
    )
    return LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = problem.temperature,
        earth_props = problem.earth_props,
        frequencies = problem.frequencies,
        Γ = problem.Γ
    )
end

let original = include("cable_18kv_1000mm2_trefoil.jl")
    case_definition(Base.Fix1(build_case, Val(:cable_18kv_1000mm2_trefoil_homogenized)),
        :cable_18kv_1000mm2_trefoil_homogenized,
        original.parameters,
        original.port_order;
        description = "18 kV 1000 mm² cables in trefoil — homogenized assembly",
        assets = ("cable_18kv_1000mm2_trefoil.jl",)
    )
end
