function known_case_problem(::Val{:two_bare_wires};
        frequencies=10.0 .^ range(-1, 7; length=101))
    library = MaterialsLibrary(add_defaults=true)
    case = :two_bare_wires
    temperature = 20.0
    core = Group(:core, Region(:core_metal, Disk(0.0425), Material(library, :copper)))
    insulation = Material(kind=:insulator, rho=1.97e14, eps_r=2.3, mu_r=1.0,
        T0=20.0, alpha=0.0)
    stack = case === :two_bare_wires ? Stack(core) :
            Stack(core, Region(:core_insulation, Shell(1.0e-3), insulation))
    design = build(CableDesign, String(case), stack)
    system = build(LineCableSystem, [design, design],
        [Pose2(0.0, -1.0), Pose2(1.0, -1.0)];
        connections=[Dict(:core=>1), Dict(:core=>2)], system_id=String(case),
        line_length=1.0)
    earth = homogeneous(rho=0.1, eps_r=1.0, mu_r=1.0)
    return LineParametersProblem(system; temperature, earth_props=earth, frequencies)
end
