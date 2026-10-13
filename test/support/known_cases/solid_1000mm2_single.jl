function known_case_problem(::Val{:solid_1000mm2_single};
        frequencies=10.0 .^ range(-1, 7; length=101))
    library = MaterialsLibrary(add_defaults=true)
    temperature = 20.0
    design = build(CableDesign, "solid_1000mm2", Stack(
        Group(:core, Region(:core_metal, Disk(sqrt(1000.0e-6 / π)),
            Material(library, :aluminum))),
        Region(:core_insulation, Shell(8.3e-3), Material(library, :xlpe))))
    system = build(LineCableSystem, design, Pose2(0.0, -1.0);
        connections=Dict(:core=>1), system_id="solid_1000mm2_single",
        line_length=1000.0)
    earth = homogeneous(rho=100.0, eps_r=10.0, mu_r=1.0)
    return LineParametersProblem(system; temperature, earth_props=earth, frequencies)
end
