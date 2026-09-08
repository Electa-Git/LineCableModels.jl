@testitem "Gauntlet / PSCAD shared grammar and constitutive export" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LineCableModels
    using .GauntletSupport
    const P=GauntletSupport.PSCADBenchmarks
    copper=Material(:conductor, 1.72e-8, 1, 1, 20, 0.004)
    semicon=Material(:semicon, 1e4, 40; tan_delta = 0.02)
    dielectric=Material(:insulator, 1e8, 2.3; tan_delta = 0.03)
    design=build(CableDesign,
        "pscad-formula-probe",
        terminal(:core,
            solid(copper, Disk(0.004)), screen(semicon; t = 0.0005),
            insulation(dielectric; t = 0.002)))
    earth=homogeneous(rho = 100.0)
    problem(height)=LineParametersProblem(
        build(LineCableSystem, [design], [Pose2(0, height)];
            connections = [Dict(:core=>1)]);
        earth_props = earth,
        temperature = 60,
        frequencies = collect(10.0 .^ range(-1, 6; length = 101)))
    overhead, underground=problem(1.0), problem(-1.0)
    for selected_problem in (overhead, underground)
        resolved=Formulation(Val(:pscad), selected_problem, Formulation(:pscad))
        @test formula_id(resolved.methods.earth_impedance) === :default
        @test formula_id(resolved.definitions.earth_impedance) === :default
        @test P.pscad_setting(resolved, selected_problem).value == 2
    end
    @test P.pscad_setting(Formulation(:pscad), overhead).field === :EarthForm2
    @test P.pscad_setting(Formulation(:pscad), underground).field === :EarthForm
    @test_throws MethodError Formulation(Val(:pscad), underground,
        Formulation(:pscad; earth_impedance = :Gary1976))
    @test_throws MethodError Formulation(Val(:pscad), underground,
        Formulation(:pscad; earth_impedance = :Carson1926))
    @test_throws ArgumentError Formulation(Val(:pscad), underground,
        Formulation(:pscad; earth_properties = :CIGRE2019))
    @test all(isfinite, compute(underground, Formulation(earth_impedance = :Pollaczek1926)).Z)
    for key in keys(Formulation(:pscad).definitions)
        selection=NamedTuple{(key,)}((Grid((formula(:default), formula(:default))),))
        space=Formulation(:pscad; selection...)
        @test space isa Gridspace{P.PSCADFormulation}
        @test length(space) == 2
        @test all(item -> isconcretetype(typeof(item)), space)
    end
    product=Formulation(:pscad; earth_impedance = Grid((:default, :Saad1996)),
        insulation_admittance = Grid((:default, :Ametani2004)))
    zipped=Formulation(:pscad; earth_impedance = Grid((:default, :Saad1996)),
        insulation_admittance = Grid((:default, :Ametani2004)), combine = :zip)
    @test length(product) == 4
    @test length(zipped) == 2

    # Native reuse is determined by exported numerical inputs plus solver setting,
    # not by whether the requested selector happened to be :default.
    @test_throws ArgumentError compute(underground, P.PSCADFormulation[])
    requested=(Formulation(:pscad),
        Formulation(:pscad; earth_impedance = :Pollaczek1926),
        Formulation(:pscad; insulation_admittance = :Ametani2004))
    resolved=map(value->Formulation(Val(:pscad), underground, value), requested)
    prepared=[P._stage_pscad_project(underground, value) for value in resolved]
    try
        @test length(unique(getproperty.(prepared, :root))) == 3
        for project in prepared
            @test P._work_parts(joinpath(project.root, "outputs")) ==
                  ["pscad", underground.system.system_id, basename(project.root)]
        end
        projects=[read(value.staged, String) for value in prepared]
        @test projects[1] == projects[2]
        @test projects[1] != projects[3]
        @test P.pscad_setting(resolved[1], underground) ==
              P.pscad_setting(resolved[2], underground)
        @test P.pscad_setting(Formulation(:pscad; earth_impedance = :Saad1996), underground) !=
              P.pscad_setting(resolved[1], underground)
    finally
        foreach(value->rm(value.root; recursive = true), prepared)
    end

    @test P.formulation_record(Formulation(:pscad)).effective.internal_impedance === nothing
    baseline=P.formulation_record(Formulation(:pscad))
    alternative=P.formulation_record(Formulation(:pscad;
        earth_impedance = formula(:default; equivalent_earth = formula(:default)), insulation_admittance = :Ametani2004))
    @test typeof(baseline) === typeof(alternative)
    result=LineParameters(PhaseDomain, zeros(ComplexF64, 1, 1, 1),
        zeros(ComplexF64, 1, 1, 1), [50.0]; details = (formulations = baseline,))
    @test computation_details(typeof(Formulation(:pscad)), result) === details(result)
end

@testitem "Gauntlet / PSCAD catalogue follows native setting dispatch" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport
    const P=GauntletSupport.PSCADBenchmarks
    const EI=LineCableModels.Engine.EarthImpedance
    for placement in (Val(:overhead), Val(:underground), Val(:mixed))
        identifiers=P.formulas(placement)
        @test allunique(identifiers)
        @test :default ∈ identifiers
        @test identifiers == P.formulas(placement)
        for identifier in EI.formulas()
            @test (identifier in identifiers) ==
                  applicable(P.pscad_setting, Val(identifier), placement)
        end
    end
    # All native settings resolve a retained mathematical selection.
    @test :Gary1976 in P.formulas(Val(:overhead))
    @test :Carson1926 in P.formulas(Val(:overhead))
    @test :Pollaczek1926 in P.formulas(Val(:underground))
    @test_throws ArgumentError P.formulas(Val(:unsupported))
end
