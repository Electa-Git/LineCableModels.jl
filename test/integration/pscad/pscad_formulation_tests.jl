@testitem "PSCAD / shared grammar and constitutive export" tags=[:integration] begin
    using LineCableModels
    const P=LineCableModels.PSCAD
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
    @test validate(overhead, Formulation(:pscad)) === overhead
    @test validate(underground, Formulation(:pscad)) === underground
    for selected_problem in (overhead, underground)
        resolved=Formulation(Val(:pscad), selected_problem, Formulation(:pscad))
        @test formula_id(resolved.methods.earth_impedance) === :default
        @test formula_id(resolved.definitions.earth_impedance) === :default
        @test all(control -> control.value == 2, P.pscad_setting(resolved, selected_problem).ground)
    end
    @test only(P.pscad_setting(Formulation(:pscad), overhead).interactions.earth_impedance).source == 1
    @test only(P.pscad_setting(Formulation(:pscad), underground).interactions.earth_impedance).source == 2
    @test_throws ArgumentError Formulation(Val(:pscad), underground,
        Formulation(:pscad; earth_impedance = :Gary1976))
    @test_throws ArgumentError Formulation(Val(:pscad), underground,
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
    root=mktempdir()
    prepared=[P._stage_pscad_project(underground, value,P.pscad_setting(value,underground),root) for value in resolved]
    try
        @test length(unique(getproperty.(prepared, :root))) == 3
        @test all(project -> dirname(project.root) == root,prepared)
        projects=[read(value.staged, String) for value in prepared]
        @test projects[1] == projects[2]
        @test projects[1] != projects[3]
        @test P.pscad_setting(resolved[1], underground).ground ==
              P.pscad_setting(resolved[2], underground).ground
        @test P.pscad_setting(Formulation(:pscad; earth_impedance = :Saad1996), underground) !=
              P.pscad_setting(resolved[1], underground)
    finally
        foreach(value->rm(value.root; recursive = true), prepared)
    end

    @test NamedTuple(Formulation(:pscad)).methods.internal_impedance.identifier === :default
    baseline=LineCableModels.computation_details(Formulation(:pscad))
    alternative=LineCableModels.computation_details(Formulation(:pscad;
        earth_impedance = formula(:default; equivalent_earth = formula(:default)), insulation_admittance = :Ametani2004))
    @test typeof(baseline) === typeof(alternative)
    result=LineParameters(PhaseDomain, zeros(ComplexF64, 1, 1, 1),
        zeros(ComplexF64, 1, 1, 1), [50.0]; details = (formulations = baseline,))
    @test computation_details(typeof(Formulation(:pscad)), result) === details(result)
end

@testitem "PSCAD / catalogue follows native setting dispatch" tags=[:integration] begin
    const P=LineCableModels.PSCAD
    const EI=LineCableModels.Engine.EarthImpedance
    for owner in (EI, LineCableModels.Engine.EarthAdmittance)
        equation = owner === EI ? P.earth_impedance : P.earth_potential_coefficient
        fallback = which(equation, Tuple{Val, Val, Val, Val, Val{:pscad}})
        for (kind, source, target) in ((:self, 1, 1), (:mutual, 1, 1),
                (:self, 2, 2), (:mutual, 2, 2), (:mutual, 1, 2), (:mutual, 2, 1))
            identifiers = P.formulas(owner, Val(kind), Val(source), Val(target))
            @test allunique(identifiers)
            @test :default in identifiers
            for identifier in owner.formulas()
                @test (identifier in identifiers) == (which(equation,
                    Tuple{Val{identifier}, Val{kind}, Val{source}, Val{target}, Val{:pscad}}) !== fallback)
            end
        end
        @test isempty(P.formulas(owner, Val(:self), Val(1), Val(2)))
        @test isempty(P.formulas(owner, Val(:mutual), Val(2), Val(3)))
    end
    @test :Gary1976 in P.formulas(EI, Val(:mutual), Val(1), Val(1))
    @test :Pollaczek1926 in P.formulas(EI, Val(:self), Val(2), Val(2))

end

@testitem "PSCAD / consumes complete homogeneous choices and preserves export settings" tags=[:integration] begin
    using EzXML
    const P = LineCableModels.PSCAD
    const E = LineCableModels.Engine
    copper = Material(:conductor, 1.72e-8, 1, 1, 20, 0.004)
    design = build(CableDesign, "mixed-native-contract", terminal(:core,
        solid(copper, Disk(0.004)), insulation(Material(:insulator, 1e14, 2.3); t = 0.002)))
    system = build(LineCableSystem, [design, design], [Pose2(0, 2), Pose2(1, -1)];
        connections = [Dict(:core => 1), Dict(:core => 2)])
    problem = LineParametersProblem(system; temperature = 60,
        earth_props = homogeneous(rho = 100.0), frequencies = [50.0])
    choices = (air = :Gary1976, earth = :Saad1996, mixed = :Lucca1994)
    selected = Formulation(:pscad; earth_impedance = choices)
    @test Formulation(Val(:pscad), problem, selected) === selected
    setting = P.pscad_setting(selected, problem)
    @test map(control -> control.value, setting.ground) == (EarthForm2 = 0, EarthForm = 3, EarthForm3 = 2)
    @test Set((r.formula, r.kind, r.source, r.target) for r in setting.interactions.earth_impedance) ==
        Set(((:Gary1976, :self, 1, 1), (:Saad1996, :self, 2, 2),
            (:Lucca1994, :mutual, 1, 2), (:Lucca1994, :mutual, 2, 1)))
    @test_throws ArgumentError Formulation(Val(:pscad), problem,
        Formulation(:pscad; earth_impedance = :Lucca1994))
    @test_throws ArgumentError Formulation(Val(:pscad), problem,
        Formulation(:pscad; earth_admittance = :Pollaczek1926))
    for invalid_choices in ((air = :Pollaczek1926, earth = :Saad1996, mixed = :Lucca1994),
            (air = :Gary1976, earth = :Carson1926, mixed = :Lucca1994))
        @test_throws ArgumentError Formulation(Val(:pscad), problem,
            Formulation(:pscad; earth_impedance = invalid_choices))
    end
    @test_throws ArgumentError Formulation(Val(:pscad), problem,
        Formulation(:pscad; earth_impedance = (
            air = formula(:Gary1976; options = (integration = (method = :quad,),)),
            earth = :Saad1996, mixed = :Lucca1994)))
    staged = P._stage_pscad_project(problem, selected, setting, mktempdir())
    try
        document = EzXML.readxml(staged.staged)
        ground = only(EzXML.findall("//User[@defn='master:Line_Ground']", document))
        fields = Dict(node["name"] => node["value"] for node in EzXML.findall("./paramlist/param", ground))
        @test all(parse(Float64, fields[string(name)]) == control.value for (name, control) in pairs(setting.ground))
        frequency = only(EzXML.findall("//User[@defn='master:Line_FrePhase_Options']", document))
        fields = Dict(node["name"] => node["value"] for node in EzXML.findall("./paramlist/param", frequency))
        @test parse(Float64, fields["enablf"]) == 1
        for cable in EzXML.findall("//User[@defn='master:Cable_Coax']", document)
            fields = Dict(node["name"] => node["value"] for node in EzXML.findall("./paramlist/param", cable))
            @test parse(Float64, fields["RHOC"]) ≈ 1.72e-8 * (1 + 0.004 * 40)
            @test parse(Float64, fields["LT1"]) == 0
        end
    finally
        rm(staged.root; recursive = true)
    end
    record = NamedTuple(selected)
    @test map(value -> value.identifier,record.requested.earth_impedance) == choices
    @test record.methods.earth_admittance.identifier === :default
    @test record.requested.earth_impedance.air.identifier === :Gary1976
    vertical = LineParametersProblem(build(LineCableSystem, [design, design],
        [Pose2(0, -1), Pose2(0, -2)]; connections = [Dict(:core => 1), Dict(:core => 2)]);
        earth_props = homogeneous(rho = 100.0), frequencies = [50.0])
    native = P.pscad_setting(Formulation(:pscad; earth_impedance = :WedepohlWilcox1973), vertical)
    @test native.ground.EarthForm.readback == "WEDEPOHL"
    @test Set(row.kind for row in native.interactions.earth_impedance) == Set((:self, :mutual))
end
