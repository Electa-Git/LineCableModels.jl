@testitem "Earth / explicit reduction and FrequencyDependent order through public compute" tags=[:unit, :parametric, :slow] setup=[
    TestFixtures, FormulaFixtures] begin
    const E=LineCableModels.Engine
    const EP=LineCableModels.Earth
    system=TestFixtures.three_phase_system()
    earth=build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 0.5), EP.EarthLayer(500.0, 20.0, 1.0)))
    problem=LineParametersProblem(system; earth_props = earth, frequencies = [50.0])
    reduced=compute(problem,
        Formulation(
            earth_impedance = formula(:default; equivalent_earth = formula(:default)),
            earth_admittance = formula(:default; equivalent_earth = formula(:default))))
    # A homogeneous formula on a layered earth uses the `:default` reduction unless one is given.
    unspecified=compute(problem, Formulation())
    @test unspecified.Z.values == reduced.Z.values && unspecified.Y.values == reduced.Y.values
    homogeneous_problem=LineParametersProblem(
        system; earth_props = homogeneous(rho = 500.0, eps_r = 20.0), frequencies = [50.0])
    @test reduced.Z.values ≈ compute(homogeneous_problem, Formulation()).Z.values
    events=Symbol[]
    law=FormulaFixtures.DispersiveEarth(scale = 50.0, events = events)
    reduction=FormulaFixtures.SquaredBottomEarth(events)
    for order in (:before, :after)
        empty!(events)
        empty!(reduction.workspaces)
        empty!(law.seen)
        sequence=order===:before ? EP.EquivalentHomogeneous.BeforeFD(reduction) :
                 EP.EquivalentHomogeneous.AfterFD(reduction)
        selected=Formulation(earth_properties = law,
            earth_impedance = formula(:default; equivalent_earth = sequence),
            earth_admittance = formula(:default; equivalent_earth = sequence))
        result=compute(problem, selected)
        expected_rho=order===:before ? 1250.0 : 625.0
        expected=compute(
            LineParametersProblem(
                system; earth_props = homogeneous(rho = expected_rho, eps_r = 20.0),
                frequencies = [50.0]),
            Formulation())
        @test result.Z.values ≈ expected.Z.values
        pairs=count(==(:ehem), events)
        @test pairs > 0
        workspace=first(reduction.workspaces)
        @test workspace isa E.LineParametersWorkspace
        @test all(w -> w === workspace, reduction.workspaces)
        @test all(record -> last(record) === workspace, law.seen)
        @test events == (order===:before ? repeat([:ehem, :fd], pairs) :
               vcat(fill(:fd, length(earth.layers)-1), fill(:ehem, pairs)))
        @test details(result).data.formulations.methods.earth_impedance.equivalent_earth.order ===
              order
        @test details(result).data.formulations.methods.earth_admittance.equivalent_earth.order ===
              order
        @test details(result).data.formulations.methods.earth_properties.identifier ===
              :DispersiveEarth
        @test details(result).data.formulations.methods.earth_impedance.equivalent_earth.rule.identifier ===
              :SquaredBottomEarth
    end
end

@testitem "Earth / missing equivalent equation fails on evaluation, not declaration" tags=[:unit, :engine] begin
    const EH = LineCableModels.Earth.EquivalentHomogeneous
    struct UnimplementedReduction <: EH.AbstractRule
        parameters::NamedTuple
        options::FormulationOptions
    end
    LineCableModels.formula_id(::UnimplementedReduction) = :UnimplementedReduction
    LineCableModels.formulation_options(::LineCableModels.FormulaMethod{
        UnimplementedReduction, typeof(EH.equivalent_material)}) = FormulationOptions()
    rule = UnimplementedReduction((;), FormulationOptions())
    pair = LineCableModels.Engine.EarthPair(1, 1, (-1.0, -1.0), 0.0, (2, 2); radius=0.01)
    model = homogeneous(rho=100.0, eps_r=10.0)
    binding = only(LineCableModels.Commons.bindings(rule, (pair,)))
    @test binding.equation.selection === rule
    @test_throws r"equivalent_material :UnimplementedReduction.*source in layer 2 and target in layer 2" rule(
        [Inf, 100.0], [1.0, 10.0], [1.0, 1.0], model, pair, 50.0; binding)
end
