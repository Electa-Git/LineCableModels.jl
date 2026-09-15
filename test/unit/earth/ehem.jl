@testitem "Earth / explicit reduction and FrequencyDependent order through public compute" tags=[:unit] setup=[TestFixtures,FormulaContractModels] begin
    const E=LineCableModels.Engine
    const EP=LineCableModels.Earth
    system=TestFixtures.three_phase_system()
    earth=build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 0.5), EP.EarthLayer(500.0, 20.0, 1.0)))
    problem=LineParametersProblem(system; earth_props = earth, frequencies = [50.0])
    @test_throws DimensionMismatch compute(problem, Formulation())
    reduced=compute(problem,
        Formulation(
            earth_impedance = formula(:default; equivalent_earth = formula(:default)),
            earth_admittance = formula(:default; equivalent_earth = formula(:default))))
    homogeneous_problem=LineParametersProblem(
        system; earth_props = homogeneous(rho = 500.0, eps_r = 20.0), frequencies = [50.0])
    @test reduced.Z.values ≈ compute(homogeneous_problem, Formulation()).Z.values
    events=Symbol[]
    law=FormulaContractModels.DispersiveEarth(scale=50.0,events=events)
    reduction=FormulaContractModels.SquaredBottomEarth(events)
    for order in (:before, :after)
        empty!(events)
        sequence=order === :before ? EP.EquivalentHomogeneous.BeforeFD(reduction) :
            EP.EquivalentHomogeneous.AfterFD(reduction)
        selected=Formulation(earth_properties=law,
            earth_impedance=formula(:default;equivalent_earth=sequence),
            earth_admittance=formula(:default;equivalent_earth=sequence))
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
        @test events == (order===:before ? repeat([:ehem, :fd], pairs) :
               vcat(fill(:fd, length(earth.layers)-1), fill(:ehem, pairs)))
        @test details(result).formulations.equivalent_earth.earth_impedance.order === order
        @test details(result).formulations.equivalent_earth.earth_admittance.order === order
        @test details(result).formulations.effective.earth_properties === :DispersiveEarth
        @test details(result).formulations.equivalent_earth.earth_impedance.identifier ===
              :SquaredBottomEarth
        @test !hasproperty(details(result).formulations,:modified)
    end
end
