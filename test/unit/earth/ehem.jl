@testitem "Earth / explicit reduction and FrequencyDependent order through public compute" tags=[:unit] setup=[TestFixtures] begin
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
    law=(
        material, f, p, options,
        workspace)->begin
        push!(events, :fd)
        EP.EarthMaterial(material.rho/2, material.eps_r, material.mu_r)
    end
    reduction=(rho, epsilon, mu, model, pair, f,
        parameters, options, workspace)->begin
        push!(events, :ehem)
        @test length(rho)==length(model.layers)==3
        @test all(layer->layer in 1:3, pair.layers)
        EP.EarthMaterial(last(rho)^2/100, last(epsilon), last(mu))
    end
    @eval LineCableModels.computation_options(
        ::LineCableModels.FormulaMethod{
            :default, typeof(EP.FrequencyDependent.earth_material)},
        ::$(typeof(law))) = (;)
    @eval LineCableModels.computation_options(
        ::LineCableModels.FormulaMethod{
            :default, typeof(EP.EquivalentHomogeneous.equivalent_material)},
        ::$(typeof(reduction))) = (;)
    for order in (:before, :after)
        empty!(events)
        selected=Formulation(
            earth_properties = formula(:default; hooks = (contribution = law,)),
            earth_impedance = formula(:default;
                equivalent_earth = formula(:default; order,
                    hooks = (contribution = reduction,))),
            earth_admittance = formula(:default;
                equivalent_earth = formula(:default; order,
                    hooks = (contribution = reduction,))))
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
        @test details(result).formulations.effective.earth_properties === :default
        @test details(result).formulations.equivalent_earth.earth_impedance.identifier ===
              :default
        @test details(result).formulations.modified.earth_properties
        @test details(result).formulations.equivalent_earth.earth_impedance.modified
    end
end
