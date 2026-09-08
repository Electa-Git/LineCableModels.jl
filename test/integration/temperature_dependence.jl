@testitem "Engine / shared temperature law reaches scalar, Gridspace, constants and export" tags=[:integration] begin
    copper = Material(:conductor,1.72e-8,1,1,20,0.004)
    dielectric = Material(:insulator,1e7,2.3,1,20,-0.003;tan_delta=0.025)
    function cable_with(metal, passive)
        build(CableDesign,"temperature-law",terminal(:core,
            core(metal;r=0.005),insulation(passive;t=0.005)))
    end
    function system_with(design)
        build(LineCableSystem,[design,design],[(0.0,-1.0),(0.1,-1.0)];
            connections=[Dict(:core=>1),Dict(:core=>2)])
    end
    design = cable_with(copper,dielectric)
    system = system_with(design)
    problem = LineParametersProblem(system;temperature=80.0,frequencies=[50.0,1000.0],earth_props=homogeneous(rho=100.0))
    calls = Ref(0)
    twice = (m,t,p,o,w) -> begin
        calls[] += 1
        2m.rho
    end
    const TD = LineCableModels.Materials.TemperatureDependent
    @eval LineCableModels.computation_options(
        ::LineCableModels.FormulaMethod{:default,typeof(TD.temperature_resistivity)},
        ::$(typeof(twice))) = (;)
    declaration = formula(:default;hooks=(contribution=twice,))
    selected = Formulation(temperature_dependence=declaration,
        insulation_admittance=:Ametani2004,options=(ideal_transposition=false,))
    reference_design = cable_with(
        Material(:conductor,2copper.rho,copper.eps_r,copper.mu_r),
        Material(:insulator,2dielectric.rho,dielectric.eps_r,dielectric.mu_r;
            tan_delta=dielectric.tan_delta))
    reference_problem = LineParametersProblem(system_with(reference_design);
        temperature=80.0,frequencies=problem.frequencies,earth_props=problem.earth_props)
    identity = Formulation(temperature_dependence=nothing,
        insulation_admittance=:Ametani2004,options=selected.options)
    reference = compute(reference_problem,identity)
    actual = compute(problem,selected)
    @test calls[] > 0
    @test actual.Z.values ≈ reference.Z.values rtol=2e-13
    @test actual.Y.values ≈ reference.Y.values rtol=2e-13
    @test details(actual).formulations.modified.temperature_dependence
    grid = Formulation(temperature_dependence=Grid((declaration,nothing)),
        insulation_admittance=:Ametani2004,options=selected.options)
    results = compute(problem,grid)
    @test results[1].Z.values == actual.Z.values
    @test results[1].Y.values == actual.Y.values
    unchanged = compute(problem,identity)
    @test results[2].Z.values == unchanged.Z.values
    @test results[2].Y.values == unchanged.Y.values
    @test !isapprox(results[1].Z.values,results[2].Z.values;rtol=1e-5)
    constants = compute(CableConstantsProblem(design;temperature=80.0),
        CableConstantsFormulation(temperature_dependence=declaration,insulation_admittance=:Ametani2004))
    reference_constants = compute(CableConstantsProblem(reference_design;temperature=80.0),
        CableConstantsFormulation(temperature_dependence=nothing,insulation_admittance=:Ametani2004))
    for request in (R,L,C,G)
        @test request(constants) ≈ request(reference_constants) rtol=2e-13
    end
    const IE = LineCableModels.ImportExport
    exported = only(IE._pscad_components(design,50.0,selected,80.0))
    expected = only(IE._pscad_components(reference_design,50.0,identity,80.0))
    @test exported.conductor.material.rho == expected.conductor.material.rho
    @test exported.dielectric.shunt_conductance ≈ expected.dielectric.shunt_conductance
    @test exported.dielectric.shunt_capacitance ≈ expected.dielectric.shunt_capacitance
    hot = LineParametersProblem(system;temperature=250.0,frequencies=[50.0],earth_props=problem.earth_props)
    @test_throws DomainError compute(hot,Formulation())
    @test all(isfinite,compute(hot,identity).Z)
    @test all(isfinite,compute(hot,selected).Z)
end
