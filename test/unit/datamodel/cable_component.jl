@testitem "DataModel / CableComponent / eager equivalent flattening" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures,
    MaterialFixtures
] begin
    conductor=ConductorGroup(Tubular(0.0, 0.006, copper_props))
    insulation=InsulatorGroup(Insulator(0.006, 0.010, insulator_props))
    component=CableComponent("core", conductor, insulation)

    @test component isa CableComponent{Float64}
    @test component.id == "core"
    @test component.conductor_group === conductor
    @test component.insulator_group === insulation
    @test component.conductor_props.rho ≈ equivalent_rho(
        conductor.resistance,
        conductor.r_ex,
        conductor.r_in
    )
    @test component.insulator_props.eps_r ≈ equivalent_eps(
        insulation.shunt_capacitance,
        insulation.r_ex,
        insulation.r_in
    )
    @test component.conductor_props.T0 == component.insulator_props.T0 == 20.0
    @test validate(component) === component

    bad_insulation=InsulatorGroup(Insulator(0.007, 0.010, insulator_props))
    @test_throws DomainError CableComponent("bad", conductor, bad_insulation)

    conductor32=convert(ConductorGroup{Float32}, conductor)
    promoted=CableComponent("promoted", conductor32, insulation)
    @test promoted isa CableComponent{Float64}
    converted=convert(CableComponent{Float32}, component)
    @test converted isa CableComponent{Float32}
    @test convert(CableComponent{Float64}, component) === component

    flattened_conductor=ConductorGroup(component)
    flattened_insulator=InsulatorGroup(component)
    @test length(flattened_conductor.layers) == 1
    @test length(flattened_insulator.layers) == 1
    @test flattened_conductor.r_in == conductor.r_in
    @test flattened_conductor.r_ex == conductor.r_ex
    @test flattened_insulator.reference_frequency == insulation.reference_frequency
end
