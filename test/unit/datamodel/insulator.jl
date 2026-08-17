@testitem "DataModel / Insulator / exact geometry and base-state values" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures,
    MaterialFixtures
] begin
    layer=Insulator(0.01, 0.015, insulator_props)

    @test layer isa Insulator{Float64}
    @test layer.shunt_capacitance ≈ shunt_capacitance(0.01, 0.015, insulator_props.eps_r)
    @test layer.shunt_conductance ≈ shunt_conductance(0.01, 0.015, insulator_props.rho)
    @test layer.shunt_conductance >= 0
    @test validate(layer) === layer
    @test !hasproperty(layer, :temperature)
    @test_throws MethodError Insulator(0.01, 0.015, insulator_props; temperature = 80.0)

    @test_throws DomainError Insulator(-0.01, 0.015, insulator_props)
    @test_throws DomainError Insulator(0.015, 0.01, insulator_props)
    @test_throws MethodError Insulator(0.01, 0.015, nothing)

    material32=convert(Material{Float32}, insulator_props)
    layer32=Insulator(Float32(0.01), Float32(0.015), material32)
    @test layer32 isa Insulator{Float32}
    @test convert(LineCableModels.DataModel.AbstractInsulatorPart{Float64}, layer32) isa
          Insulator{Float64}
end
