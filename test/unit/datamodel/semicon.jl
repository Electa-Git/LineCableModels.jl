@testitem "DataModel / Semicon / exact geometry and base-state values" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures,
    MaterialFixtures
] begin
    layer=Semicon(0.01, 0.012, semicon_props)

    @test layer isa Semicon{Float64}
    @test layer.shunt_capacitance ≈ shunt_capacitance(0.01, 0.012, semicon_props.eps_r)
    @test layer.shunt_conductance ≈ shunt_conductance(0.01, 0.012, semicon_props.rho)
    @test validate(layer) === layer
    @test !hasproperty(layer, :temperature)
    @test_throws MethodError Semicon(0.01, 0.012, semicon_props; temperature = 80.0)

    @test_throws DomainError Semicon(-0.01, 0.012, semicon_props)
    @test_throws DomainError Semicon(0.012, 0.01, semicon_props)
    @test_throws MethodError Semicon(0.01, 0.012, nothing)

    material32=convert(Material{Float32}, semicon_props)
    layer32=Semicon(Float32(0.01), Float32(0.012), material32)
    @test layer32 isa Semicon{Float32}
    @test convert(LineCableModels.DataModel.AbstractInsulatorPart{Float64}, layer32) isa
          Semicon{Float64}
end
