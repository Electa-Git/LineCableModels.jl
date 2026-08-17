@testitem "DataModel / Strip / exact geometry and base-state values" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures,
    MaterialFixtures
] begin
    layer=Strip(0.01, 0.012, 0.05, 10.0, copper_props)

    @test layer isa Strip{Float64}
    @test layer.thickness ≈ 0.002
    @test layer.cross_section ≈ layer.thickness*layer.width
    @test layer.resistance > 0
    @test layer.pitch_length > 0
    @test validate(layer) === layer
    @test !hasproperty(layer, :temperature)
    @test_throws MethodError Strip(
        0.01,
        0.012,
        0.05,
        10.0,
        copper_props;
        temperature = 80.0
    )

    @test_throws DomainError Strip(-0.01, 0.012, 0.05, 10.0, copper_props)
    @test_throws DomainError Strip(0.012, 0.01, 0.05, 10.0, copper_props)
    @test_throws DomainError Strip(0.01, 0.012, -0.05, 10.0, copper_props)
    @test_throws ArgumentError Strip(
        0.01,
        0.012,
        0.05,
        10.0,
        copper_props;
        lay_direction = 0
    )

    material32=convert(Material{Float32}, copper_props)
    layer32=Strip(
        Float32(0.01),
        Float32(0.012),
        Float32(0.05),
        Float32(10),
        material32
    )
    @test layer32 isa Strip{Float32}
    @test convert(LineCableModels.DataModel.AbstractConductorPart{Float64}, layer32) isa
          Strip{Float64}
end
