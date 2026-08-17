@testitem "DataModel / Tubular / exact geometry and base-state values" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures,
    MaterialFixtures
] begin
    layer=Tubular(0.01, 0.02, copper_props)

    @test layer isa Tubular{Float64}
    @test layer.r_in == 0.01
    @test layer.r_ex == 0.02
    @test layer.cross_section ≈ π*(0.02^2-0.01^2)
    @test layer.resistance ≈ tubular_resistance(
        0.01,
        0.02,
        copper_props.rho,
        copper_props.alpha,
        copper_props.T0,
        copper_props.T0
    )
    @test layer.gmr ≈ tubular_gmr(0.02, 0.01, copper_props.mu_r)
    @test validate(layer) === layer
    @test !hasproperty(layer, :temperature)
    @test_throws MethodError Tubular(0.01, 0.02, copper_props; temperature = 80.0)

    @test_throws DomainError Tubular(-0.01, 0.02, copper_props)
    @test_throws DomainError Tubular(0.02, 0.01, copper_props)
    @test_throws DomainError Tubular(0.01, Inf, copper_props)
    @test_throws MethodError Tubular("0.01", 0.02, copper_props)

    material32=convert(Material{Float32}, copper_props)
    layer32=Tubular(Float32(0.01), Float32(0.02), material32)
    @test layer32 isa Tubular{Float32}
    @test convert(LineCableModels.DataModel.AbstractConductorPart{Float64}, layer32) isa
          Tubular{Float64}
    mixed=Tubular(Float32(0.01), 0.02, material32)
    @test mixed isa Tubular{Float64}
end
