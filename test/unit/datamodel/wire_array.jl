@testitem "DataModel / CircStrands / packing and exact geometry" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures,
    MaterialFixtures
] begin
    solid=CircStrands(0.0, 0.001, 0.001, 1, 0.0, copper_props)
    layer=CircStrands(0.001, 0.003, 0.001, 6, 10.0, copper_props)

    @test solid isa CircStrands{Float64, Int}
    @test solid.cross_section ≈ π*0.001^2
    @test layer.num_wires == maxfill(CircStrands, 0.001, 0.001) == 6
    @test layer.r_ex == layer.r_in+2layer.radius_wire
    @test layer.cross_section ≈ 6π*0.001^2
    @test layer.resistance > 0
    @test validate(layer) === layer
    @test !hasproperty(layer, :temperature)
    @test_throws MethodError CircStrands(
        0.001,
        0.003,
        0.001,
        6,
        10.0,
        copper_props;
        temperature = 80.0
    )

    @test_throws DomainError CircStrands(0.001, 0.004, 0.001, 6, 10.0, copper_props)
    @test_throws DomainError CircStrands(0.001, 0.003, 0.001, 7, 10.0, copper_props)
    @test_throws ArgumentError CircStrands(
        0.001,
        0.003,
        0.001,
        6,
        10.0,
        copper_props;
        lay_direction = 0
    )
    @test_throws ArgumentError maxfill(Tubular, 0.001, 0.001)
    @test !isdefined(LineCableModels.DataModel, :WireArray)
    @test !isdefined(LineCableModels.DataModel, :MaxFill)

    material32=convert(Material{Float32}, copper_props)
    layer32=CircStrands(
        Float32(0.001),
        Float32(0.003),
        Float32(0.001),
        6,
        Float32(10),
        material32
    )
    converted=convert(LineCableModels.DataModel.AbstractConductorPart{Float64}, layer32)
    @test layer32 isa CircStrands{Float32, Int}
    @test converted isa CircStrands{Float64, Int}
    @test isapprox(
        converted.r_ex,
        converted.r_in+2converted.radius_wire;
        atol = sqrt(eps(Float64)),
        rtol = sqrt(eps(Float64))
    )
end
