@testitem "DataModel / CablePosition / mapping and geometry" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures
] begin
    design=TestFixtures.mv_cable_design()
    mapping=Dict("core"=>1, "sheath"=>0, "jacket"=>0)
    position=CablePosition(design, 0.0, -1.0, mapping)

    @test position isa CablePosition{Float64}
    @test position.design_data === design
    @test position.conn == [1, 0, 0]
    @test validate(position) === position
    @test outer_radius(design) > 0

    @test_throws KeyError CablePosition(design, 0.0, -1.0, Dict("missing"=>1))
    @test_throws DomainError CablePosition(design, 0.0, 0.0, mapping)
    @test_throws DomainError CablePosition(design, 0.0, -outer_radius(design)/2, mapping)

    position32=convert(CablePosition{Float32}, position)
    @test position32 isa CablePosition{Float32}
    @test position32.conn == position.conn
    @test position32.conn !== position.conn
    @test convert(CablePosition{Float64}, position) === position

    default_mapping=CablePosition(design, 0.0, -1.0)
    @test default_mapping.conn == [1, 0, 0]
end
