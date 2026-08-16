@testitem "DataModel / RectStrands / area and electrical invariants" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestNumerics
] begin
    material=Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    r_in=0.01
    thickness=0.001
    width=0.002
    num_wires=12
    layer=RectStrands(
        r_in,
        thickness,
        width,
        num_wires,
        12.0,
        material;
        temperature = 40.0,
        lay_direction = -1
    )

    metallic_area=num_wires*thickness*width
    annular_area=π*(layer.r_ex^2-layer.r_in^2)
    @test TestNumerics.isapprox_scaled(annular_area, metallic_area)
    @test layer.shape.cross_section == metallic_area
    @test layer.shape.thickness == layer.r_ex - layer.r_in
    @test layer.shape.num_wires == num_wires
    @test layer.shape.lay_direction == -1
    @test layer.resistance > 0
    @test layer.gmr > layer.r_in
    @test layer.gmr < layer.r_ex

    warmer=RectStrands(
        r_in, thickness, width, num_wires, 12.0, material; temperature = 80.0)
    @test warmer.resistance > layer.resistance

    max_fill=RectStrands(r_in, thickness, width, MaxFill(), 12.0, material)
    @test max_fill.shape.num_wires == floor(Int, 2π * r_in / width)
end

@testitem "DataModel / RectStrands / representative invalid geometry" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport
] begin
    material=Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    valid_tail=(0.001, 0.002, 12, 12.0, material)

    @test_throws ArgumentError RectStrands(-0.01, valid_tail...)
    @test_throws ArgumentError RectStrands(0.01, -0.001, 0.002, 12, 12.0, material)
    @test_throws ArgumentError RectStrands(0.01, 0.001, -0.002, 12, 12.0, material)
    @test_throws ArgumentError RectStrands(0.01, 0.001, 0.002, 100, 12.0, material)
    @test_throws ArgumentError RectStrands(
        0.01,
        0.001,
        0.002,
        12,
        12.0,
        material;
        lay_direction = 0
    )
    @test_throws ErrorException RectStrands(0.01, 0.001, 0.002, 1, 12.0, material)
end
