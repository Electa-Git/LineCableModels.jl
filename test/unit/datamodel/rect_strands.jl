@testitem "DataModel / RectStrands / packing and area-preserving geometry" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures,
    MaterialFixtures
] begin
    r_in=0.01
    thickness=0.001
    width=0.002
    count=12
    r_ex=sqrt(r_in^2+count*width*thickness/π)
    layer=RectStrands(
        r_in,
        r_ex,
        thickness,
        width,
        count,
        12.0,
        copper_props
    )

    @test layer isa RectStrands{Float64}
    @test layer.cross_section ≈ count*width*thickness
    @test layer.thickness == layer.shape.thickness == thickness
    @test layer.width == layer.shape.width == width
    @test :width in propertynames(layer)
    @test maxfill(RectStrands, r_in, width) == floor(Int, 2π*r_in/width)
    @test validate(layer) === layer
    @test !hasproperty(Tubular(0.0, 0.01, copper_props), :width)

    @test_throws DomainError RectStrands(
        r_in,
        r_ex+0.001,
        thickness,
        width,
        count,
        12.0,
        copper_props
    )
    @test_throws DomainError RectStrands(
        r_in,
        r_ex,
        thickness,
        width,
        maxfill(RectStrands, r_in, width)+1,
        12.0,
        copper_props
    )

    material32=convert(Material{Float32}, copper_props)
    rin32=Float32(r_in)
    thick32=Float32(thickness)
    width32=Float32(width)
    rex32=sqrt(rin32^2+count*width32*thick32/(one(Float32)*π))
    layer32=RectStrands(
        rin32,
        rex32,
        thick32,
        width32,
        count,
        Float32(12),
        material32
    )
    converted=convert(LineCableModels.DataModel.AbstractConductorPart{Float64}, layer32)
    @test layer32 isa RectStrands{Float32}
    @test converted isa RectStrands{Float64}
    @test isapprox(
        converted.r_ex,
        sqrt(converted.r_in^2 +
             converted.num_wires*converted.width*converted.thickness/π);
        atol = sqrt(eps(Float64)),
        rtol = sqrt(eps(Float64))
    )
end
