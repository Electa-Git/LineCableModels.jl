@testitem "DataModel / coated wires / a coat fills its outer circle and adds two layers to the blueprint" tags=[:unit, :datamodel] setup=[CoatedArmour] begin
    const DM = LineCableModels.DataModel
    coat_material = Material(:insulator, 1.0e14, 2.3)
    design = CoatedArmour.armoured(; coat_material)
    regions = design.geometry.regions
    wires = filter(placed -> placed.terminal === :armor, regions)
    coats = filter(placed -> placed.source.tag === :armor_coat, regions)
    @test length(wires) == length(coats) == 36
    @test all(placed -> placed.terminal === nothing, coats)
    centers(shapes) = [coordinate for shape in shapes for coordinate in (shape.at.x, shape.at.y)]
    @test all(placed -> [placed.primitive.ri, placed.primitive.ro] ≈ [1.5e-3, 1.7e-3], coats)
    @test centers(placed.primitive for placed in coats) ≈ centers(placed.primitive for placed in wires)
    # A coated wire leaves one hole in the fill, the coat's outer circle.
    fill = only(filter(placed -> placed.source.tag === :armor_matrix_fill, regions))
    @test all(hole -> hole isa Disk && hole.r ≈ 1.7e-3, fill.primitive.holes)
    @test centers(fill.primitive.holes) ≈ centers(placed.primitive for placed in wires)
    @test area(fill.primitive) ≈ pi * ((23.2e-3 + 1.7e-3)^2 - 21.5e-3^2) - 36pi * 1.7e-3^2
    # Coats overlap on a ring whose bare wires would fit, and a coat must lie in its fill.
    @test_throws "ring members overlap" CoatedArmour.armoured(n = 45)
    @test_throws "enclosed geometry must fit inside the enclosure boundary" CoatedArmour.armoured(
        span = 1.5e-3)

    # The tube spans the metal. The coat material lies on both sides of it, as thick as the
    # coat.
    pitch = 23.2e-3
    blueprint = DM.flatten(design)
    sheath, armor = blueprint[2], blueprint[3]
    inner, outer = last(sheath.layers), first(armor.layers)
    @test [inner.r_in, inner.r_ex] ≈ [pitch - 1.7e-3, pitch - 1.5e-3]
    @test [armor.conductor.r_in, armor.conductor.r_ex] ≈ [pitch - 1.5e-3, pitch + 1.5e-3]
    @test [outer.r_in, outer.r_ex] ≈ [pitch + 1.5e-3, pitch + 1.7e-3]
    @test inner.material == outer.material == coat_material
    @test length(armor.layers) == 2
    @test last(armor.layers).r_in ≈ outer.r_ex
    equivalent = homogenize(design)
    layers = Dict(placed.source.tag => placed for placed in equivalent.geometry.regions)
    radii(placed) = [placed.primitive.ri, placed.primitive.ro]
    @test radii(layers[:sheath_dielectric_2]) ≈ [pitch - 1.7e-3, pitch - 1.5e-3]
    @test radii(layers[:armor_dielectric_1]) ≈ [pitch + 1.5e-3, pitch + 1.7e-3]
    @test layers[:armor_dielectric_1].source.material.eps_r == 2.3
end

@testitem "Engine / coated wires / the coaxial Z of a coated ring is that of its bare wires" tags=[:unit, :engine] setup=[CoatedArmour] begin
    coated = compute(CoatedArmour.problem(CoatedArmour.armoured()))
    # Bare wires of the same radius at the same pitch, with matrix in place of the coat.
    bare = compute(CoatedArmour.problem(CoatedArmour.armoured(coat = nothing, gap = 0.2e-3)))
    @test coated.Z ≈ bare.Z rtol=1e-14
    # The coat's permittivity enters the shunt admittance.
    @test !isapprox(coated.Y, bare.Y; rtol = 1e-3)
    # Float32 inputs build and compute.
    single = CoatedArmour.armoured(T = Float32, coat = 0.2f-3,
        coat_material = Material(:insulator, 1.0f14, 2.3f0), enclosed = false)
    result = compute(CoatedArmour.problem(single; frequencies = Float32[50, 1e4]))
    @test all(isfinite, result.Z) && all(isfinite, result.Y)
end

@testitem "ParametricBuilder / a Grid coat thickness builds, flattens and computes" tags=[:unit, :parametric] setup=[CoatedArmour] begin
    const DM = LineCableModels.DataModel
    designs = CoatedArmour.armoured(coat = Grid((0.1e-3, 0.2e-3)), enclosed = false)
    @test designs isa Gridspace{CableDesign}
    @test length(designs) == 2
    for (coat, design) in zip((0.1e-3, 0.2e-3), designs)
        armor = DM.flatten(design)[3]
        # The ring takes its pitch from the outline of the coated wire.
        pitch = 21.5e-3 + 1.5e-3 + coat
        @test armor.conductor.r_in ≈ pitch - 1.5e-3
        @test first(armor.layers).r_ex ≈ pitch + 1.5e-3 + coat
    end
    results = compute(CoatedArmour.problem(designs), Formulation())
    @test length(results) == 2
    @test all(result -> all(isfinite, result.Y), results)
end

@testitem "DataModel / Measurements / an uncertain coat permittivity reaches the blueprint and Y" tags=[:unit, :measurements] setup=[CoatedArmour] begin
    using Measurements
    const DM = LineCableModels.DataModel
    design = CoatedArmour.armoured(
        coat_material = Material(:insulator, 1.0e14, measurement(2.3, 0.1)))
    blueprint = DM.flatten(design)
    @test uncertainty(last(blueprint[2].layers).material.eps_r) ≈ 0.1
    @test uncertainty(first(blueprint[3].layers).material.eps_r) ≈ 0.1
    result = compute(CoatedArmour.problem(design))
    @test uncertainty(imag(result.Y[3, 3, 1])) > 0
end
