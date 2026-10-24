@testitem "DataModel / lifting / one family method lifts every positional constructor" tags=[:unit, :datamodel] begin
    const DM = LineCableModels.DataModel
    const EA = LineCableModels.Earth
    copper = Material(:conductor, 1.7241e-8)

    # Dispatch sends a positional call with a Grid argument to the family method of its
    # type, and each point calls the scalar constructor.
    sectors = Sector(Grid((2pi / 3, pi / 2)), 0.003, 0.02)
    @test sectors isa Gridspace{Sector}
    @test [sector.span for sector in sectors] == [2pi / 3, pi / 2]
    strips = DM.BentStrip(Grid((0.01, 0.012)), 0.02, 1.0)
    @test [strip.ri for strip in strips] == [0.01, 0.012]
    layers = EA.EarthLayer(Grid((10.0, 100.0)), 1.0, 1.0, Inf)
    @test [layer.rho for layer in layers] == [10.0, 100.0]
    materials = Material(:conductor, Grid((1.7e-8, 2.8e-8)))
    @test [material.rho for material in materials] == [1.7e-8, 2.8e-8]
    regions = Region(:core, Disk(Grid((0.01, 0.02))), copper; combine = :product)
    @test [region.primitive.r for region in regions] == [0.01, 0.02]

    # The scalar call is the scalar method, and a call outside every method still fails.
    @test @inferred(Sector(2pi / 3, 0.003, 0.02)) isa Sector{Float64}
    @test @inferred(Disk(0.01)) isa Disk{Float64}
    @test_throws MethodError Disk("radius")
    # A raw positional group declaration types its pattern, path and compaction.
    member = Region(:wire, Disk(0.001), copper)
    @test DM.Group(:core, Pose2(0, 0, 0), member, nothing, nothing, nothing) isa DM.Group
    grouped = DM.Group(:core, Grid((Pose2(0, 0, 0), Pose2(0.1, 0, 0))), member, nothing,
        nothing, nothing)
    @test grouped isa Gridspace{DM.Group}
    @test length(grouped) == 2
end

@testitem "DataModel / lifting / formations and homogenize lift once" tags=[:unit, :datamodel] begin
    copper = Material(:conductor, 1.7241e-8)
    xlpe = Material(:insulator, 1.97e14, 2.5)
    designs = build(CableDesign, "lifted", Stack(Group(:core,
        Region(:core, Disk(Grid((5e-3, 6e-3))), copper)), Region(:ins, Shell(2e-3), xlpe)))

    # A space of designs homogenizes point by point.
    homogenized = homogenize(designs)
    @test homogenized isa Gridspace{CableDesign}
    @test [design.cable_id for design in homogenized] == fill("lifted_equivalent", 2)
    # One method takes a design or a space of designs, and lifts its keyword too.
    @test homogenize(first(designs)) isa CableDesign
    renamed = homogenize(first(designs); new_id = Grid(("first", "second")))
    @test [design.cable_id for design in renamed] == ["first", "second"]

    # A formation zips its own parameters.
    paired = trefoil(first(designs); spacing = Grid((0.1, 0.2)), φ0 = Grid((0.0, 0.1)),
        connections = (core = (1, 2, 3),), combine = :zip)
    @test paired isa Gridspace{Vector}
    @test length(paired) == 2
end
