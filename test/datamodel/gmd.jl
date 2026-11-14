using LineCableModels.BaseGeometry: Contour, gmd
using LineCableModels.DataModel: ArcPolygon

@testset "Shape geometric mean distance" begin
    # Float64 bounds measured with gmd at τ = 1e-6 on Julia 1.12.7, x86_64 Linux.
    disk = Circle(3e-3)
    for (A, B) in ((disk, disk),
        (Annulus(1e-3, 3e-3), Annulus(1e-3, 3e-3)),
        (Annulus(2.7e-3, 3e-3), Annulus(2.7e-3, 3e-3)),
        (Annulus(4e-3, 5e-3), disk),
        (Annulus(4e-3, 5e-3), Rectangle(3e-3, 2e-3, Placement(0.5e-3, 0.3e-3))),
        (Annulus(4e-3, 5e-3), ArcPolygon(Contour(Rectangle(3e-3, 2e-3)).outer)),
        (Annulus(4e-3, 5e-3), Rectangle(3e-3, 2e-3, Placement(8e-3, 0))),
        (disk, Rectangle(3e-3, 2e-3, Placement(8e-3, 0))),
        (disk, Circle(2e-3, Placement(8e-3, 1e-3))))
        @test abs(log(gmd(A, B)) - log(gmd(Contour(A), Contour(B)))) <= 1e-9
        @test gmd(A, B) ≈ gmd(B, A) rtol = 1e-12
    end
    @test isempty(Test.detect_ambiguities(LineCableModels.BaseGeometry,
        LineCableModels.DataModel; recursive = true))

    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1e14, 2.5)
    bar = Region(:core_bar, Rectangle(10e-3, 5e-3), copper)
    plated = build(CableDesign, "plated",
        Terminal(:core, Stack(bar, Region(:core_plating, Offset(0.5e-3), copper))),
        Region(:core_insulation, Offset(1e-3), polymer))
    bare = build(CableDesign, "bare", Terminal(:core, bar),
        Region(:core_insulation, Offset(1.5e-3), polymer))
    @test only(LineCableModels.flatten(plated).conductors).material.mu_r ≈ 2.050 rtol = 1e-3
    @test only(LineCableModels.flatten(bare).conductors).material.mu_r ≈ 2.423 rtol = 1e-3
    shapes = [plated.regions[i].shape for i in only(plated.terminals).regions]
    shares = area.(shapes) ./ sum(area, shapes)
    weighted = sum(shares[i] * shares[j] * log(gmd(shapes[i], shapes[j]))
    for i in eachindex(shapes), j in eachindex(shapes))
    outline = Contour(last(shapes))
    outer = Contour(outline.outer, empty(outline.holes))
    @test abs(weighted - log(gmd(outer, outer))) <= 1e-9

    # Each copy is a rigid image. Separate own terms differ only by numerical roundoff.
    reference = Dict(REFERENCES)["coaxial_rectangular_strands"]
    design = only(reference.system.designs)
    every = deepcopy(design)
    empty!(every.repeats)
    repeated, separate = LineCableModels.flatten(design), LineCableModels.flatten(every)
    for (A, B) in zip(repeated.conductors, separate.conductors)
        @test A.r_in ≈ B.r_in rtol = 1e-12
        @test A.r_ex ≈ B.r_ex rtol = 1e-12
        @test A.material.rho ≈ B.material.rho rtol = 1e-12
        @test A.material.mu_r ≈ B.material.mu_r rtol = 1e-12
    end
end
