@testset "DataModel geometry" begin
    @testset "mapping scales real fields" begin
        placement = Placement(1.0, -2.0, 0.25)
        mapped = map(x -> 2x, placement)
        @test (mapped.x, mapped.y, mapped.φ) == (2.0, -4.0, 0.5)
        circle = map(x -> 2x, Circle(1.5e-3, placement))
        @test circle.r == 3.0e-3
        @test (circle.at.x, circle.at.y, circle.at.φ) == (2.0, -4.0, 0.5)
    end

    @testset "an offset around a part is an annulus" begin
        insulation = LineCableModels.resolve(Circle(1.5e-3), Offset(1.0e-3))
        @test (r_in(insulation), r_ex(insulation)) == (1.5e-3, 2.5e-3)
        jacket = LineCableModels.resolve(insulation, Offset(2.0e-3))
        @test r_ex(jacket) ≈ 4.5e-3
        @test_throws DomainError LineCableModels.resolve(insulation, Annulus(2.0e-3, 3.0e-3))
    end

    @testset "a circular pattern lays round wires on its circle" begin
        course = placements(CircularPattern(6; r = 3.0e-3), Circle(1.5e-3))
        @test length(course.placements) == 6
        @test all(placement -> hypot(placement.x, placement.y) ≈ 3.0e-3, course.placements)
        @test_throws DomainError placements(CircularPattern(7; r = 3.0e-3), Circle(1.5e-3))
    end

    @testset "a circular pattern at an offset lies against the outline inside" begin
        pattern = CircularPattern(6; r = Offset(0))
        @test LineCableModels.resolve(Circle(1.5e-3), pattern, Circle(1.5e-3)).r ≈ 3.0e-3
        @test LineCableModels.resolve(Circle(1.5e-3), pattern, Rectangle(3.0e-3, 3.0e-3)).r ≈
              3.0e-3
        clearance = CircularPattern(12; r = Offset(1.0e-3))
        @test LineCableModels.resolve(Circle(4.5e-3), clearance, Circle(1.5e-3)).r ≈ 7.0e-3
        @test_throws DomainError Offset(-1.0e-3)
    end

    @testset "a circular pattern bends rectangular strands into annular sectors of their area" begin
        strand = Rectangle(3.0e-3, 3.0e-3)
        sector = placements(CircularPattern(12; r = 6.0e-3), strand).shape
        @test r_in(sector) ≈ 4.5e-3 && r_ex(sector) ≈ 7.5e-3
        @test sector.span ≈ 0.5
        @test area(sector) ≈ area(strand)
        annulus = placements(CircularPattern(1; r = 3.0e-3), Rectangle(6π * 1.0e-3, 1.0e-3)).shape
        @test annulus.span ≈ 2π
        @test all(abs.(centroid(annulus)) .< 1.0e-15)
    end

    @testset "a placement locates a cross-section in its frame" begin
        wire = LineCableModels.resolve(Placement(0, -1), Circle(1.5e-3, Placement(3.0e-3, 0)))
        @test all(centroid(wire) .≈ (3.0e-3, -1.0))
        @test support(wire) ≈ hypot(3.0e-3, 1.0) + 1.5e-3
    end

    @testset "the nearest distance of a shape from the origin" begin
        @test distance(Circle(1.0e-3, Placement(3.0e-3, 0))) ≈ 2.0e-3
        @test distance(Circle(1.0e-3)) == 0
        @test distance(Annulus(2.0e-3, 3.0e-3)) == 2.0e-3
        sector = placements(CircularPattern(18; r = 9.0e-3), Rectangle(3.0e-3, 3.0e-3)).shape
        @test distance(sector) ≈ 7.5e-3
    end

    @testset "a cross-section with holes" begin
        hole = Circle(1.0e-3, Placement(2.0e-3, 0))
        shape = LineCableModels.DifferenceShape(Circle(4.0e-3), [hole])
        @test area(shape) ≈ 15π * 1.0e-6
        @test all(centroid(shape) .≈ (-2.0e-3 / 15, 0))
        @test (support(shape), distance(shape)) == (4.0e-3, 0)
        moved = LineCableModels.resolve(Placement(1.0e-3, 0), shape)
        @test all(centroid(moved) .≈ (1.0e-3 - 2.0e-3 / 15, 0))
    end

    @testset "a lay lengthens the helical path" begin
        helix = Helix(LayRatio(12))
        @test overlength(helix, 3.0e-3) ≈ sqrt(1 + (π / 12)^2)
        @test overlength(helix, 0.0) == 1
        @test pitch(helix, 3.0e-3) ≈ 72.0e-3
    end

    @testset "Float32 and mixed inputs" begin
        @test LineCableModels.resolve(Circle(1.5f-3), Offset(1.0f-3)) isa Annulus{Float32}
        @test Annulus(1, 2.0f0) isa Annulus{Float32}
        @test placements(CircularPattern(6; r = 3.0f-3), Circle(1.5f-3)).shape isa
              Circle{Float32}
    end
end
