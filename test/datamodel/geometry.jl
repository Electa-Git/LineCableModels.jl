using LineCableModels.DataModel: Sector, ArcPolygon, Curve, Segment, Arc

@testset "DataModel geometry" begin
    @testset "mapping scales real fields" begin
        placement = Placement(1.0, -2.0, 0.25)
        mapped = map(x -> 2x, placement)
        @test (mapped.x, mapped.y, mapped.φ) == (2.0, -4.0, 0.5)
        circle = map(x -> 2x, Circle(1.5e-3, placement))
        @test circle.r == 3.0e-3
        @test (circle.at.x, circle.at.y, circle.at.φ) == (2.0, -4.0, 0.5)
        pattern = CircularPattern(6; r = Offset(1.0e-3), φ0 = 0.25, span = π,
            gap_frac = 0.1)
        mapped = map(x -> 2x, pattern)
        @test (mapped.n, mapped.r.t, mapped.φ0, mapped.span, mapped.gap_frac) ==
              (6, 2.0e-3, 0.5, 2π, 0.2)
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
        shape = LineCableModels.PerforatedShape(Circle(4.0e-3), [hole])
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
        @test pitch(helix, 0.0) === Inf
    end

    @testset "mixed inputs" begin
        @test Annulus(1, 2.0) isa Annulus{Float64}
    end

    @testset "the five sector configurations and their bounds" begin
        at = Placement(0, 0)
        S1 = Sector(5.0e-3, 10.0e-3, π / 3, 0, 0, at)
        S2 = Sector(5.0e-3, 10.0e-3, π / 2, 0, 1.0e-3, at)
        S3 = Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, at)
        S4 = Sector(0, 0.03115, π / 3, 0, 0.04 * 0.03115, at)
        S5 = Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 0, at)
        for sector in (S1, S2, S3, S4, S5)
            @test sector isa Sector{Float64}
        end
        let sector = S3
            @test_throws DomainError Sector(-1.0e-3, sector.ro, sector.span,
                sector.apex, sector.fillet, at)
            @test_throws DomainError Sector(Inf, sector.ro, sector.span,
                sector.apex, sector.fillet, at)
            @test_throws DomainError Sector(sector.ri, sector.ri, sector.span,
                sector.apex, sector.fillet, at)
            @test_throws DomainError Sector(sector.ri, Inf, sector.span,
                sector.apex, sector.fillet, at)
            @test_throws DomainError Sector(sector.ri, sector.ro, 0,
                sector.apex, sector.fillet, at)
            @test_throws DomainError Sector(sector.ri, sector.ro, 2π + 0.1,
                sector.apex, sector.fillet, at)
            @test_throws DomainError Sector(sector.ri, sector.ro, sector.span,
                -1.0e-3, sector.fillet, at)
            @test_throws DomainError Sector(sector.ri, sector.ro, sector.span,
                Inf, sector.fillet, at)
            @test_throws DomainError Sector(sector.ri, sector.ro, sector.span,
                sector.ro * cos(sector.span / 2) + 1.0e-3, sector.fillet, at)
            @test_throws DomainError Sector(sector.ri, sector.ro, sector.span,
                sector.apex, -1.0e-3, at)
            @test_throws DomainError Sector(sector.ri, sector.ro, sector.span,
                sector.apex, Inf, at)
            @test_throws DomainError Sector(sector.ri, sector.ro, sector.span,
                sector.apex, sector.ro, at)
            @test_throws DomainError Sector(0, sector.ro, sector.span,
                0, 1.01sector.ro / (csc(sector.span / 2) + 1), at)
            @test_throws DomainError LineCableModels.resolve(sector,
                Offset(sector.apex * sin(sector.span / 2) + 1.0e-3))
        end
    end

    @testset "sharp sectors retain their closed forms" begin
        at = Placement(0, 0)
        S1 = Sector(5.0e-3, 10.0e-3, π / 3, 0, 0, at)
        @test area(S1) ≈ S1.span * (S1.ro^2 - S1.ri^2) / 2 rtol = 1e-12
        d = 4sin(S1.span / 2) * (S1.ro^3 - S1.ri^3) /
            (3S1.span * (S1.ro^2 - S1.ri^2))
        @test collect(centroid(S1)) ≈ [d, 0] rtol = 1e-12
        @test perimeter(S1) ≈ S1.span * (S1.ro + S1.ri) + 2(S1.ro - S1.ri) rtol = 1e-12
        S5 = Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 0, at)
        α = S5.span / 2
        g = S5.apex * sin(α)
        a, r = S5.apex, S5.ro
        A = α * (r^2 - a^2) -
            (r^2 * asin(g / r) + g * sqrt(r^2 - g^2) -
             a^2 * asin(g / a) - g * sqrt(a^2 - g^2))
        M = 2sin(α) * ((r^2 - g^2)^(3 / 2) - (a^2 - g^2)^(3 / 2)) / 3 -
            g * cos(α) * (r^2 - a^2)
        @test area(S5) ≈ A rtol = 1e-12
        @test collect(centroid(S5)) ≈ [M / A, 0] rtol = 1e-12
        @test perimeter(S5) ≈ 2r * (α - asin(g / r)) +
                              2(sqrt(r^2 - g^2) - a * cos(α)) rtol = 1e-12
        placed = Placement(0.02, -0.01)
        d = hypot(placed.x, placed.y)
        @test distance(LineCableModels.resolve(placed, S1)) ≈
              max(0, d - S1.ro, S1.ri - d) rtol = 1e-12
        @test distance(LineCableModels.resolve(placed, S5)) ≈
              max(0, d - S5.ro, S5.apex - d) rtol = 1e-12
    end

    @testset "rounded sectors agree with independent chords" begin
        at = Placement(0, 0)
        S2 = Sector(5.0e-3, 10.0e-3, π / 2, 0, 1.0e-3, at)
        S3 = Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, at)
        S4 = Sector(0, 0.03115, π / 3, 0, 0.04 * 0.03115, at)
        for (sector, A, x, y, P) in (
            (S2, 5.8070608729313431e-05, 0.0070197713396683358,
                -1.7908672166583439e-79, 0.031887982680126824),
            (S3, 8.3590913574623174e-05, 0.006094058395038585,
                -7.9958233757167357e-80, 0.035513403959788118),
            (S4, 0.00050628577182981569, 0.019860545479483984,
                -2.3901619677833067e-79, 0.092091223162218932))
            @test area(sector) ≈ A rtol = 1e-6
            @test collect(centroid(sector)) ≈ [x, y] rtol = 1e-6
            @test perimeter(sector) ≈ P rtol = 1e-6
        end
        for (φ, value) in zip((0, π / 2, π, 3π / 2, 1.0),
            (0.01, 0.0072856175163482457, -0.0013867513459481288,
                0.0072856175163482457, 0.0097329367889796402))
            @test support(S3, φ) ≈ value rtol = 1e-6
        end
    end

    @testset "the fillets change the wedge first moments" begin
        at = Placement(0, 0)
        S3 = Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, at)
        S5 = Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 0, at)
        upper_center = (0.0062270827117790768, 0.0057856236397351443)
        upper_back = (0.0073259796609165609, 0.0068066160467472286)
        upper_side = (0.0049280446061024188, 0.0065356236397351443)
        tip_center = (0.0028867513459481288, 0.0)
        tip_upper = (0.0015877132402714709, 0.00075)
        tip_lower = (0.0015877132402714709, -0.00075)
        lower_center = (0.0062270827117790768, -0.0057856236397351443)
        lower_side = (0.0049280446061024188, -0.0065356236397351443)
        lower_back = (0.0073259796609165609, -0.0068066160467472286)
        α = S3.span / 2
        β = α - asin(S3.apex * sin(α) / S3.ro)
        γ = atan(upper_back[2], upper_back[1])
        vertex_upper = (S3.ro * cos(β), S3.ro * sin(β))
        vertex_lower = (vertex_upper[1], -vertex_upper[2])
        upper = ArcPolygon(Curve{Float64}[
        Arc(Placement(0, 0, γ), S3.ro, β - γ),
        Segment(vertex_upper, upper_side),
        Arc(Placement(upper_center..., α + π / 2), S3.fillet, γ - α - π / 2)])
        tip = ArcPolygon(Curve{Float64}[
        Segment(tip_upper, (S3.apex, 0.0)),
        Segment((S3.apex, 0.0), tip_lower),
        Arc(Placement(tip_center..., 3π / 2 - α), S3.fillet, S3.span - π)])
        lower = ArcPolygon(Curve{Float64}[
        Arc(Placement(0, 0, -β), S3.ro, β - γ),
        Arc(Placement(lower_center..., -γ), S3.fillet, γ - α - π / 2),
        Segment(lower_side, vertex_lower)])
        moments = area(S5) .* centroid(S5)
        for corner in (upper, tip, lower)
            moments = moments .- area(corner) .* centroid(corner)
        end
        @test collect(area(S3) .* centroid(S3)) ≈ collect(moments) rtol = 1e-12
        @test centroid(S3) != centroid(S5)
    end

    @testset "a wedge offset satisfies Steiner's formula" begin
        S3 = Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(0, 0))
        t = 0.5e-3
        offset = LineCableModels.resolve(S3, Offset(t))
        expected = Sector(max(S3.ri - t, 0), S3.ro + t, S3.span,
            S3.apex - t / sin(S3.span / 2), S3.fillet + t, S3.at)
        for name in fieldnames(typeof(expected))
            @test getfield(offset, name) == getfield(expected, name)
        end
        @test area(offset) ≈ area(S3) + perimeter(S3) * t + π * t^2 rtol = 1e-12
    end

    @testset "a placement preserves rounded sector moments" begin
        S3 = Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(0, 0))
        at = Placement(0.02, -0.01, π / 5)
        placed = LineCableModels.resolve(at, S3)
        @test area(placed) ≈ area(S3) rtol = 1e-12
        @test collect(centroid(placed)) ≈ collect(at(centroid(S3))) rtol = 1e-12
    end

    @testset "an arc polygon bounds a rounded rectangle" begin
        w, h, ρ = 4.0e-3, 2.0e-3, 0.5e-3
        shape = ArcPolygon(Curve{Float64}[
        Segment((-w / 2, -h / 2 - ρ), (w / 2, -h / 2 - ρ)),
        Arc(Placement(w / 2, -h / 2, -π / 2), ρ, π / 2),
        Segment((w / 2 + ρ, -h / 2), (w / 2 + ρ, h / 2)),
        Arc(Placement(w / 2, h / 2, 0), ρ, π / 2),
        Segment((w / 2, h / 2 + ρ), (-w / 2, h / 2 + ρ)),
        Arc(Placement(-w / 2, h / 2, π / 2), ρ, π / 2),
        Segment((-w / 2 - ρ, h / 2), (-w / 2 - ρ, -h / 2)),
        Arc(Placement(-w / 2, -h / 2, π), ρ, π / 2)])
        @test area(shape) ≈ w * h + 2ρ * (w + h) + π * ρ^2 rtol = 1e-12
        @test perimeter(shape) ≈ 2(w + h) + 2π * ρ rtol = 1e-12
        @test hypot(centroid(shape)...) <= 1e-12 * hypot(w, h)
        @test distance(shape) == 0
        placed = LineCableModels.resolve(Placement(0.01, 0.0), shape)
        @test distance(placed) ≈ 0.01 - w / 2 - ρ rtol = 1e-12
        @test support(placed) ≈ hypot(0.01 + w / 2, h / 2) + ρ rtol = 1e-12
    end

    @testset "rectangles contain regions and shapes have directional supports" begin
        w, h = 4.0e-3, 2.0e-3
        rectangle = Rectangle(w, h)
        @test support(rectangle) ≈ hypot(w, h) / 2 rtol = 1e-12
        @test distance(rectangle) == 0
        @test perimeter(rectangle) ≈ 2(w + h) rtol = 1e-12
        @test centroid(rectangle) == (0, 0)
        at = Placement(0.01, 0.005)
        placed = LineCableModels.resolve(at, rectangle)
        @test support(placed) ≈ hypot(at.x, at.y) + hypot(w, h) / 2 rtol = 1e-12
        @test distance(placed) ≈ hypot(0.008, 0.004) rtol = 1e-12
        @test perimeter(placed) ≈ 2(w + h) rtol = 1e-12
        @test centroid(placed) == (at.x, at.y)
        circle = Circle(1.0e-3, at)
        annulus = Annulus(0.5e-3, 1.0e-3, at)
        for φ in (0, π / 2, π, 3π / 2, 1.0)
            projection = at.x * cos(φ) + at.y * sin(φ)
            @test support(placed, φ) ≈
                  projection + w * abs(cos(φ)) / 2 +
                  h * abs(sin(φ)) / 2 rtol = 1e-12
            @test support(circle, φ) ≈ projection + circle.r rtol = 1e-12
            @test support(annulus, φ) ≈ projection + annulus.ro rtol = 1e-12
        end
    end
end
