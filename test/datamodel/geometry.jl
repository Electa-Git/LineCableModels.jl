using LineCableModels.DataModel: Sector, Ellipse, Polygon, ArcPolygon
using LineCableModels.BaseGeometry: Curve, Segment, Arc, EllipticArc, Contour, parallel,
                                    istangent,
                                    clearance,
                                    dimension, turn, lowest, vertex, apsides, winding,
                                    moment

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

    @testset "a circular layer has the outline of the generic offset" begin
        place = Placement(0.013, -0.7, 0.4)
        probed = (Circle(2.0e-3), Circle(2.0e-3, place), Annulus(1.0e-3, 3.0e-3, place))
        for inner in probed, t in (1.0e-4, 7.0e-4, 2.5e-3)

            layer = LineCableModels.resolve(inner, Offset(t))
            general = invoke(LineCableModels.resolve,
                Tuple{LineCableModels.AbstractShape, Offset}, inner, Offset(t))
            @test layer isa Annulus
            a, b = Contour(layer), Contour(general)
            @test length.([a.outer, a.holes...]) == length.([b.outer, b.holes...]) == [4, 4]
            curves = [(c, d)
                      for (p, q) in zip([a.outer, a.holes...], [b.outer, b.holes...])
                      for (c, d) in zip(p, q)]
            @test all((c.at.x, c.at.y, c.r) === (d.at.x, d.at.y, d.r) for (c, d) in curves)
            # The parallel takes its angles from `atan`. Its fourth arc starts at -π/2, and
            # the circle's fourth arc at 3π/2. At the origin, the parallel's first arc starts
            # at -2.4e-16 rad, and its span is one ulp above π/2.
            @test all(isapprox(rem2pi(c.at.φ - d.at.φ, RoundNearest), 0; atol = eps(2π)) &&
                      isapprox(c.span, d.span; atol = eps(2π)) for (c, d) in curves)
            # The general layer's area comes from Green's integrals, the annulus's from its
            # closed form. They differ by 1.6e-15 relative at most, at 0.1 mm on the placed
            # circle.
            @test area(layer) ≈ area(general) rtol = 4.0e-15
            @test dimension(layer) === dimension(general)
        end
        # Every typed offset method other than the generic one takes a probed type first.
        typed = [Base.unwrap_unionall(m.sig).parameters[2]
                 for m in methods(LineCableModels.resolve)
                 if length(Base.unwrap_unionall(m.sig).parameters) == 3 &&
                    Base.unwrap_unionall(m.sig).parameters[3] <: Offset &&
                    Base.unwrap_unionall(m.sig).parameters[2] !=
                    LineCableModels.AbstractShape]
        @test typed ==
              [Union{(Base.typename(typeof(inner)).wrapper for inner in probed)...}]
    end

    @testset "the generic offset rests on the outer boundary" begin
        layer = invoke(
            LineCableModels.resolve, Tuple{LineCableModels.AbstractShape, Offset},
            Annulus(1.0e-3, 3.0e-3), Offset(2.5e-3))
        @test area(layer) ≈ π * ((3.0e-3 + 2.5e-3)^2 - (3.0e-3)^2) rtol = 1.0e-12
        bar = LineCableModels.PerforatedShape(Rectangle(12.0e-3, 6.0e-3),
            [Circle(0.5e-3, Placement(2.0e-3, 0))])
        for t in (0.4e-3, 1.0e-3)
            layer = LineCableModels.resolve(bar, Offset(t))
            @test area(layer) ≈ 36.0e-3 * t + π * t^2 rtol = 1.0e-12
        end
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

    @testset "a circular pattern of radius zero turns a sector about its center" begin
        segment = Sector(0, 31.15e-3, π / 3, 0, 1.25e-3, Placement(0, 0))
        course = placements(CircularPattern(6; r = 0.0), segment)
        @test course.shape === segment
        @test all(placement -> placement.x == 0 && placement.y == 0, course.placements)
        @test [placement.φ for placement in course.placements] ≈ (0:5) .* (π / 3)
        @test_throws DomainError placements(CircularPattern(6; r = 1.0e-3), segment)
        @test_throws DomainError placements(CircularPattern(7; r = 0.0), segment)
        wide = Sector(0, 31.15e-3, π / 3 + 1.0e-9, 0, 0, Placement(0, 0))
        @test length(placements(CircularPattern(6; r = 0.0), wide).placements) == 6
    end

    @testset "a placement locates a cross-section in its frame" begin
        wire = LineCableModels.resolve(Placement(0, -1), Circle(1.5e-3, Placement(3.0e-3, 0)))
        @test all(centroid(wire) .≈ (3.0e-3, -1.0))
        @test support(wire) ≈ hypot(3.0e-3, 1.0) + 1.5e-3
    end

    @testset "the inverse of a placement undoes it" begin
        for p in (Placement(1.0e-3, -2.0e-3, 0.3), Placement(-0.25, 0.5, 2.0),
            Placement(0.04, 0.0, -π / 3))
            for q in (inv(p) * p, p * inv(p))
                @test [q.x, q.y, q.φ] ≈ [0, 0, 0] atol = 1e-15
            end
        end
    end

    @testset "the nearest distance of a shape from the origin" begin
        @test distance(Circle(1.0e-3, Placement(3.0e-3, 0))) ≈ 2.0e-3
        @test distance(Circle(1.0e-3)) == 0
        @test distance(Annulus(2.0e-3, 3.0e-3)) == 2.0e-3
        sector = placements(CircularPattern(18; r = 9.0e-3), Rectangle(3.0e-3, 3.0e-3)).shape
        @test distance(sector) ≈ 7.5e-3
        @test distance(Ellipse(3e-3, 2e-3, Placement(0, 10e-3))) ≈ 8e-3 rtol = 1e-12
        @test distance(Ellipse(3e-3, 2e-3, Placement(0, 10e-3, π / 2))) ≈ 7e-3 rtol = 1e-12
    end

    @testset "the distance of a sector is its contour's nearest point" begin
        at = Placement(0, 0)
        for shape in (Sector(8.5e-3, 11.5e-3, 0.3, 0, 0, at),
                Sector(0, 31.15e-3, π / 3, 0, 1.25e-3, at),
                Sector(10.0e-3, 10.2e-3, π, 0, 0, at),
                Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, at),
                Sector(0, 10.0e-3, 2π / 3, 1.0e-3, 0, at)),
            place in (Placement(0, 10.0), Placement(0.02, -0.01, π / 5))

            placed = LineCableModels.resolve(place, shape)
            curves, parameters = Contour(placed).outer, range(0, 1; length = 20_001)
            sample = minimum(hypot(curve(s)...) for curve in curves, s in parameters)
            @test distance(placed) <= sample
            @test distance(placed) ≈ sample rtol = 1e-10
        end
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

    # The origin lies in a hole centered on it, in a hole off it, in the region, or far
    # outside the cross-section. Far outside, the distance is that of the outer cross-section.
    @testset "the distance of a cross-section with holes is the distance of its contour" begin
        far = Circle(3.0e-3, Placement(10.0e-3, 0))
        for (outer, hole, nearest) in (
            (Circle(12.0e-3), Circle(8.0e-3), 8.0e-3),
            (Circle(12.0e-3), Circle(3.0e-3, Placement(1.0e-3, 0)), 2.0e-3),
            (Circle(12.0e-3), Circle(3.0e-3, Placement(6.0e-3, 0)), 0.0),
            (far, Circle(1.0e-3, Placement(10.0e-3, 0)), distance(far)))
            shape = LineCableModels.PerforatedShape(outer, [hole])
            @test distance(shape) ≈ nearest rtol = 1e-12
            @test distance(shape) == distance(Contour(shape))
        end
    end

    @testset "disjoint shapes have the area-weighted centroid of their centroids" begin
        shapes = [Circle(1.0e-3), Circle(2.0e-3, Placement(3.0e-3, -3.0e-3))]
        @test all(centroid(shapes) .≈ (2.4e-3, -2.4e-3))
        @test all(centroid(shapes[2:2]) .≈ centroid(shapes[2]))
    end

    # The centroid takes the first moments about the centroid of `first(shapes)`. A circle off the
    # origin keeps its centroid to the bit, and at a measured radius a zero spread.
    @testset "the centroid of one shape is the centroid of that shape" begin
        at = Placement(-28.0e-3, 16.0e-3)
        @test centroid([Circle(19.5e-3, at)]) === centroid(Circle(19.5e-3, at))
        r = measurement(19.5e-3, 1.0e-4)
        x, y = centroid([Circle(r, at)])
        @test (x, y) == centroid(Circle(r, at))
        @test Measurements.uncertainty(x) == 0 && Measurements.uncertainty(y) == 0
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
            t = sector.apex * sin(sector.span / 2) + 1.0e-3
            layer = LineCableModels.resolve(sector, Offset(t))
            @test area(layer) ≈ perimeter(sector) * t + π * t^2 rtol = 1e-12
            @test perimeter(boundary(layer)) ≈ perimeter(sector) + 2π * t rtol = 1e-12
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
        # Off the origin, each is nearest at a corner: the end of the inner arc, and the apex.
        placed = Placement(0.02, -0.01)
        corner = (S1.ri * cos(S1.span / 2), S1.ri * sin(S1.span / 2))
        @test distance(LineCableModels.resolve(placed, S1)) ≈
              hypot(placed.x + corner[1], placed.y + corner[2]) rtol = 1e-12
        @test distance(LineCableModels.resolve(placed, S5)) ≈
              hypot(placed.x + S5.apex, placed.y) rtol = 1e-12
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
        layer = LineCableModels.resolve(S3, Offset(t))
        @test layer.outer isa ArcPolygon
        @test area(layer) ≈ perimeter(S3) * t + π * t^2 rtol = 1e-12
        @test perimeter(boundary(layer)) ≈ perimeter(S3) + 2π * t rtol = 1e-12
    end

    @testset "a layer around a sector rests on its contour's parallel" begin
        at = Placement(0, 0)
        for (shape, t) in ((Sector(8.5e-3, 11.5e-3, 0.3, 0, 0, at), 0.5e-3),
            (Sector(0, 31.15e-3, π / 3, 0, 1.25e-3, at), 0.1e-3),
            (Sector(10.0e-3, 10.2e-3, π, 0, 0, at), 0.5e-3),
            (Sector(10.0e-3, 10.2e-3, 2π - 0.1, 0, 0, at), 0.3e-3))
            layer = LineCableModels.resolve(shape, Offset(t))
            @test layer.outer isa ArcPolygon
            @test area(layer) ≈ perimeter(shape) * t + π * t^2 rtol = 1e-12
            @test perimeter(boundary(layer)) ≈ perimeter(shape) + 2π * t rtol = 1e-12
            @test @inferred(parallel(Contour(shape), t)) isa Contour{Float64}
        end
        strand = Sector(8.5e-3, 11.5e-3, 0.3, 0, 0, at)
        @test LineCableModels.resolve(strand, Offset(0.99 * 8.5e-3)) isa
              LineCableModels.PerforatedShape
        @test_throws DomainError LineCableModels.resolve(strand, Offset(8.5e-3))
        @test_throws "the outward parallel at distance 0.0085 m self-intersects at curve 3 of path 1" LineCableModels.resolve(
            strand, Offset(8.5e-3))
        tape = Sector(10.0e-3, 10.2e-3, 2π - 0.1, 0, 0, at)
        @test_throws DomainError LineCableModels.resolve(tape, Offset(1.0e-3))
        @test_throws "self-intersects at curve 1 of path 1 and curve 4 of path 1" LineCableModels.resolve(
            tape, Offset(1.0e-3))
        face = EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, 0.3, 2.0, 0.0)
        @test @inferred(turn(face, 0.5e-3, face(0), face(1))) ==
              turn(parallel(face, 0.5e-3))
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
    @testset "ellipse moments and directional support" begin
        ellipse = Ellipse(3.0e-3, 2.0e-3, Placement(0.01, -0.02, π / 7))
        @test_throws DomainError Ellipse(0, 2.0e-3)
        @test_throws DomainError Ellipse(3.0e-3, -1.0e-3)
        @test_throws DomainError Ellipse(Inf, 2.0e-3)
        @test area(ellipse) ≈ 6π * 1.0e-6 rtol = 1e-12
        @test collect(centroid(ellipse)) ≈ [0.01, -0.02] rtol = 1e-12
        @test perimeter(ellipse) ≈ 0.012 * 1.3221199657742158 rtol = 1e-12
        @test perimeter(Ellipse(2.0e-3, 3.0e-3)) ≈ perimeter(ellipse) rtol = 1e-12
        for φ in (0, π / 2, π, 3π / 2, 1.0)
            @test support(ellipse, φ) ≈
                  0.01cos(φ) - 0.02sin(φ) +
                  sqrt(9.0e-6 * cos(φ - π / 7)^2 + 4.0e-6 * sin(φ - π / 7)^2) rtol = 1e-12
        end
        contour = LineCableModels.BaseGeometry.Contour(ellipse)
        @test area(contour) ≈ area(ellipse) rtol = 1e-12
        @test collect(centroid(contour)) ≈ collect(centroid(ellipse)) rtol = 1e-12
        @test perimeter(contour) ≈ perimeter(ellipse) rtol = 1e-12
    end

    @testset "an elliptic arc bounds a half ellipse" begin
        shape = ArcPolygon(Curve{Float64}[
        EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, 0, π),
        Segment((-3.0e-3, 0.0), (3.0e-3, 0.0))])
        @test area(shape) ≈ 3π * 1.0e-6 rtol = 1e-12
        @test collect(centroid(shape)) ≈ [0, 8.0e-3 / (3π)] rtol = 1e-12
        @test perimeter(shape) ≈ 6.0e-3 + 0.006 * 1.3221199657742158 rtol = 1e-12
    end

    @testset "polygon moments and directional support" begin
        @test_throws DomainError Polygon([(0, 0), (1, 0)])
        @test_throws DomainError Polygon([(0, 0), (0, 1), (1, 0)])
        @test_throws DomainError Polygon([(0, 0), (3, 2), (0, 2), (2, 0)])
        pentagon = [
            (0.0, 0.0), (0.004, 0.0), (0.005, 0.003), (0.002, 0.005), (-0.001, 0.003)]
        L = [(0.0, 0.0), (0.004, 0.0), (0.004, 0.001),
            (0.001, 0.001), (0.001, 0.003), (0.0, 0.003)]
        for (vertices, A, c, P, d, projections) in (
            (pentagon, 21.0e-6, [0.002, 46.0e-3 / 21],
                (4 + 2sqrt(10) + 2sqrt(13)) * 1.0e-3, sqrt(370) * 1.0e-3,
                (0.005, 0.005, 0.001, 0.0, 0.002cos(1) + 0.005sin(1))),
            (L, 6.0e-6, [0.0015, 0.001], 0.014, 0.02,
                (0.004, 0.003, 0.003sin(π), 0.0, 0.001cos(1) + 0.003sin(1))))
            polygon = Polygon(vertices)
            @test area(polygon) ≈ A rtol = 1e-12
            @test collect(centroid(polygon)) ≈ c rtol = 1e-12
            @test perimeter(polygon) ≈ P rtol = 1e-12
            @test distance(Polygon(vertices, Placement(-0.0005, -0.0005))) == 0
            @test distance(Polygon(vertices, Placement(0.02, 0.0))) ≈ d rtol = 1e-12
            for (φ, projection) in zip((0, π / 2, π, 3π / 2, 1.0), projections)
                @test support(polygon, φ) ≈ projection rtol = 1e-12
            end
        end
    end

    @testset "the polygon's constructors infer on vertices of one type" begin
        vertices = [(0.0, 0.0), (4.0e-3, 0.0), (4.0e-3, 3.0e-3), (0.0, 3.0e-3)]
        at = Placement(1.0e-3, -2.0e-3, 0.3)
        polygon = Polygon(vertices)
        @test @inferred(Polygon(vertices, at)) isa Polygon{Float64}
        @test @inferred(LineCableModels.resolve(at, polygon)) isa Polygon{Float64}
        @test @inferred(polygon ⊆ Circle(10.0e-3))
        # Literals that mix integers and floats build the same polygon.
        for mixed in ([(0, 0.0), (4.0e-3, 0), (4.0e-3, 3.0e-3), (0, 3.0e-3)],
            [(0, 0.0), (4.0e-3, 0.0), (4.0e-3, 3.0e-3), (0, 3.0e-3)])
            @test Polygon(mixed) isa Polygon{Float64} && Polygon(mixed).vertices == vertices
        end
    end

    @testset "contour offsets satisfy Steiner's formulas" begin
        pentagon = Polygon([
            (0.0, 0.0), (0.004, 0.0), (0.005, 0.003), (0.002, 0.005), (-0.001, 0.003)])
        for (shape, A, P) in ((Rectangle(4.0e-3, 2.0e-3), 8.0e-6, 0.012),
            (pentagon, 21.0e-6, (4 + 2sqrt(10) + 2sqrt(13)) * 1.0e-3))
            layer = LineCableModels.resolve(shape, Offset(0.5e-3))
            @test layer isa LineCableModels.PerforatedShape
            @test area(layer) ≈ P * 0.5e-3 + π * 0.5e-3^2 rtol = 1e-12
            @test perimeter(boundary(layer)) ≈ P + 2π * 0.5e-3 rtol = 1e-12
            @test area(boundary(layer)) ≈ A + P * 0.5e-3 + π * 0.5e-3^2 rtol = 1e-12
            second = LineCableModels.resolve(boundary(layer), Offset(0.25e-3))
            @test area(boundary(second)) ≈ A + P * 0.75e-3 + π * 0.75e-3^2 rtol = 1e-12
        end
        # Each of the five convex right corners of the L adds π t² / 4, and its reflex corner
        # removes t² tan(π / 4).
        L = Polygon([(0.0, 0.0), (0.004, 0.0), (0.004, 0.001),
            (0.001, 0.001), (0.001, 0.003), (0.0, 0.003)])
        t = 0.5e-3
        layer = LineCableModels.resolve(L, Offset(t))
        @test area(boundary(layer)) ≈ 6.0e-6 + 0.014 * t + (5π / 4 - 1) * t^2 rtol = 1e-12
        # The rectangle less the upper half of an ellipse has cusps at its two lower corners.
        # Its parallel is simple, and a count of 400 000 random points gives the area of the
        # dilation as 1.932462e-5 ± 2.0e-8 m².
        clockwise = ArcPolygon(Curve{Float64}[
        Segment((3.0e-3, 0.0), (3.0e-3, 3.0e-3)),
        Segment((3.0e-3, 3.0e-3), (-3.0e-3, 3.0e-3)),
        Segment((-3.0e-3, 3.0e-3), (-3.0e-3, 0.0)),
        EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, π, -π)])
        layer = LineCableModels.resolve(clockwise, Offset(t))
        P = 12.0e-3 + perimeter(Ellipse(3.0e-3, 2.0e-3)) / 2
        @test area(boundary(layer)) ≈ (18 - 3π) * 1.0e-6 + P * t + π * t^2 rtol = 1e-12
        @test abs(area(boundary(layer)) - 1.932462e-5) <= 3 * 2.0e-8
    end
    @testset "an elliptic arc carries its normal offset" begin
        t = 0.5e-3
        @test_throws DomainError EllipticArc(
            Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2, -1.4e-3)
        @test_throws DomainError EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2, Inf)
        arc = EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2, t)
        for (s, point) in ((0, [3.5e-3, 0]),
            (1 / 2, [3.0e-3 / √2 + t * 2 / √13, 2.0e-3 / √2 + t * 3 / √13]),
            (1, [0, 2.5e-3]))
            @test collect(arc(s)) ≈ point rtol = 1e-12
        end
    end

    @testset "the parallel of an ellipse satisfies Steiner's formulas" begin
        t = 0.5e-3
        at = Placement(0.01, -0.02, π / 7)
        for (a, b, P) in ((3.0e-3, 2.0e-3, 0.012 * 1.3221199657742158),
            (2.0e-3, 3.0e-3, 0.012 * 1.3221199657742158),
            (2.0e-3, 2.0e-3, 4π * 1.0e-3))
            layer = LineCableModels.resolve(Ellipse(a, b, at), Offset(t))
            @test layer isa LineCableModels.PerforatedShape
            parallel = boundary(layer)
            @test parallel isa ArcPolygon
            @test length(parallel.outer) == 4
            @test all(curve -> curve isa EllipticArc && curve.offset == t, parallel.outer)
            @test area(parallel) ≈ π * a * b + P * t + π * t^2 rtol = 1e-12
            @test perimeter(parallel) ≈ P + 2π * t rtol = 1e-12
            @test collect(centroid(parallel)) ≈ [0.01, -0.02] rtol = 1e-12
            @test area(layer) ≈ P * t + π * t^2 rtol = 1e-12
            points = Tuple{Float64, Float64}[]
            for k in 0:4095
                θ = 2π * k / 4096
                ρ = hypot(b * cos(θ), a * sin(θ))
                x = a * cos(θ) + t * b * cos(θ) / ρ
                y = b * sin(θ) + t * a * sin(θ) / ρ
                push!(points,
                    (0.01 + x * cos(π / 7) - y * sin(π / 7),
                        -0.02 + x * sin(π / 7) + y * cos(π / 7)))
            end
            tessellation = Polygon(points)
            @test area(tessellation) ≈ π * a * b + P * t + π * t^2 rtol = 1e-6
            @test perimeter(tessellation) ≈ P + 2π * t rtol = 1e-6
            @test collect(centroid(tessellation)) ≈ [0.01, -0.02] rtol = 1e-6
            for φ in (0, π / 2, π, 3π / 2, 1.0, 2.5)
                @test support(parallel, φ) ≈
                      0.01cos(φ) - 0.02sin(φ) +
                      sqrt(a^2 * cos(φ - π / 7)^2 + b^2 * sin(φ - π / 7)^2) + t rtol = 1e-12
            end
        end
    end

    @testset "the parallel of a half ellipse" begin
        t = 0.5e-3
        for (a, b, P) in ((3.0e-3, 2.0e-3, 0.012 * 1.3221199657742158),
            (2.0e-3, 3.0e-3, 0.012 * 1.3221199657742158),
            (2.0e-3, 2.0e-3, 4π * 1.0e-3))
            half = ArcPolygon(Curve{Float64}[
            EllipticArc(Placement(0, 0), a, b, 0, π), Segment((-a, 0.0), (a, 0.0))])
            parallel = boundary(LineCableModels.resolve(half, Offset(t)))
            @test area(parallel) ≈ π * a * b / 2 + (2a + P / 2) * t + π * t^2 rtol = 1e-12
            @test perimeter(parallel) ≈ 2a + P / 2 + 2π * t rtol = 1e-12
            points = Tuple{Float64, Float64}[]
            for k in 0:2048
                θ = π * k / 2048
                ρ = hypot(b * cos(θ), a * sin(θ))
                push!(points, (a * cos(θ) + t * b * cos(θ) / ρ,
                    b * sin(θ) + t * a * sin(θ) / ρ))
            end
            for k in 1:1024
                θ = π + π / 2 * k / 1024
                push!(points, (-a + t * cos(θ), t * sin(θ)))
            end
            for k in 0:1022
                θ = 3π / 2 + π / 2 * k / 1023
                push!(points, (a + t * cos(θ), t * sin(θ)))
            end
            tessellation = Polygon(points)
            @test area(tessellation) ≈ area(parallel) rtol = 1e-6
            @test collect(centroid(tessellation)) ≈ collect(centroid(parallel)) rtol = 1e-6
        end
    end

    @testset "a second offset adds the distances" begin
        E = Ellipse(3.0e-3, 2.0e-3, Placement(0.01, -0.02, π / 7))
        second = boundary(LineCableModels.resolve(
            boundary(LineCableModels.resolve(E, Offset(0.3e-3))), Offset(0.2e-3)))
        single = boundary(LineCableModels.resolve(E, Offset(0.5e-3)))
        @test area(second) ≈ area(single) rtol = 1e-15
        @test [curve.offset for curve in second.outer] ==
              [curve.offset for curve in single.outer]
    end

    @testset "a layer contains what lies in its band" begin
        band = LineCableModels.resolve(Rectangle(4.0e-3, 2.0e-3), Offset(1.0e-3))
        @test Circle(0.4e-3, Placement(0, 1.5e-3)) ⊆ band
        @test !(Circle(0.4e-3, Placement(0, 1.2e-3)) ⊆ band)
        @test !(Circle(0.4e-3, Placement(0, 1.8e-3)) ⊆ band)
        band = LineCableModels.PerforatedShape(Rectangle(10.0e-3, 10.0e-3),
            AbstractShape{Float64}[Circle(2.0e-3)])
        @test Circle(1.0e-3, Placement(3.5e-3, 0)) ⊆ band
        @test !(Circle(1.0e-3, Placement(2.5e-3, 0)) ⊆ band)
    end

    # Each distance of the curve geometry runs along the outward normal: positive out of the
    # region, negative into it. A path runs counterclockwise about the region it bounds.
    @testset "a curve's direction and turn" begin
        segment = Segment((0.0, 0.0), (3.0e-3, 4.0e-3))
        arc = Arc(Placement(0, 0, π / 6), 2.0e-3, π / 2)
        reversed = Arc(Placement(0, 0, 0), 2.0e-3, -π / 2)
        elliptic = EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2, 0.5e-3)
        for t in (0, 1 / 2, 1)
            @test LineCableModels.BaseGeometry.direction(segment, t) ≈ atan(4, 3) rtol = 1e-12
            @test LineCableModels.BaseGeometry.direction(arc, t) ≈ 2π / 3 + t * π / 2 rtol = 1e-12
            @test LineCableModels.BaseGeometry.direction(reversed, t) ≈ -π / 2 - t * π / 2 rtol = 1e-12
        end
        @test LineCableModels.BaseGeometry.direction(elliptic, 0) ≈ π / 2 rtol = 1e-12
        @test LineCableModels.BaseGeometry.direction(elliptic, 1) ≈ π rtol = 1e-12
        @test turn(segment).fraction == 0
        @test turn(segment).length ≈ 5.0e-3 rtol = 1e-12
        @test turn(arc).fraction ≈ 1 / 4 rtol = 1e-12
        @test turn(arc).length ≈ π * 1.0e-3 rtol = 1e-12
        @test turn(elliptic).fraction ≈ 1 / 4 rtol = 1e-12
        @test turn(elliptic).length ≈ 0.003 * 1.3221199657742158 + 0.5e-3 * π / 2 rtol = 1e-12
    end

    @testset "the parallel of a curve" begin
        moved = parallel(Segment((0.0, 0.0), (4.0e-3, 0.0)), 0.5e-3)
        @test collect(moved.p) ≈ [0, -0.5e-3] atol = 1e-15
        @test collect(moved.q) ≈ [4.0e-3, -0.5e-3] atol = 1e-15
        @test parallel(Arc(Placement(0, 0, 0), 2.0e-3, π / 2), 0.5e-3).r ≈ 2.5e-3 rtol = 1e-12
        @test parallel(Arc(Placement(0, 0, 0), 2.0e-3, π / 2), -0.5e-3).r ≈ 1.5e-3 rtol = 1e-12
        @test parallel(Arc(Placement(0, 0, 0), 2.0e-3, -π / 2), 0.5e-3).r ≈ 1.5e-3 rtol = 1e-12
        elliptic = EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2, 0.5e-3)
        @test parallel(elliptic, -0.2e-3).offset ≈ 0.3e-3 rtol = 1e-12
        # The smallest radius of curvature of this quarter is b² / a = 4/3 mm, at θ = 0.
        @test EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2, -1.3e-3).offset ==
              -1.3e-3
        @test_throws DomainError EllipticArc(
            Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2, -1.4e-3)
    end

    @testset "the inward parallel of an ellipse satisfies Steiner's formulas" begin
        a, b, d = 3.0e-3, 2.0e-3, -1.2e-3
        P = 0.012 * 1.3221199657742158
        inner = ArcPolygon(parallel(Contour(Ellipse(a, b, Placement(0.01, -0.02, π / 7))), d).outer)
        @test all(curve -> curve isa EllipticArc && curve.offset == d, inner.outer)
        @test area(inner) ≈ π * a * b + P * d + π * d^2 rtol = 1e-12
        @test perimeter(inner) ≈ P + 2π * d rtol = 1e-12
        @test collect(centroid(inner)) ≈ [0.01, -0.02] rtol = 1e-12
        points = Tuple{Float64, Float64}[]
        for k in 0:4095
            θ = 2π * k / 4096
            ρ = hypot(b * cos(θ), a * sin(θ))
            x = a * cos(θ) + d * b * cos(θ) / ρ
            y = b * sin(θ) + d * a * sin(θ) / ρ
            push!(points,
                (0.01 + x * cos(π / 7) - y * sin(π / 7),
                    -0.02 + x * sin(π / 7) + y * cos(π / 7)))
        end
        tessellation = Polygon(points)
        @test area(tessellation) ≈ area(inner) rtol = 1e-6
        @test perimeter(tessellation) ≈ perimeter(inner) rtol = 1e-6
        @test collect(centroid(tessellation)) ≈ collect(centroid(inner)) rtol = 1e-6
    end

    # The parallels of two consecutive curves meet at one point at a converging corner, a convex
    # corner inward or a concave one outward. They also meet at a tangent junction, where the
    # shared point moves along the shared normal. Where they diverge, each end moves along its
    # own normal.
    @testset "the junctions of two parallels" begin
        at = Placement(0, 0)
        @testset "a sharp quarter sector at the depth $e m" for e in
                                                                (1.0e-3, 2.0e-3, 3.0e-3)
            r = 10.0e-3
            outer, inward, outward = Contour(Sector(0, r, π / 2, 0, 0, at)).outer
            s = sqrt((r - e)^2 - e^2)
            p1 = (s * cos(π / 4) + e * cos(π / 4), -s * sin(π / 4) + e * sin(π / 4))
            p2 = (s * cos(π / 4) + e * cos(π / 4), s * sin(π / 4) - e * sin(π / 4))
            p3 = (e * sqrt(2), 0.0)
            for (before, after, point) in ((outward, outer, p1), (outer, inward, p2),
                (inward, outward, p3))
                junction = parallel(before, after, -e)
                @test collect(junction[1]) ≈ collect(point) rtol = 1e-12
                @test junction[1] == junction[2]
                p = vertex(parallel(before, -e), parallel(after, -e), before(1))
                @test collect(p) ≈ collect(point) rtol = 1e-12
            end
            @test turn(outward, -e, p3, p1).length ≈ s - e rtol = 1e-12
            @test turn(outer, -e, p1, p2).length ≈ (r - e) * (π / 2 - 2atan(e / s)) rtol = 1e-12
        end
        @testset "a bent strand at the depth $e m" for e in (0.5e-3, 1.0e-3, 1.4e-3)
            ri, ro = 8.5e-3, 11.5e-3
            outer, inward, inner, outward = Contour(Sector(ri, ro, 0.3, 0, 0, at)).outer
            so, si = sqrt((ro - e)^2 - e^2), sqrt((ri + e)^2 - e^2)
            p1 = (so * cos(0.15) + e * sin(0.15), -so * sin(0.15) + e * cos(0.15))
            p2 = (so * cos(0.15) + e * sin(0.15), so * sin(0.15) - e * cos(0.15))
            p3 = (si * cos(0.15) + e * sin(0.15), si * sin(0.15) - e * cos(0.15))
            p4 = (si * cos(0.15) + e * sin(0.15), -si * sin(0.15) + e * cos(0.15))
            for (before, after, point) in ((outward, outer, p1), (outer, inward, p2),
                (inward, inner, p3), (inner, outward, p4))
                junction = parallel(before, after, -e)
                @test collect(junction[1]) ≈ collect(point) rtol = 1e-12
                @test junction[1] == junction[2]
                p = vertex(parallel(before, -e), parallel(after, -e), before(1))
                @test collect(p) ≈ collect(point) rtol = 1e-12
            end
            @test turn(outward, -e, p4, p1).length ≈ so - si rtol = 1e-12
            @test turn(inward, -e, p2, p3).length ≈ so - si rtol = 1e-12
            @test turn(outer, -e, p1, p2).length ≈ (ro - e) * (0.3 - 2asin(e / (ro - e))) rtol = 1e-12
            @test turn(inner, -e, p3, p4).length ≈ (ri + e) * (0.3 - 2asin(e / (ri + e))) rtol = 1e-12
        end
        @testset "a Milliken segment at the depth $e m" for e in (0.3e-3, 0.6e-3, 0.9e-3)
            R, ρ = 31.15e-3, 1.25e-3
            outer, back, side, tip = Contour(Sector(0, R, π / 3, 0, ρ, at)).outer[1:4]
            a, b = side(0), side(1)
            p1 = outer(1) .* ((R - e) / R)
            p2 = a .+ (e / ρ) .* ((back.at.x, back.at.y) .- a)
            p3 = b .+ (e / ρ) .* ((tip.at.x, tip.at.y) .- b)
            for (junction, point) in ((parallel(outer, back, -e), p1),
                (parallel(back, side, -e), p2), (parallel(side, tip, -e), p3))
                @test collect(junction[1]) ≈ collect(point) rtol = 1e-12
                @test junction[1] == junction[2]
            end
            @test turn(side, -e, p2, p3).length ≈ hypot((b .- a)...) rtol = 1e-12
            @test turn(back, -e, p1, p2).length ≈ (ρ - e) * abs(back.span) rtol = 1e-12
        end
        @testset "a reflex corner inward and a convex corner outward diverge" begin
            e = 1.0e-3
            outer, inward, outward = Contour(Sector(0, 10.0e-3, 4π / 3, 0, 0, at)).outer
            p, q = parallel(inward, outward, -e)
            @test collect(p) ≈ [e * sqrt(3) / 2, e / 2] rtol = 1e-12
            @test collect(q) ≈ [e * sqrt(3) / 2, -e / 2] rtol = 1e-12
            bottom, right = Contour(Rectangle(4.0e-3, 2.0e-3)).outer[1:2]
            p, q = parallel(bottom, right, 0.5e-3)
            @test collect(p) ≈ [2.0e-3, -1.5e-3] rtol = 1e-12
            @test collect(q) ≈ [2.5e-3, -1.0e-3] rtol = 1e-12
        end
        @testset "two arcs at a corner" begin
            e = 0.5e-3
            a = Arc(Placement(-1.0e-3, 0, -π / 3), 2.0e-3, 2π / 3)
            b = Arc(Placement(1.0e-3, 0, 2π / 3), 2.0e-3, 2π / 3)
            p = vertex(parallel(a, -e), parallel(b, -e), a(1))
            @test collect(p) ≈ [0, sqrt(1.25) * 1.0e-3] atol = 1e-15
            @test parallel(a, b, -e) == (p, p)
        end
        @testset "an elliptic arc at a corner" begin
            arc = EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, 0, π)
            chord = Segment((-3.0e-3, 0.0), (3.0e-3, 0.0))
            p, q = parallel(arc, chord, 0.5e-3)
            @test collect(p) ≈ [-3.5e-3, 0] atol = 1e-15
            @test collect(q) ≈ [-3.0e-3, -0.5e-3] atol = 1e-15
            @test_throws DomainError parallel(arc, chord, -0.5e-3)
            @test_throws DomainError vertex(
                parallel(arc, -0.5e-3), parallel(chord, -0.5e-3),
                arc(1))
        end
    end

    @testset "the parallel of a contour" begin
        t = 0.5e-3
        hand = Curve{Float64}[
        Segment((-2.0e-3, -1.5e-3), (2.0e-3, -1.5e-3)),
        Arc(
            Placement(2.0e-3, -1.0e-3, -π / 2), t, π / 2),
        Segment(
            (2.5e-3, -1.0e-3), (2.5e-3, 1.0e-3)),
        Arc(
            Placement(2.0e-3, 1.0e-3, 0), t, π / 2),
        Segment(
            (2.0e-3, 1.5e-3), (-2.0e-3, 1.5e-3)),
        Arc(
            Placement(-2.0e-3, 1.0e-3, π / 2), t, π / 2),
        Segment(
            (-2.5e-3, 1.0e-3), (-2.5e-3, -1.0e-3)),
        Arc(
            Placement(-2.0e-3, -1.0e-3, π), t, π / 2)]
        outline = parallel(Contour(Rectangle(4.0e-3, 2.0e-3)), t).outer
        @test length(outline) == length(hand)
        for (curve, expected) in zip(outline, hand), s in (0, 1 / 2, 1)

            @test collect(curve(s)) ≈ collect(expected(s)) atol = 1e-15
        end
        annulus = parallel(Contour(Annulus(1.0e-3, 3.0e-3, Placement(0, 0))), -t)
        @test all(arc -> arc.r ≈ 2.5e-3, annulus.outer)
        @test all(arc -> arc.r ≈ 1.5e-3, only(annulus.holes))
        reflex = Contour(Sector(0, 10.0e-3, 4π / 3, 0, 0, Placement(0, 0)))
        @test length(parallel(reflex, -1.0e-3).outer) == length(reflex.outer) + 1
    end

    @testset "the inward parallel of a sector wider than π" begin
        r = 10.0e-3
        sector = Sector(0, r, 4π / 3, 0, 0, Placement(0, 0))
        @test dimension(sector) ≈ r rtol = 1e-12
        for e in (1.0e-3, 2.5e-3, 4.5e-3)
            R, s = r - e, sqrt((r - e)^2 - e^2)
            outline = parallel(Contour(sector), -e)
            @test length(outline.outer) == 4
            @test area(outline) ≈ R^2 * (4π / 3 - 2atan(e, s)) / 2 - e * s - π / 6 * e^2 rtol = 1e-12
        end
    end

    @testset "the inward parallel refuses a piece that vanishes" begin
        r = 10.0e-3
        reflex = Contour(Sector(0, r, 4π / 3, 0, 0, Placement(0, 0)))
        @test parallel(reflex, -0.4999r) isa Contour{Float64}
        for e in (0.5001r, 0.51r)
            @test_throws DomainError parallel(reflex, -e)
            @test_throws "has degenerated at curve 2 of path 1: the distance across is $(dimension(reflex)) m" parallel(
                reflex, -e)
        end
        L = Contour(Polygon([(0.0, 0.0), (0.004, 0.0), (0.004, 0.001),
            (0.001, 0.001), (0.001, 0.003), (0.0, 0.003)]))
        @test parallel(L, -0.4e-3) isa Contour{Float64}
        @test_throws DomainError parallel(L, -0.6e-3)
        @test_throws "has degenerated at curve 2 of path 1: the distance across is 0.001 m" parallel(
            L, -0.6e-3)
    end

    @testset "the clearance of two curves" begin
        lower = Segment((0.0, 0.0), (4.0e-3, 0.0))
        upper = Segment((5.0e-3, 1.0e-3), (1.0e-3, 1.0e-3))
        p, q, d, s, t = clearance(lower, upper)
        @test collect(p) ≈ [2.5e-3, 0] atol = 1e-15
        @test collect(q) ≈ [2.5e-3, 1.0e-3] atol = 1e-15
        @test d ≈ 1.0e-3 rtol = 1e-12
        @test collect((s, t)) ≈ [0.625, 0.625] rtol = 1e-12
        segment = Segment((-3.0e-3, -2.0e-3), (3.0e-3, -2.0e-3))
        p, q, d = clearance(segment, Arc(Placement(0, 0, π), 1.0e-3, π))
        @test collect(p) ≈ [0, -2.0e-3] atol = 1e-15
        @test collect(q) ≈ [0, -1.0e-3] atol = 1e-15
        @test d ≈ 1.0e-3 rtol = 1e-12
        inner = Arc(Placement(0, 0, 0), 2.0e-3, π / 2)
        outer = Arc(Placement(0, 0, π / 4), 3.0e-3, π / 2)
        # Every pair at one angle of the common range [π/4, π/2] is nearest, and which one
        # returns depends on the floating-point error. The bounds are two spacings of
        # floating-point numbers at the radii and at the angle.
        p, q, d, s, t = clearance(inner, outer)
        φ = atan(p[2], p[1])
        @test collect(p) ≈ 2.0e-3 .* [cos(φ), sin(φ)] rtol = 2eps()
        @test collect(q) ≈ 3.0e-3 .* [cos(φ), sin(φ)] rtol = 2eps()
        @test d ≈ 3.0e-3 - 2.0e-3 rtol = 2eps()
        @test π / 4 - 2eps(float(π)) <= φ <= π / 2 + 2eps(float(π))
        quarter = EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2)
        moved = EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2, 0.5e-3)
        p, q, d = clearance(quarter, moved)
        @test collect(p) ≈ [3.0e-3 / √2, 2.0e-3 / √2] rtol = 1e-12
        @test collect(q) ≈ [3.0e-3 / √2 + 0.5e-3 * 2 / √13, 2.0e-3 / √2 + 0.5e-3 * 3 / √13] rtol = 1e-12
        @test d === 0.5e-3 - 0.0
        opposite = EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, π, π / 2, 0.5e-3)
        p, q, d, s, t = clearance(quarter, opposite)
        @test (s, t) == (0, 1)
        @test p == quarter(0)
        @test q == opposite(1)
        @test d ≈ hypot(3.0e-3, 2.5e-3) rtol = 1e-12
        another = EllipticArc(Placement(0, 0), 3.0e-3, 2.5e-3, π, π / 2)
        @test_throws DomainError clearance(quarter, another)
    end

    # The distance across the metal is twice the depth at which the inward parallel of the
    # outline first degenerates. The bent strand's inner arc closes between its two sides. A
    # sharp sector's sides close on its arc at its inscribed radius. A fully filleted outline
    # holds a disk of its fillet's radius, and an ellipse twice its smallest radius of curvature.
    @testset "the distance across the metal" begin
        at = Placement(0, 0)
        shapes = (
            ("bar", Rectangle(30.0e-3, 12.0e-3), 12.0e-3),
            ("bent strand", Sector(8.5e-3, 11.5e-3, 0.3, 0, 0, at),
                2 * 8.5e-3 * sin(0.15) / (1 - sin(0.15))),
            ("bent tape", Sector(14.95e-3, 15.05e-3, 5π / 3, 0, 0, at), 0.1e-3),
            ("rounded wedge",
                Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, at), 3.0e-3),
            ("Milliken segment", Sector(0, 31.15e-3, π / 3, 0, 1.25e-3, at), 2.5e-3),
            ("sharp quarter", Sector(0, 10.0e-3, π / 2, 0, 0, at),
                2 * 10.0e-3 * sin(π / 4) / (1 + sin(π / 4))),
            ("sharp third", Sector(0, 10.0e-3, 2π / 3, 0, 0, at),
                2 * 10.0e-3 * sin(π / 3) / (1 + sin(π / 3))),
            ("perforated bar",
                LineCableModels.PerforatedShape(Rectangle(30.0e-3, 12.0e-3),
                    LineCableModels.AbstractShape{Float64}[Circle(2.0e-3)]), 4.0e-3),
            ("ellipse", Ellipse(3.0e-3, 2.0e-3), 2 * 2.0e-3^2 / 3.0e-3),
            ("round ellipse", Ellipse(2.1e-3, 2.0e-3), 2 * 2.0e-3^2 / 2.1e-3),
            ("duct wall",
                LineCableModels.resolve(Ellipse(30.0e-3, 20.0e-3), Offset(3.0e-3)),
                3.0e-3))
        @testset "at y = $y m" for y in (0.0, 10.0)
            place = Placement(0.0, y)
            circle = LineCableModels.resolve(place, Circle(1.5e-3))
            annulus = LineCableModels.resolve(place, Annulus(20.5e-3, 21.5e-3, Placement(0, 0)))
            @test dimension(circle) === 2 * 1.5e-3
            # The distance across a placed annulus differs from the difference of its radii by at
            # most the spacing of floating-point numbers at its coordinate y.
            @test dimension(annulus) ≈ 21.5e-3 - 20.5e-3 atol = eps(y)
            @testset "$id" for (id, shape, across) in shapes
                @test dimension(LineCableModels.resolve(place, shape)) ≈ across rtol = 1e-12
            end
            half = ArcPolygon(Curve{Float64}[
            EllipticArc(Placement(0, y), 3.0e-3, 2.0e-3, 0, π),
            Segment((-3.0e-3, y), (3.0e-3, y))])
            @test_throws DomainError dimension(half)
        end
    end

    @testset "the lowest point of an ellipse" begin
        ellipse = Ellipse(9.0e-3, 3.0e-3)
        for (shape, low) in (
            (ellipse, [0, -3.0e-3]),
            (Ellipse(9.0e-3, 3.0e-3, Placement(0, 0, π / 2)), [0, -9.0e-3]),
            (Ellipse(9.0e-3, 3.0e-3, Placement(0, 10)), [0, 10 - 3.0e-3]),
            (LineCableModels.resolve(ellipse, Offset(1.0e-3)), [0, -4.0e-3]))
            @test collect(lowest(Contour(shape))) ≈ low atol = 1e-15
        end
        quarter = EllipticArc(Placement(0, 0), 9.0e-3, 3.0e-3, 0, π / 2)
        @test collect(lowest(quarter)) ≈ [9.0e-3, 0] atol = 1e-15
    end

    @testset "the lowest point of an outline agrees with its downward support" begin
        at = Placement(0, 0)
        shapes = (
            ("circle", Circle(1.5e-3)),
            ("annulus", Annulus(20.5e-3, 21.5e-3, at)),
            ("bar", Rectangle(30.0e-3, 12.0e-3)),
            ("bent strand", Sector(8.5e-3, 11.5e-3, 0.3, 0, 0, at)),
            ("rounded wedge", Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, at)),
            ("Milliken segment", Sector(0, 31.15e-3, π / 3, 0, 1.25e-3, at)),
            ("ellipse", Ellipse(9.0e-3, 3.0e-3)),
            ("pentagon",
                Polygon([(0.0, 0.0), (0.004, 0.0), (0.005, 0.003),
                    (0.002, 0.005), (-0.001, 0.003)])),
            ("half ellipse",
                ArcPolygon(Curve{Float64}[
                EllipticArc(at, 3.0e-3, 2.0e-3, 0, π), Segment(
                    (-3.0e-3, 0.0), (3.0e-3, 0.0))])),
            ("perforated bar",
                LineCableModels.PerforatedShape(Rectangle(30.0e-3, 12.0e-3),
                    LineCableModels.AbstractShape{Float64}[Circle(2.0e-3)])),
            ("duct wall",
                LineCableModels.resolve(Ellipse(30.0e-3, 20.0e-3), Offset(3.0e-3))))
        @testset "at ($(place.x), $(place.y)) m" for place in
                                                     (at, Placement(0.3, -1.2, 0.7))
            @testset "$id" for (id, shape) in shapes
                placed = LineCableModels.resolve(place, shape)
                @test lowest(Contour(placed))[2] ≈ -support(placed, 3π / 2) atol = 1e-15
            end
        end
    end

    @testset "the lowest point of a strand of a course in a fill" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        polymer = Material(:insulator, Inf)
        design = build(CableDesign, "mixed_courses",
            Terminal(:core, Region(:core_metal, Circle(0.01), copper)),
            Terminal(:screen,
                Enclosure(Region(:fill, Offset(0.004), polymer),
                    Stack(
                        Repeat(CircularPattern(6; r = Offset(0)),
                            Region(:wires, Circle(0.001), copper)),
                        Repeat(CircularPattern(3; r = Offset(0)),
                            Region(:strands, Rectangle(0.004, 0.002), copper))))))
        system = build(CableSystem, [design], [Placement(0.0, -1.0)];
            connections = [(core = 1, screen = 2)], id = "mixed_courses", line_length = 1.0)
        sector = system.regions[8].shape
        @test collect(lowest(Contour(sector))) ≈
              [0.014cos(2 / 13), -1.0 - 0.014sin(2 / 13)]
    end

    @testset "the unary support is the farthest point" begin
        rectangle = LineCableModels.resolve(Placement(0.01, 0), Rectangle(4.0e-3, 2.0e-3))
        @test support(rectangle) ≈ hypot(0.012, 0.001) rtol = 1e-14
        sector = LineCableModels.resolve(Placement(0.02, 0, π / 2),
            Sector(8.5e-3, 11.5e-3, 0.3, 0, 0, Placement(0, 0)))
        @test support(sector) ≈
              sqrt(0.02^2 + 11.5e-3^2 + 2 * 0.02 * 11.5e-3 * sin(0.15)) rtol = 1e-14
        ellipse = Ellipse(3.0e-3, 2.0e-3, Placement(0, 0.01))
        @test support(ellipse) ≈ 0.012 rtol = 1e-14
        @test support(Contour(ellipse)) ≈ 0.012 rtol = 1e-14
        arc = EllipticArc(Placement(0, 0.01), 3.0e-3, 2.0e-3, π, π / 2)
        @test support(arc) ≈ hypot(3.0e-3, 0.01) rtol = 1e-14
    end

    @testset "a hole that touches its outer leaves no metal" begin
        shape = LineCableModels.PerforatedShape(Rectangle(0.25, 0.125),
            LineCableModels.AbstractShape{Float64}[Circle(0.03125, Placement(0, 0.03125))])
        @test dimension(shape) == 0
    end

    @testset "an elliptic arc's least offset is its ellipse's convex bound" begin
        limit = -dimension(Ellipse(9.0e-3, 3.0e-3)) / 2
        @test EllipticArc(Placement(0, 0), 9.0e-3, 3.0e-3, 0, 2π, nextfloat(limit)).offset ==
              nextfloat(limit)
        @test_throws DomainError EllipticArc(
            Placement(0, 0), 9.0e-3, 3.0e-3, 0, 2π, prevfloat(limit))
        base = EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, π, -π)
        least = -dimension(reverse(base)) / 2
        @test EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, π, -π, nextfloat(least)).offset ==
              nextfloat(least)
        @test_throws DomainError EllipticArc(
            Placement(0, 0), 3.0e-3, 2.0e-3, π, -π, prevfloat(least))
    end

    @testset "a small hole far from the sides keeps the bar's dimension" begin
        shape = LineCableModels.PerforatedShape(Rectangle(30.0e-3, 12.0e-3),
            LineCableModels.AbstractShape{Float64}[Circle(0.5e-3)])
        @test dimension(shape) ≈ 6.0e-3 - 0.5e-3 rtol = 1e-12
    end

    @testset "the area and centroid of a layer infer" begin
        t = 0.5e-3
        inner = Rectangle(4.0e-3, 2.0e-3, Placement(0.01, 0.02))
        layer = LineCableModels.resolve(inner, Offset(t))
        @test @inferred(area(layer)) ≈ 2 * (4.0e-3 + 2.0e-3) * t + π * t^2 rtol = 1e-12
        @test collect(@inferred(centroid(layer))) ≈ [0.01, 0.02] rtol = 1e-12
    end

    @testset "a rectangle placed by its constructor" begin
        at = Placement(0.01, -0.02, π / 7)
        placed = Rectangle(4.0e-3, 2.0e-3, at)
        moved = LineCableModels.resolve(at, Rectangle(4.0e-3, 2.0e-3))
        @test area(placed) == area(moved)
        @test centroid(placed) == centroid(moved)
        @test support(placed) == support(moved)
        @test support(placed, 0.3) == support(moved, 0.3)
    end

    @testset "the common angular range of two arcs" begin
        a = Arc(Placement(0, 0), 2.0e-3, π / 2)
        b = Arc(Placement(0, 0, π / 4), 3.0e-3, π / 2)
        @test extrema(LineCableModels.BaseGeometry.overlap(a, b)) == (π / 4, π / 2)
        c = Arc(Placement(0, 0, π / 2), 3.0e-3, π / 4)
        @test isempty(LineCableModels.BaseGeometry.overlap(
            Arc(Placement(0, 0), 2.0e-3, π /
                                         4), c))
        e = EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2)
        f = EllipticArc(Placement(0, 0), 3.0e-3, 2.0e-3, π / 3, π / 2)
        @test extrema(LineCableModels.BaseGeometry.overlap(e, f)) == (π / 3, π / 2)
    end

    @testset "the apsides of a curve about a point" begin
        @test apsides(Segment((1.0, -1.0), (1.0, 2.0)), (0.0, 0.0)) ==
              (near = 1 / 3, far = 1.0)
        @test apsides(Arc(Placement(3.0, 0.0, π / 2), 1.0, π), (0.0, 0.0)).near == 0.5
        @test @inferred(distance(Arc(Placement(3.0, 0.0, π / 2), 1.0, π))) == 2.0
        @test @inferred(distance(Segment((1.0, -1.0), (1.0, 2.0)))) ≈ 1.0
        @test @inferred(distance(EllipticArc(Placement(0, 10e-3), 3e-3, 2e-3, 0, 2π))) ≈
              8e-3 rtol = 1e-12
    end

    @testset "the distance across shapes with elliptic arcs" begin
        D = Contour(
            Curve{Float64}[EllipticArc(Placement(0, 0), 20e-3, 5e-3, 0, π),
            Segment((-20e-3, 0.0), (-20e-3, -1e-3)),
            Segment((-20e-3, -1e-3), (20e-3, -1e-3)),
            Segment((20e-3, -1e-3), (20e-3, 0.0))],
            Vector{Curve{Float64}}[])
        @test dimension(D) == 0.002
        @test dimension(LineCableModels.PerforatedShape(Rectangle(40e-3, 20e-3),
            [Ellipse(8e-3, 4e-3)])) == 0.006
        @test dimension(LineCableModels.PerforatedShape(Ellipse(30e-3, 20e-3),
            [Circle(10e-3)])) == 0.01
        for shape in (
            LineCableModels.PerforatedShape(Ellipse(30e-3, 20e-3), [Ellipse(10e-3, 15e-3)]),
            LineCableModels.PerforatedShape(Rectangle(80e-3, 20e-3),
                [Ellipse(8e-3, 4e-3, Placement(-20e-3, 0)),
                    Ellipse(8e-3, 4e-3, Placement(20e-3, 0))]))
            @test_throws DomainError dimension(shape)
            @test_throws "needs one ellipse" dimension(shape)
        end
    end

    @testset "the directional support of a shape is its contour's" begin
        for shape in (Circle(3e-3, Placement(1e-3, 2e-3)),
            Annulus(1e-3, 3e-3, Placement(-2e-3, 1e-3)),
            Ellipse(4e-3, 2e-3, Placement(1e-3, -2e-3, 0.3)),
            Polygon([(0.0, 0.0), (4e-3, 0.0), (0.0, 3e-3)], Placement(1e-3, 1e-3, 0.2)))
            for φ in (0, π / 2, π, 3π / 2, 0.7)
                @test support(shape, φ) == support(Contour(shape), φ)
            end
        end
    end

    @testset "the winding about a point and a point in a contour" begin
        @test winding(Segment((1.0, 0.0), (0.0, 1.0)), (0.0, 0.0)) ≈ 0.25
        @test winding(Arc(Placement(0, 0), 1e-3, π), (0, 0.5e-3)) ≈
              (atan(-1, -0.75) + 2π) / 2π
        @test winding(Arc(Placement(0, 0), 1e-3, 2π), (0, 0)) ≈ 1
        @test winding(Arc(Placement(0, 0), 1e-3, 2π), (2e-3, 0)) ≈ 0 atol = 1e-15
        @test (0.0, 0.0) ∈ Contour(Circle(1e-3))
        @test (2e-3, 0.0) ∉ Contour(Circle(1e-3))
        @test (0.0, 0.0) ∉ Contour(Annulus(1e-3, 2e-3))
        @test (1.5e-3, 0.0) ∈ Contour(Annulus(1e-3, 2e-3))
        arc = EllipticArc(Placement(0, 0), 10e-3, 1e-3, π / 2 - 0.1, 0.2, -1.5e-3)
        m = (arc(1 / 2) .+ (arc(0) .+ arc(1)) ./ 2) ./ 2
        u, v = arc(0) .- m, arc(1) .- m
        @test winding(arc, m) ≈
              (atan(u[1] * v[2] - u[2] * v[1], u[1] * v[1] + u[2] * v[2]) + 2π) / 2π
        duct = LineCableModels.resolve(Ellipse(30e-3, 20e-3), Offset(3e-3))
        @test (31.5e-3, 0.0) ∈ Contour(duct)
        @test (0.0, 0.0) ∉ Contour(duct)
        @test (34e-3, 0.0) ∉ Contour(duct)
        @test @inferred((0.0, 0.0) ∈ Contour(Circle(1e-3)))
        @test @inferred((31.5e-3, 0.0) ∈ Contour(duct))
    end

    @testset "a content apart from a hole" begin
        @test isdisjoint(Circle(0.5e-3, Placement(5e-3, 0.0)), Circle(2e-3))
        @test !isdisjoint(Circle(0.5e-3, Placement(2e-3, 0.0)), Circle(2e-3))
        @test isdisjoint(Circle(0.5e-3, Placement(5e-3, 0.0)), Rectangle(2e-3, 2e-3))
        @test !isdisjoint(Circle(0.5e-3, Placement(1e-3, 0.0)), Rectangle(2e-3, 2e-3))
    end

    # A margin k places a content k tolerances clear of the circle, or into it when negative.
    @testset "a shape touches a circle within the tolerance of containment" begin
        r, w = 0.017, 0.5e-3
        tolerance = sqrt(eps(r)) * max(1.0, r)
        A = Circle(r, Placement(0.01, -0.02))
        for (k, touches) in ((0.0, true), (-0.9, true), (0.9, true), (-1.2, false),
            (1.2, false))
            inside = Circle(w, Placement(0.01 + r - w - k * tolerance, -0.02))
            outside = Circle(w, Placement(0.01 + r + w + k * tolerance, -0.02))
            @test istangent(inside, A) == touches
            @test istangent(outside, A) == touches
        end
        sector = Sector(0.01, 0.012, π / 3, 0.0, 0.0, Placement(0.01, -0.02, 0.25))
        @test istangent(sector, Circle(0.01, Placement(0.01, -0.02)))
        @test istangent(sector, Circle(0.012, Placement(0.01, -0.02)))
        @test !istangent(sector, Circle(0.011, Placement(0.01, -0.02)))
        @test @inferred(istangent(Circle(w, Placement(r - w, 0.0)), Circle(r)))
        @test @inferred(istangent(sector, Circle(0.012, Placement(0.01, -0.02))))
        # Uncertainties larger than the tolerance: the nominal values decide.
        @test istangent(Circle(w ± 1e-6, Placement(r - w, 0.0)), Circle(r ± 1e-6))
        @test !istangent(Circle(w ± 1e-6, Placement(r - w + 1.2tolerance, 0.0)),
            Circle(r ± 1e-6))
    end

    # The radii of this annulus lie in two binades. An overlap of 1.5e-9 m exceeds the
    # tolerance of the inner radius alone.
    @testset "an annulus compares a content with each circle within its own tolerance" begin
        A = Annulus(0.0155, 0.017)
        hole, rim = Circle(A.ri), Circle(A.ro)
        inner = sqrt(eps(A.ri)) * max(1.0, A.ri)
        outer = sqrt(eps(A.ro)) * max(1.0, A.ro)
        w = 0.5e-3
        for (k, admitted) in ((-1.2, false), (-0.9, true))
            at_hole = Circle(w, Placement(A.ri + w + k * inner, 0.0))
            at_rim = Circle(w, Placement(A.ro - w - k * outer, 0.0))
            @test (at_hole ⊆ A) == admitted
            @test istangent(at_hole, hole) == admitted
            @test (at_rim ⊆ A) == admitted
            @test istangent(at_rim, rim) == admitted
        end
        into = Circle(w, Placement(A.ri + w - 1.5e-9, 0.0))
        @test !(into ⊆ A)
        @test !istangent(into, hole)
    end

    @testset "the area and moment of a curve about a point and of a path" begin
        segment = Segment((1.0, 0.0), (0.0, 1.0))
        @test @inferred(area(segment, (0.0, 0.0))) == 0.5
        @test collect(@inferred(moment(segment, (0.0, 0.0)))) ≈ [1 / 6, 1 / 6] rtol = 1e-12
        circle = Arc(Placement(0, 0), 1.0, 2π)
        @test @inferred(area(circle, (1.0, 0.0))) ≈ π rtol = 1e-12
        @test collect(@inferred(moment(circle, (1.0, 0.0)))) ≈ [-π, 0.0] rtol = 1e-12
        ellipse = EllipticArc(Placement(0, 0), 3.0, 2.0, 0, 2π, 0.5)
        A = 6π + 0.5 * perimeter(Ellipse(3.0, 2.0)) + 0.25π
        @test @inferred(area(ellipse, (3.5, 0.0))) ≈ A rtol = 1e-12
        @test collect(@inferred(moment(ellipse, (3.5, 0.0)))) ≈ [-3.5A, 0.0] rtol = 1e-12
        path = Contour(Rectangle(4.0e-3, 2.0e-3, Placement(0.01, 0.02))).outer
        @test @inferred(area(path)) ≈ 8.0e-6 rtol = 1e-12
        @test collect(@inferred(moment(path))) ≈ [0.01, 0.02] .* 8.0e-6 rtol = 1e-12
    end

    @testset "a wedge far from the origin keeps its area and centroid" begin
        S0 = Sector(0, 10e-3, 2π / 3, 1e-3 / sin(π / 3), 1.5e-3, Placement(0.0, 0.0))
        S10 = Sector(0, 10e-3, 2π / 3, 1e-3 / sin(π / 3), 1.5e-3, Placement(0.0, 10.0))
        @test area(S10) ≈ area(S0) rtol = 1e-12
        @test all(abs.(centroid(S10) .- (0.0, 10.0) .- centroid(S0)) .<= 1e-14)
        C0 = Contour(LineCableModels.resolve(S0, Offset(0.5e-3)))
        C10 = Contour(LineCableModels.resolve(S10, Offset(0.5e-3)))
        @test area(C10) ≈ area(C0) rtol = 1e-12
        @test all(abs.(centroid(C10) .- (0.0, 10.0) .- centroid(C0)) .<= 1e-14)
    end

    @testset "a curve reversed twice is itself" begin
        for curve in (Segment((0.0, 0.0), (1.0e-3, 2.0e-3)),
            Arc(Placement(1.0e-3, 0, 0.25), 2.0e-3, 0.5),
            EllipticArc(Placement(0, 0, 0.25), 3.0e-3, 2.0e-3, 0.25, 0.5, 0.5e-3))
            @test reverse(reverse(curve)) == curve
        end
    end
end
