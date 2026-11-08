using LineCableModels.BaseGeometry: Segment, Arc, EllipticArc, Curve, Contour, draw,
                                    clearance, apsides, lowest

struct RecordingTarget
    curves::Vector{Union{Segment, Arc}}
end

function LineCableModels.BaseGeometry.draw(target::RecordingTarget, segment::Segment)
    push!(target.curves, segment)
    return segment
end

function LineCableModels.BaseGeometry.draw(target::RecordingTarget, arc::Arc)
    push!(target.curves, arc)
    return arc
end

@testset "BaseGeometry" begin
    @testset "circle quadrants ignore the placement rotation" begin
        for φ in (0.0, 0.7, -π / 3)
            circle = Circle(2.0e-3, Placement(1.0e-3, -3.0e-3, φ))
            contour = Contour(circle)
            @test length(contour.outer) == 4
            @test isempty(contour.holes)
            for (arc, angle) in zip(contour.outer, (0.0, π / 2, Float64(π), 3π / 2))
                @test arc isa Arc
                @test arc.r == circle.r
                @test arc.span == π / 2
                @test (arc.at.x, arc.at.y, arc.at.φ) ==
                      (circle.at.x, circle.at.y, angle)
                @test all(isapprox.(arc(0),
                    (circle.at.x + circle.r * cos(angle),
                        circle.at.y + circle.r * sin(angle)); atol = 1e-15))
            end
        end
    end

    @testset "annulus quadrants" begin
        annulus = Annulus(1.0e-3, 2.0e-3, Placement(3.0e-3, -1.0e-3, 0.7))
        contour = Contour(annulus)
        @test length(contour.outer) == 4
        @test length(contour.holes) == 1
        @test length(only(contour.holes)) == 4
        for (path, radius) in ((contour.outer, annulus.ro),
            (only(contour.holes), annulus.ri))
            for (arc, angle) in zip(path, (0.0, π / 2, Float64(π), 3π / 2))
                @test arc isa Arc
                @test arc.r == radius
                @test arc.span == π / 2
                @test (arc.at.x, arc.at.y, arc.at.φ) ==
                      (annulus.at.x, annulus.at.y, angle)
            end
        end
    end

    @testset "sector paths and rotated corners" begin
        for φ in (0.0, 0.6)
            sector = LineCableModels.Sector(
                1.0e-3, 3.0e-3, 1.2, 0, 0, Placement(2.0e-3, -4.0e-3, φ))
            contour = Contour(sector)
            @test length(contour.outer) == 4
            @test isempty(contour.holes)
            @test contour.outer[1] isa Arc
            @test contour.outer[2] isa Segment
            @test contour.outer[3] isa Arc
            @test contour.outer[4] isa Segment
            @test contour.outer[1].span == sector.span
            @test contour.outer[3].span == -sector.span
            corners = ((sector.ro, φ - sector.span / 2),
                (sector.ro, φ + sector.span / 2),
                (sector.ri, φ + sector.span / 2),
                (sector.ri, φ - sector.span / 2))
            for (k, (radius, angle)) in enumerate(corners)
                @test all(isapprox.(contour.outer[k](0),
                    (sector.at.x + radius * cos(angle),
                        sector.at.y + radius * sin(angle)); atol = 1e-15))
                @test all(isapprox.(contour.outer[k](1),
                    contour.outer[mod1(k + 1, 4)](0); atol = 1e-15))
            end
        end
        sector = LineCableModels.Sector(0.0, 3.0e-3, 1.2, 0, 0, Placement(0.0, 0.0))
        contour = Contour(sector)
        @test length(contour.outer) == 3
        @test isempty(contour.holes)
        @test contour.outer[1] isa Arc
        @test contour.outer[2] isa Segment
        @test contour.outer[3] isa Segment
        @test contour.outer[2](1) == (0.0, 0.0)
        for k in eachindex(contour.outer)
            @test all(isapprox.(contour.outer[k](1),
                contour.outer[mod1(k + 1, 3)](0); atol = 1e-15))
        end
        for ri in (0.0, 1.0e-3)
            at = Placement(2.0e-3, -4.0e-3, 0.6)
            full = Contour(LineCableModels.Sector(ri, 3.0e-3, 2π, 0, 0, at))
            expected = iszero(ri) ? Contour(Circle(3.0e-3, at)) :
                       Contour(Annulus(ri, 3.0e-3, at))
            @test full.outer == expected.outer
            @test full.holes == expected.holes
        end
    end

    @testset "perforated contours and the three-core fill" begin
        at = Placement(0.0, 0.0)
        shape = LineCableModels.PerforatedShape(Annulus(1.0e-3, 6.0e-3, at),
            [Circle(0.5e-3, Placement(3.0e-3, 0.0))])
        contour = Contour(shape)
        @test contour.outer == Contour(shape.outer).outer
        @test length(contour.holes) == 2
        @test contour.holes[1] == only(Contour(shape.outer).holes)
        @test contour.holes[2] == Contour(only(shape.holes)).outer

        copper = Material(:conductor, 1.7241e-8)
        polymer = Material(:insulator, 1.0e14, 2.3)
        air = material(kind = :insulator, rho = Inf)
        strands = Lay(Helix(LayRatio(12)),
            Repeat(CircularPattern(6; r = 1.0e-3), Region(:core_wires, Circle(0.5e-3), copper)))
        core = Stack(
            Terminal(:core,
                Enclosure(Region(:core_fill, Circle(1.5e-3), air),
                    Stack(Region(:core_centre, Circle(0.5e-3), copper), strands))),
            Region(:core_insulation, Offset(1.0e-3), polymer),
            Terminal(:screen, Region(:core_screen, Offset(0.5e-3), copper)))
        cores = Repeat(CircularPattern(3; r = 3.5e-3), core)
        pipe = Region(:pipe_fill, Circle(6.5e-3), polymer)
        design = build(CableDesign, "three_cores", Enclosure(pipe, cores))
        fill = design.regions[last(design.fills).region].shape
        contour = Contour(fill)
        @test contour.outer == Contour(fill.outer).outer
        @test length(contour.holes) == 3
        for (hole, angle) in zip(contour.holes, (0, 2π / 3, 4π / 3))
            @test length(hole) == 4
            for (arc, start) in zip(hole, (0.0, π / 2, Float64(π), 3π / 2))
                @test arc isa Arc
                @test arc.at.x ≈ 3.5e-3 * cos(angle)
                @test arc.at.y ≈ 3.5e-3 * sin(angle)
                @test arc.r ≈ 3.0e-3
                @test arc.at.φ == start
                @test arc.span == π / 2
            end
        end
    end

    @testset "curve evaluation and floating types" begin
        at = Placement(2.0e-3, -1.0e-3, 0.7)
        for span in (1.2, -1.2)
            arc = Arc(at, 3.0e-3, span)
            for t in (0, 0.5, 1)
                angle = at.φ + t * span
                @test all(isapprox.(arc(t),
                    (at.x + arc.r * cos(angle), at.y + arc.r * sin(angle)); atol = 1e-15))
            end
        end
        segment = Segment((1.0e-3, 2.0e-3), (5.0e-3, -2.0e-3))
        for (t, point) in ((0, (1.0e-3, 2.0e-3)),
            (0.25, (2.0e-3, 1.0e-3)), (1, (5.0e-3, -2.0e-3)))
            @test all(isapprox.(segment(t), point; atol = 1e-15))
        end
        @test Curve{Float64} === Union{Segment{Float64}, Arc{Float64}, EllipticArc{Float64}}
        sector = LineCableModels.Sector(1.0e-3, 3.0e-3, 1.2, 0, 0, at)
        @test eltype(Contour(sector).outer) === Curve{Float64}
        @test Segment((0, 1.0), (2.0, 3)).q == (2.0, 3.0)
        @test Arc(Placement(0.0, 0.0), 1, 0.5).span == 0.5
    end

    @testset "signed contour areas" begin
        for shape in (Circle(3.0e-3, Placement(1.0e-3, -2.0e-3, 0.4)),
            Annulus(1.0e-3, 3.0e-3, Placement(1.0e-3, -2.0e-3)),
            LineCableModels.Sector(1.0e-3, 3.0e-3, 1.2, 0, 0, Placement(0.0, 0.0)),
            LineCableModels.Sector(0.0, 3.0e-3, 1.2, 0, 0, Placement(0.0, 0.0)),
            LineCableModels.Sector(1.0e-3, 3.0e-3, 1.2, 0, 0,
                Placement(2.0e-3, -4.0e-3, 0.6)))
            contour = Contour(shape)
            signed = Float64[]
            for path in (contour.outer, contour.holes...)
                points = [curve(t) for curve in path
                          for t in range(0, 1; length = 10_000)]
                push!(signed,
                    sum(points[k][1] * points[mod1(k + 1, length(points))][2] -
                        points[mod1(k + 1, length(points))][1] * points[k][2]
                    for k in eachindex(points)) / 2)
            end
            @test all(>(0), signed)
            @test isapprox(first(signed) - sum(signed[2:end]), area(shape); rtol = 1e-6)
        end
    end

    @testset "invalid curves and paths" begin
        open = Curve{Float64}[
        Segment(
            (0.0, 0.0), (1.0, 0.0))]
        empty = Curve{Float64}[]
        holes = Vector{Curve{Float64}}[]
        @test_throws DomainError Contour(open, holes)
        @test_throws DomainError Contour(empty, holes)
        @test_throws DomainError Segment((0.0, 0.0), (0.0, 0.0))
        @test_throws DomainError Arc(Placement(0.0, 0.0), 0.0, 1.0)
        @test_throws DomainError Arc(Placement(0.0, 0.0), 1.0, 0.0)
        @test_throws DomainError Arc(Placement(0.0, 0.0), 1.0, 2π + 0.1)
    end

    @testset "draw a design and a system" begin
        problem = only(problem
        for (name, problem) in REFERENCES
        if name == "coaxial_circular_strands")
        for subject in (only(problem.system.designs), problem.system)
            target = RecordingTarget(Union{Segment, Arc}[])
            result = draw(target, subject)
            contours = [Contour(placed.shape) for placed in subject.regions]
            expected = [curve for contour in contours
                        for path in (contour.outer, contour.holes...) for curve in path]
            @test length(result) == length(subject.regions)
            @test target.curves == expected
            @test length(target.curves) == sum(length(contour.outer) +
                      sum(length, contour.holes; init = 0) for contour in contours)
            for (entry, contour) in zip(result, contours)
                @test entry isa NamedTuple{(:outer, :holes)}
                @test entry.outer == contour.outer
                @test entry.holes == contour.holes
            end
        end
    end

    @testset "the apsides, support and lowest point of a curve along a direction" begin
        ellipse = EllipticArc(Placement(0, 0), 2.0e-3, 1.0e-3, 0, 2π)
        @test @inferred(apsides(ellipse, 0)) == (near = 0.0, far = 0.5)
        @test @inferred(apsides(ellipse, π / 2)) == (near = 0.25, far = 0.75)
        @test @inferred(apsides(Arc(Placement(0, 0), 1.0e-3, π / 2), π / 4)).near ≈ 0.5
        @test @inferred(apsides(Segment((0.0, 0.0), (1.0e-3, 1.0e-3)), π)) ==
              (near = 0.0, far = 1.0)
        segment = Segment((1.0e-3, -2.0e-3), (-3.0e-3, 4.0e-3))
        arc = Arc(Placement(1.0e-3, 2.0e-3, 0.0), 3.0e-3, π)
        elliptic = EllipticArc(
            Placement(1.0e-3, -2.0e-3, 0.3), 4.0e-3, 2.0e-3, 0, 2π, 0.5e-3)
        for φ in (0.4, 1.0, 2.5)
            c, s = cos(φ), sin(φ)
            h = hypot(4.0e-3 * cos(φ - 0.3), 2.0e-3 * sin(φ - 0.3))
            @test @inferred(support(segment, φ)) ≈
                  max(1.0e-3 * c - 2.0e-3 * s, -3.0e-3 * c + 4.0e-3 * s) rtol = 1e-12
            @test @inferred(support(arc, φ)) ≈ 1.0e-3 * c + 2.0e-3 * s + 3.0e-3 rtol = 1e-12
            @test @inferred(support(elliptic, φ)) ≈
                  1.0e-3 * c - 2.0e-3 * s + h + 0.5e-3 rtol = 1e-12
        end
        @test @inferred(lowest(segment)) == (1.0e-3, -2.0e-3)
        @test collect(@inferred(lowest(Arc(Placement(1.0e-3, 2.0e-3, π), 3.0e-3, π)))) ≈
              [1.0e-3, -1.0e-3] atol = 1e-15
        h = hypot(4.0e-3 * cos(3π / 2 - 0.3), 2.0e-3 * sin(3π / 2 - 0.3))
        @test @inferred(lowest(elliptic))[2] ≈ -2.0e-3 - h - 0.5e-3 atol = 1e-15
    end

    @testset "the clearance of a segment or an arc and an elliptic arc" begin
        upper = EllipticArc(Placement(0, 0), 20e-3, 5e-3, 0, π)
        p, q, d, s, t = @inferred clearance(Segment((-30e-3, 0.0), (30e-3, 0.0)), upper)
        @test d == 0
        @test t ∈ (0.0, 1.0)
        p, q, d = clearance(Segment((0.0, 0.0), (0.0, 10e-3)), upper)
        @test d == 0
        @test all(abs.(p .- (0.0, 5e-3)) .<= 1e-15)
        @test all(abs.(q .- (0.0, 5e-3)) .<= 1e-15)
        above = Segment((-10e-3, 6e-3), (10e-3, 6e-3))
        p, q, d, s, t = clearance(above, upper)
        @test d ≈ 1e-3 rtol = 1e-12
        @test t ≈ 0.5
        @test @inferred(clearance(upper, above)) == (q, p, d, t, s)
        half = Arc(Placement(0, 12e-3, -π), 4e-3, π)
        p, q, d, s, t = @inferred clearance(half, upper)
        @test d ≈ 3e-3 rtol = 1e-12
        @test all(abs.(q .- (0.0, 5e-3)) .<= 1e-15)
        @test @inferred(clearance(upper, half)) == (q, p, d, t, s)
        @test clearance(Arc(Placement(0, 12e-3, -π), 8e-3, π), upper)[3] == 0
    end

    @testset "the clearance of two elliptic arcs of one ellipse" begin
        at = Placement(-0.00059794295776892294, -0.00069360917295612293, -1.8835078656211701)
        a, b = 0.025493804607145745, 0.0034161657108289691
        e = EllipticArc(at, a, b, -1.3060388229981847, 0.29299243356178317,
            0.00035255212983780062)
        f = EllipticArc(at, a, b, -3.1025237020383027, -2.933366926794454,
            0.0004067968507438287)
        p, q, d, s, t = @inferred clearance(e, f)
        @test d ≈ 0.006529847843022441 rtol = 1e-12
        @test (p, q) == (e(s), f(t))
        @test hypot((p .- q)...) ≈ d rtol = 1e-12
        deep = EllipticArc(Placement(0, 0), 20e-3, 5e-3, 0.2, π - 0.4, -2e-3)
        other = EllipticArc(Placement(0, 0), 20e-3, 5e-3, π + 0.2, π - 0.4, 0.0)
        @test_throws DomainError clearance(deep, other)
        @test_throws "least radius of curvature" clearance(deep, other)
        @test_throws "needs one ellipse" clearance(other,
            EllipticArc(Placement(1e-3, 0), 20e-3, 5e-3, 0.2, π - 0.4))
    end

    @testset "the clearance of every pair of curves" begin
        curves = (
            Segment((-30e-3, 9e-3), (30e-3, 9e-3)), Arc(Placement(0, 12e-3, -π), 4e-3, π),
            EllipticArc(Placement(0, 0), 20e-3, 5e-3, 0, π))
        widths = [0 0 4e-3; 0 0 3e-3; 4e-3 3e-3 0]
        for (i, a) in enumerate(curves), (j, b) in enumerate(curves)

            @test clearance(a, b)[3] ≈ widths[i, j] atol = 1e-15
        end
    end

    @testset "a bracket halved to adjacent parameters" begin
        lo, hi = @inferred LineCableModels.BaseGeometry.bisect(x -> x < 0.3, 0.0, 1.0)
        @test lo < 0.3 <= hi
        @test nextfloat(lo) == hi
    end

    # A fresh process reads the inference before any other call has compiled the curves.
    @testset "the directional support of a contour infers in a fresh process" begin
        script = "using LineCableModels; " *
                 "using LineCableModels.BaseGeometry: Contour, support; " *
                 "using LineCableModels.DataModel: Sector; " *
                 "println(only(Base.return_types(support, (Contour{Float64}, Float64)))); " *
                 "println(only(Base.return_types(support, (Sector{Float64}, Float64))))"
        command = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script`
        @test readlines(command) == ["Float64", "Float64"]
    end
end
