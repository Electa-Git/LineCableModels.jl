@testitem "DataModel / explicit placement geometry and area constraints" tags=[:unit, :datamodel] begin
    # A placement is a geometry operation: the caller retains its primitives.
    ellipse = Ellipse(0.2e-3, 0.1e-3)
    ring = Ring(5; r=2e-3, φ0=0.27, span=pi)
    positions = placements(ring, ellipse, nothing).poses
    @test length(positions) == 5
    for (index, position) in enumerate(positions)
        angle = 0.27 + (index - 1) * pi / 4
        @test position.x ≈ 2e-3 * cos(angle)
        @test position.y ≈ 2e-3 * sin(angle)
        @test position.φ ≈ angle
    end
    @test_throws DomainError placements(Ring(20; r=0.2e-3), ellipse, nothing)
    @test_throws ArgumentError placements(Ring(5), ellipse, nothing)
    @test_throws ArgumentError placements(Ring(capacity(); r=2e-3), ellipse, nothing)

    unmoved = only(placements(nothing, ellipse, nothing).poses)
    @test (unmoved.x, unmoved.y, unmoved.φ) == (0, 0, 0)
    explicit = [Pose2(1e-3, 2e-3, 0.1), Pose2(-2e-3, -1e-3, 0.2)]
    copied = placements(explicit, ellipse, nothing)
    @test copied.poses == explicit
    @test copied.shape === ellipse
    @test_throws ArgumentError placements(Pose2[], ellipse, nothing)

    wire = Disk(0.2e-3)
    for span in (pi, 2pi), gap in (0.0, 0.15)
        available = Ring(capacity(); r=2e-3, span, gap_frac=gap)
        count = capacity(available, wire, nothing)
        packed = placements(Ring(count; r=2e-3, span, gap_frac=gap), wire, nothing).poses
        @test length(packed) == count
        @test all(hypot(a.x - b.x, a.y - b.y) >= 2wire.r * (1 + gap) - 1e-12
            for (index, a) in enumerate(packed) for b in packed[(index + 1):end])
        @test_throws DomainError placements(
            Ring(count + 1; r=2e-3, span, gap_frac=gap), wire, nothing)
    end
    @test_throws ArgumentError capacity(Ring(capacity()), wire, nothing)

    strip = Rectangle(0.3e-3, 0.1e-3)
    inner = 1e-3
    radius = inner + strip.h / 2
    for fraction in (0.65, 1.0)
        factor = FillFactor(fraction)
        course = placements(Ring(12; r=radius, φ0=0.19), strip, factor)
        outer = sqrt(inner^2 + 12area(strip) / (pi * fraction))
        @test length(course.poses) == 12
        @test course.shape isa LineCableModels.DataModel.BentStrip
        @test area(course.shape) ≈ area(strip)
        @test r_in(course.shape) ≈ inner
        @test r_ex(course.shape) ≈ outer
        @test 12area(course.shape) / (pi * (outer^2 - inner^2)) ≈ fraction
        @test [pose.φ for pose in course.poses] ≈
            [0.19 + (index - 1) * 2pi / 12 for index in 1:12]
        count = capacity(Ring(capacity(); r=radius), strip, factor)
        available_area = pi * ((radius + strip.h / 2)^2 - inner^2)
        @test count * area(strip) <= fraction * available_area
        @test (count + 1) * area(strip) > fraction * available_area
    end
    @test_throws ArgumentError placements(Ring(12), strip, FillFactor(1))
    @test_throws ArgumentError placements(
        Ring(capacity(); r=radius), strip, FillFactor(1))
    @test_throws DomainError capacity(
        Ring(capacity(); r=radius, gap_frac=0.1), strip, FillFactor(1))
    @test capacity(Ring(capacity(); r=0), strip, FillFactor(0.8)) == 0
    @test_throws ArgumentError capacity(Ring(capacity()), strip, FillFactor(0.8))
end

@testitem "DataModel / contextual capacity respects resolved member extents" tags=[:unit, :parametric] begin
    sector = Sector(span=2pi / 3, r_base=0.6e-3, r_back=5e-3, fillet=0.2e-3)
    @test capacity(Ring(capacity(); r=0), sector, nothing) == 3
    origins = placements(Ring(3; r=0), sector, nothing).poses
    @test length(origins) == 3
    @test all(pose -> iszero(pose.x) && iszero(pose.y), origins)
    @test getproperty.(origins, :φ) ≈ [0.0, 2pi / 3, 4pi / 3]
    copper = Material(kind=:conductor, rho=1.72e-8)
    dielectric = Material(kind=:insulator, rho=Inf, eps_r=2.3)
    geometry = resolve(EmptyBoundary(),
        terminal(:core, core(copper; r=0.5e-3), insulation(dielectric; t=0.2e-3)))
    for item in (sector, geometry), span in (pi, 2pi), gap in (0.0, 0.1)
        shape = item isa Sector ? resolve(EmptyBoundary(), item) : boundary(item)
        width = item isa Sector ? support(shape, pi / 2) + support(shape, -pi / 2) :
            2support(shape)
        radius = 20e-3
        available = Ring(capacity(); r=radius, span, gap_frac=gap)
        count = capacity(available, item, nothing)
        minimum_angle = 2asin(width * (1 + gap) / (2radius))
        expected = floor(Int, span / minimum_angle) + (span == 2pi ? 0 : 1)
        @test count == expected
        placed = placements(Ring(count; r=radius, span, gap_frac=gap), item, nothing).poses
        @test length(placed) == count
        @test all(hypot(a.x - b.x, a.y - b.y) >= width * (1 + gap) - 1e-12
            for (index, a) in enumerate(placed) for b in placed[(index + 1):end])
        @test_throws DomainError placements(
            Ring(count + 1; r=radius, span, gap_frac=gap), item, nothing)
    end
end

@testitem "DataModel / contextual ring dispatch preserves inferred geometry" tags=[:unit, :datamodel] begin
    const DM = LineCableModels.DataModel
    copper = Material(kind=:conductor, rho=1.72e-8)
    for T in (Float32, Float64, BigFloat)
        part = DM.Region(:wire, DM.Disk(T(0.001)), copper)
        child = DM.resolve(DM.EmptyBoundary(), part)
        for context in (DM.EmptyBoundary(), DM.Disk(T(0.01)))
            fixed = Ring(6; r=T(0.02))
            resolved = @inferred DM._contextual_pattern(fixed, part, child, nothing, context)
            @test resolved.n == 6
            @test resolved.r == fixed.r
            available = Ring(capacity(); r=T(0.02))
            filled = @inferred DM._contextual_pattern(available, part, child, nothing, context)
            @test filled.n == capacity(available, part.primitive, nothing)
            derived = @inferred DM._contextual_pattern(Ring(6), part, child, nothing, context)
            @test derived.r == (context isa DM.EmptyBoundary ? zero(T) : support(context)) + part.primitive.r
        end
        @test (@inferred DM._contextual_pattern(nothing, part, child, nothing, DM.EmptyBoundary())) === nothing
    end
end

@testitem "DataModel / a ring course bends rectangles and keeps every other member" tags=[:unit, :datamodel] begin
    const DM = LineCableModels.DataModel
    copper = Material(kind=:conductor, rho=1.72e-8)
    wire = Disk(0.2e-3)
    @test placements(Ring(6; r=2e-3), wire, nothing).shape === wire
    @test placements(Ring(6; r=2e-3), wire, nothing).poses isa Vector{Pose2{Float64}}

    # Every rectangle on a ring becomes the same strip, bent about the ring center.
    rectangle = Rectangle(0.4e-3, 0.2e-3)
    geometry = DM.resolve(DM.EmptyBoundary(),
        Group(:g, Region(:w, rectangle, copper); pattern=Ring(6; r=2e-3)))
    shapes = [placed.primitive for placed in geometry.regions]
    @test length(shapes) == 6
    @test all(shape -> shape isa DM.BentStrip, shapes)
    @test all(shape -> r_in(shape) ≈ 1.9e-3 && r_ex(shape) ≈ 2.1e-3, shapes)
    @test all(shape -> area(shape) ≈ area(rectangle), shapes)
end

@testitem "DataModel / minimum radius distinguishes centered holes from displaced boundaries" tags=[:unit, :datamodel] begin
    const DM = LineCableModels.DataModel
    for shape in (Annulus(1.0, 2.0), DM.BentStrip(1.0, 2.0, 0.5))
        @test (@inferred DM._minimum_radius(shape)) == 1.0
    end
    for moved in (Disk(1.0, Pose2(3.0, 0.0)), Ellipse(1.0, 0.5, Pose2(3.0, 0.0)),
            Annulus(0.5, 1.0, Pose2(3.0, 0.0)))
        @test (@inferred DM._minimum_radius(moved)) ≈ 2.0
    end
    @test (@inferred DM._minimum_radius(Disk(1.0))) == 0.0
    @test r_in(Disk(1.0, Pose2(3.0, 0.0))) == 0.0
end

@testitem "DataModel / lifting / an invalid scalar placement input throws MethodError" tags=[:unit, :datamodel] begin
    # Each lift accepts the scalar type of its constructor or a finite source, so an input
    # of another type reaches no method.
    @test_throws MethodError Pose2("a", 0.0, 0.0)
    @test_throws MethodError Pose2(0.0, :b)
    @test_throws MethodError LayRatio("a")
    @test_throws MethodError Pitch(:p)
    @test_throws MethodError LayAngle(1 + 2im)
    @test_throws MethodError Helix(10.0)
    @test_throws MethodError Helix(LayRatio(10.0); dir = true)
    @test_throws MethodError FillFactor("a")
    @test_throws MethodError Ring(2.5)
    @test Pose2(Grid((0.0, 1.0)), 2.0) isa Gridspace{Pose2}
    @test LayRatio(Grid((10.0, 12.0))) isa Gridspace{LayRatio}
    @test Helix(LayRatio(10.0); dir = Grid((1, -1))) isa Gridspace{Helix}
    @test Ring(Grid((6, 12)); r = 1e-3) isa Gridspace{Ring}
end
