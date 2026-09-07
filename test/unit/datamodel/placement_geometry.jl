@testitem "DataModel / explicit placement geometry and area constraints" tags=[:unit] begin
    # A placement is a geometry operation: the caller retains its primitives.
    ellipse = Ellipse(0.2e-3, 0.1e-3)
    ring = Ring(5; r=2e-3, φ0=0.27, span=pi)
    positions = placements(ring, ellipse, nothing)
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

    unmoved = only(placements(nothing, ellipse, nothing))
    @test (unmoved.x, unmoved.y, unmoved.φ) == (0, 0, 0)
    explicit = [Pose2(1e-3, 2e-3, 0.1), Pose2(-2e-3, -1e-3, 0.2)]
    copied = placements(explicit, ellipse, nothing)
    @test copied == explicit
    @test copied !== explicit
    @test_throws ArgumentError placements(Pose2[], ellipse, nothing)

    wire = Disk(0.2e-3)
    for span in (pi, 2pi), gap in (0.0, 0.15)
        available = Ring(capacity(); r=2e-3, span, gap_frac=gap)
        count = capacity(available, wire, nothing)
        packed = placements(Ring(count; r=2e-3, span, gap_frac=gap), wire, nothing)
        @test length(packed) == count
        @test all(hypot(a.x - b.x, a.y - b.y) >= 2wire.r * (1 + gap) - 1e-12
            for (index, a) in enumerate(packed) for b in packed[(index + 1):end])
        @test_throws DomainError placements(
            Ring(count + 1; r=2e-3, span, gap_frac=gap), wire, nothing)
    end

    strip = Rectangle(0.3e-3, 0.1e-3)
    inner = 1e-3
    radius = inner + strip.h / 2
    for fraction in (0.65, 1.0)
        factor = FillFactor(fraction)
        placements_value = placements(Ring(12; r=radius, φ0=0.19), strip, factor)
        outer = sqrt(inner^2 + 12area(strip) / (pi * fraction))
        @test length(placements_value) == 12
        resolved = [value.primitive for value in placements_value]
        @test all(value -> area(value) ≈ area(strip), resolved)
        @test all(value -> r_in(value) ≈ inner, resolved)
        @test all(value -> r_ex(value) ≈ outer, resolved)
        @test sum(area, resolved) / (pi * (outer^2 - inner^2)) ≈ fraction
        @test [value.pose.φ for value in placements_value] ≈
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
end

@testitem "DataModel / contextual capacity respects resolved member extents" tags=[:unit] begin
    sector = Sector(span=2pi / 3, r_base=0.6e-3, r_back=5e-3, fillet=0.2e-3)
    @test capacity(Ring(capacity(); r=0), sector, nothing) == 3
    origins = placements(Ring(3; r=0), sector, nothing)
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
        placed = placements(Ring(count; r=radius, span, gap_frac=gap), item, nothing)
        @test length(placed) == count
        @test all(hypot(a.x - b.x, a.y - b.y) >= width * (1 + gap) - 1e-12
            for (index, a) in enumerate(placed) for b in placed[(index + 1):end])
        @test_throws DomainError placements(
            Ring(count + 1; r=radius, span, gap_frac=gap), item, nothing)
    end
    pattern = Fill(r=3e-3, φ=pi / 6)
    @test placements(pattern, geometry, nothing) ==
        placements(pattern, Disk(support(boundary(geometry))), nothing)
end

@testitem "DataModel / polar placement preserves explicit radial and angular intent" tags=[:unit] begin
    wire = Disk(0.1e-3)
    for origin in (0.0, 0.5e-3), span in (pi, 2pi)
        pattern = Polar(; nr=3, nφ=6, r0=origin, dr=1e-3, φ0=0.27, span)
        poses = placements(pattern, wire, nothing)
        @test length(poses) == (iszero(origin) ? 13 : 18)
        @test poses == placements(pattern, wire, nothing)
        for shell in 0:2
            radius = origin + shell * 1e-3
            row = filter(pose -> isapprox(hypot(pose.x, pose.y), radius; atol=1e-12), poses)
            @test length(row) == (iszero(radius) ? 1 : 6)
            step = span / (span == 2pi ? 6 : 5)
            for (index, pose) in enumerate(row)
                angle = 0.27 + (index - 1) * step
                @test hypot(pose.x - radius * cos(angle), pose.y - radius * sin(angle)) <= 1e-12
                @test pose.φ ≈ angle
            end
        end
        @test all(hypot(a.x - b.x, a.y - b.y) >= 2wire.r - 1e-12
            for (index, a) in enumerate(poses) for b in poses[(index + 1):end])
    end
end
