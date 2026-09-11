@testitem "DataModel / stranding / complete-course capacity at roundoff boundaries" tags=[:unit] begin
    import LineCableModels.DataModel as DM

    for F in (Float32, Float64), courses in 1:12, strand_area in F.((1e-8, 1.0, 1e8))
        inventory = F(3courses * (courses + 1))
        for ratio in (inventory, inventory - 4eps(inventory), inventory + 4eps(inventory))
            @test DM.course_count(ratio * strand_area, strand_area) == courses
        end
        # A real area deficit must still exclude the incomplete outer course.
        insufficient = (inventory - 256eps(inventory)) * strand_area
        if courses == 1
            @test_throws DomainError DM.course_count(insufficient, strand_area)
        else
            @test DM.course_count(insufficient, strand_area) == courses - 1
        end
    end

    # The 320/380 kV gauntlet cases have exactly enough copper for 127 wires.
    # Before the fix, the inverse returned 5.999999999999999 and lost 36 wires.
    radius = 3.66e-3 / 2
    wire = Disk(radius)
    boundary = Disk(sqrt(127) * radius)
    sites = DM.circular_courses(boundary, wire, wire, true)
    @test length(sites) == 126
    @test [count(site -> site.course == k, sites) for k in 1:6] == 6 .* (1:6)

    copper = Material(kind=:conductor, rho=1.7241e-8)
    design = @cable "full-127-wire-core" begin
        @terminal :core begin
            stranded(copper; shape=wire, boundary=boundary, compact=true)
        end
    end
    wires = filter(region -> region.terminal === :core, design.geometry.regions)
    @test length(wires) == 127
    @test all(region -> isapprox(area(region), area(wire); rtol=5e-12), wires)
    @test sum(area, wires) ≈ area(boundary) rtol=5e-12
    @test all(region -> DM.support(region.primitive) <= boundary.r * (1 + 5e-6), wires)
end

@testitem "DataModel / stranding / clipped disks conserve area inside cells" tags=[:unit] begin
    using Measurements: measurement, value, uncertainty
    import LineCableModels.DataModel as DM

    cell = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
    for target in (pi / 4, 3.0, 3.9, 4.0)
        strand = DM.area_preserving_strand(cell, target)
        @test DM.signed_polygon_area(strand) ≈ target rtol=5e-13
        @test all(point -> maximum(abs, point) <= 1 + 1e-12, strand)
    end
    circular = DM.area_preserving_strand(cell, pi / 4)
    radii = [hypot(point...) for point in circular]
    @test maximum(radii) - minimum(radii) < 1e-12
    @test maximum(radii) ≈ 0.5 rtol=3e-4
    @test DM.area_preserving_strand(cell, 4.0) == cell
    @test_throws DomainError DM.area_preserving_strand(cell, 4.1)

    # A tolerance crossing must never extrapolate a polygon edge.
    near_edge = [(-1.0, 1e-15), (1.0, 1e-12), (0.0, -1.0)]
    clipped = DM.clip_halfplane(near_edge, (0.0, 1.0), 0.0)
    @test all(point -> -1 <= point[1] <= 1, clipped)

    target = measurement(3.0, 0.03)
    propagated = DM.signed_polygon_area(DM.area_preserving_strand(cell, target))
    @test value(propagated) ≈ value(target) rtol=5e-13
    @test uncertainty(propagated) ≈ uncertainty(target) rtol=1e-6
end

@testitem "DataModel / stranding / sector bundle identity, containment and symmetry" tags=[:unit] begin
    import LineCableModels.DataModel as DM

    function vertices(shape)
        c, s = cos(shape.at.φ), sin(shape.at.φ)
        return [(shape.at.x + c * p[1] - s * p[2],
                 shape.at.y + s * p[1] + c * p[2]) for p in shape.points]
    end

    # Separating-axis check uses the actual conductor polygons, independently
    # of the power diagram which allocated their space.
    function separated(left, right, tolerance)
        for polygon in (left, right), index in eachindex(polygon)
            p, q = polygon[index], polygon[mod1(index + 1, length(polygon))]
            normal = (q[2] - p[2], p[1] - q[1])
            norm = hypot(normal...)
            norm <= tolerance && continue
            a = [point[1] * normal[1] + point[2] * normal[2] for point in left]
            b = [point[1] * normal[1] + point[2] * normal[2] for point in right]
            min(maximum(a) - minimum(b), maximum(b) - minimum(a)) <=
                tolerance * norm && return true
        end
        return false
    end

    boundaries = (
        (Sector(deg2rad(118), 1.5e-3, 8e-3, 0.6e-3), 0.55e-3, 0.7423123065),
        (Sector(2pi / 3, 0.6e-3, 5e-3, 0.2e-3), 0.35e-3, 0.6661730464),
        (Sector(pi / 2, 0.6e-3, 5e-3, 0.2e-3), 0.35e-3, 0.8852644112)
    )
    for (definition, radius, expected_fill) in boundaries
        boundary = DM.resolve(DM.EmptyBoundary(), definition)
        for a in (radius, sqrt(0.93 * area(boundary) / (37pi)))
            members, shapes, courses = DM.sector_courses(boundary, Disk(a))
            @test courses == 3
            @test length(shapes) == 37
            @test [count(m -> m.course == k, members) for k in 0:3] == [1, 6, 12, 18]
            @test all(isapprox.(first(members).site,
                DM.polygon_centroid(DM.sector_polygon(boundary))))
            @test all(k -> [m.member for m in members if m.course == k] == 1:(6k), 1:3)
            @test maximum(abs.(area.(shapes) ./ (pi * a^2) .- 1)) < 5e-12
            @test 37pi * a^2 / area(boundary) ≈
                  (a == radius ? expected_fill : 0.93) rtol=1e-9
            polygons = vertices.(shapes)
            tolerance = 1e-9 * definition.r_back
            @test all(polygon -> all(p -> DM.clearance(boundary, p) >= -tolerance,
                                    polygon), polygons)
            @test all(separated(polygons[i], polygons[j], tolerance)
                      for i in eachindex(polygons) for j in (i + 1):length(polygons))

            for k in 0:3
                indices = findall(m -> m.course == k, members)
                @test all(indices) do i
                    centre = DM.centroid(shapes[i])
                    any(indices) do j
                        mirror = DM.centroid(shapes[j])
                        hypot(centre[1] - mirror[1], centre[2] + mirror[2]) < tolerance
                    end
                end
            end
        end

        pose = Pose2(0.012, -0.008, 0.43)
        original, shapes, _ = DM.sector_courses(boundary, Disk(radius))
        moved, rotated, _ = DM.sector_courses(DM.resolve(pose, boundary), Disk(radius))
        @test getproperty.(moved, :course) == getproperty.(original, :course)
        @test all(eachindex(shapes)) do i
            expected = DM.centroid(DM.resolve(pose, shapes[i]))
            actual = DM.centroid(rotated[i])
            hypot(expected[1] - actual[1], expected[2] - actual[2]) < 1e-9
        end
    end
end
