@testitem "DataModel / stranding / complete-course capacity at roundoff boundaries" tags=[:unit, :parametric] begin
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

    # The 320/380 kV cable examples have exactly enough copper for 127 wires.
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
    wires = filter(placed -> placed.terminal === :core, design.geometry.regions)
    @test length(wires) == 127
    @test all(placed -> isapprox(area(placed), area(wire); rtol=5e-12), wires)
    @test sum(area, wires) ≈ area(boundary) rtol=5e-12
    @test all(placed -> DM.support(placed.primitive) <= boundary.r * (1 + 5e-6), wires)

    # Natural wires fill the first course by fit and add six per course. No two wires overlap.
    strand = Disk(0.5e-3)
    for (radius, counts) in ((0.5e-3, [6, 12, 18]), (0.25e-3, [4, 10, 16]),
        (1.0e-3, [9, 15, 21]))
        sites = DM.circular_courses(Disk(4.0e-3), Disk(radius), strand, false)
        @test [count(member -> member.course == k, sites) for k in 1:3] == counts
        disks = [((0.0, 0.0), radius); [(member.site, strand.r) for member in sites]]
        @test all(hypot((left[1] .- right[1])...) >= (left[2] + right[2]) * (1 - 1e-12)
                  for (index, left) in enumerate(disks) for right in disks[(index + 1):end])
    end
end

@testitem "DataModel / stranding / declared course counts" tags=[:unit, :parametric] begin
    import LineCableModels.DataModel as DM
    copper = Material(kind=:conductor, rho=1.7241e-8)
    oil = Material(kind=:insulator, rho=1e13, eps_r=2.2)
    wire, boundary = Disk(0.5e-3), Disk(3.5e-3)

    # A declared schedule lies on the natural radii, and the fill takes the room left over.
    e17 = build(CableDesign, "declared-schedule",
        terminal(:core, stranded(copper; shape=wire, n=(5, 11, 16), boundary)))
    wires = filter(placed -> placed.terminal === :core, e17.geometry.regions)
    e17_courses = [only(placed.placement.patterns).pattern.course for placed in wires]
    @test [count(==(course), e17_courses) for course in 0:3] == [1, 5, 11, 16]
    @test all(hypot(centroid(placed.primitive)...) ≈ (2course) * wire.r
              for (placed, course) in zip(wires, e17_courses) if course > 0)
    centers = [centroid(placed.primitive) for placed in wires]
    @test all(hypot((left .- right)...) >= 2wire.r * (1 - 1e-12)
              for (index, left) in enumerate(centers) for right in centers[(index + 1):end])
    @test only(e17.tree.member.items).contents.member.items[end].pattern == [5, 11, 16]

    # With equal center and strand sizes, a base count of six is the inferred inventory.
    for (compact, disk) in ((false, boundary), (true, Disk(sqrt(37) * wire.r)))
        @test build(CableDesign, "base-count",
            terminal(:core, stranded(copper; shape=wire, n=6, compact, boundary=disk))
        ).geometry == build(CableDesign, "base-count",
            terminal(:core, stranded(copper; shape=wire, compact, boundary=disk))).geometry
    end
    compacted = build(CableDesign, "declared-compacted", terminal(:core,
        stranded(copper; shape=wire, n=(5, 10), compact=true, boundary=Disk(sqrt(18) * wire.r))))
    @test [count(placed -> only(placed.placement.patterns).pattern.course == course,
               compacted.geometry.regions[1:16]) for course in 0:2] == [1, 5, 10]
    @test all(placed -> area(placed.primitive) ≈ area(wire), compacted.geometry.regions[1:16])

    # After a part, the declared courses surround its outline.
    around = build(CableDesign, "declared-around", filler(oil, Disk(6e-3)),
        terminal(:core, stranded(copper; shape=Disk(1.5e-3), n=(15, 20),
            boundary=Disk(12e-3))))
    around_wires = filter(placed -> placed.terminal === :core, around.geometry.regions)
    @test [count(placed -> only(placed.placement.patterns).pattern.course == course,
               around_wires) for course in 1:2] == [15, 20]
    @test all(placed -> hypot(centroid(placed.primitive)...) >= 7.5e-3 * (1 - 1e-12),
        around_wires)

    # Declared rectangular courses that fill the disk boundary give the inferred core.
    rod, rod_center = Rectangle(0.75e-3, 0.35e-3), Disk(0.55e-3)
    @test build(CableDesign, "declared-rods", terminal(:core, stranded(copper;
        center=rod_center, shape=rod, n=[6, 8, 11, 14, 17, 20, 23, 25, 18], boundary))
    ).geometry == build(CableDesign, "declared-rods", terminal(:core, stranded(copper;
        center=rod_center, shape=rod, boundary))).geometry

    @test_throws "ring members overlap" build(CableDesign, "over-full",
        terminal(:core, stranded(copper; shape=wire, n=(7, 12), boundary)))
    @test_throws "reach beyond the disk boundary" build(CableDesign, "beyond",
        terminal(:core, stranded(copper; shape=wire, n=(6, 12, 18, 24), boundary)))
    @test_throws "exceeds its authoritative boundary area" build(CableDesign, "overflow",
        terminal(:core, stranded(copper; shape=wire, n=(6, 12, 18), compact=true,
            boundary=Disk(sqrt(30) * wire.r))))
    @test_throws "more than one strand of the disk boundary unfilled" build(CableDesign,
        "few-rods", terminal(:core, stranded(copper; center=rod_center, shape=rod,
            n=[6, 8, 11], boundary)))
    @test_throws "exceed the disk boundary" build(CableDesign, "many-rods",
        terminal(:core, stranded(copper; center=rod_center, shape=rod,
            n=[6, 8, 11, 14, 17, 20, 23, 25, 19], boundary)))
    @test_throws "n must be a positive base count" stranded(copper; shape=wire, n=0, boundary)

    # In a non-circular boundary, declared courses follow the center wire, and each wire keeps
    # its area. The sector holds 55.5 wire areas, so a base count of 5 gives four courses
    # (1 + 50 wires) and not five (1 + 75).
    for (core_boundary, strand, n, counts) in (
        (Rectangle(10e-3, 6e-3), wire, [5, 11, 17], [1, 5, 11, 17]),
        (Sector(2pi / 3, 0.6e-3, 5e-3, 0.2e-3), Disk(0.35e-3), 5, [1, 5, 10, 15, 20]))
        design = build(CableDesign, "declared-non-circular",
            terminal(:core, stranded(copper; shape=strand, n, boundary=core_boundary)))
        strands = filter(placed -> placed.terminal === :core, design.geometry.regions)
        courses = [only(placed.placement.patterns).pattern.course for placed in strands]
        @test [count(==(course), courses) for course in 0:(length(counts) - 1)] == counts
        @test all(placed -> area(placed.primitive) ≈ area(strand), strands)
    end
    @test_throws "exceed the area of the boundary, which holds 76" build(CableDesign,
        "over-full-rectangle", terminal(:core, stranded(copper; shape=wire,
            n=[6, 12, 18, 24, 30], boundary=Rectangle(10e-3, 6e-3))))
end

@testitem "DataModel / stranding / non-circular boundaries" tags=[:unit, :parametric] begin
    import LineCableModels.DataModel as DM
    copper = Material(kind=:conductor, rho=1.7241e-8)
    oil = Material(kind=:insulator, rho=1e13, eps_r=2.2)
    wire = Disk(0.5e-3)
    rectangle, ellipse = Rectangle(10e-3, 6e-3), Ellipse(4e-3, 2.5e-3)
    # A five-pointed star about the origin, star-shaped about its centroid and not convex.
    star = Polygon([(r * cos(k * pi / 5), r * sin(k * pi / 5))
                    for (k, r) in zip(0:9, Iterators.cycle((5e-3, 2.5e-3)))])

    # Each geometric boundary contains a center wire and the complete 6k courses that its area
    # admits. An area of 76.4 wires gives four courses (61 wires), and 40 or 46.8 give three
    # (37 wires).
    designs = [build(CableDesign, "non-circular",
                   terminal(:core, stranded(copper; shape=wire, boundary=core_boundary)))
               for core_boundary in (rectangle, ellipse, star)]
    wires = [filter(placed -> placed.terminal === :core, design.geometry.regions)
             for design in designs]
    for (core_boundary, design, strands, counts) in zip((rectangle, ellipse, star), designs,
        wires, ([1, 6, 12, 18, 24], [1, 6, 12, 18], [1, 6, 12, 18]))
        courses = [only(placed.placement.patterns).pattern.course for placed in strands]
        @test [count(==(course), courses) for course in 0:(length(counts) - 1)] == counts
        @test all(placed -> isapprox(area(placed.primitive), area(wire); rtol=1e-12), strands)
        # The fill takes the area of the geometric boundary that the declared metal leaves.
        fills = filter(placed -> placed.terminal === nothing, design.geometry.regions)
        @test sum(area, fills) ≈ area(core_boundary) - sum(counts) * area(wire) rtol=1e-12
    end

    # Every strand vertex lies inside its geometric boundary.
    @test all(wires[1]) do placed
        all(placed.primitive.points) do point
            x, y = placed.primitive.at(point)
            abs(x) <= 5e-3 * (1 + 1e-12) && abs(y) <= 3e-3 * (1 + 1e-12)
        end
    end
    @test all(wires[2]) do placed
        all(placed.primitive.points) do point
            x, y = placed.primitive.at(point)
            (x / 4e-3)^2 + (y / 2.5e-3)^2 <= 1 + 1e-12
        end
    end
    # A star vertex at angle θ lies on the inner side of the star edge that spans θ.
    @test all(wires[3]) do placed
        all(placed.primitive.points) do point
            x, y = placed.primitive.at(point)
            k = floor(Int, mod(atan(y, x), 2pi) / (pi / 5))
            first_point = star.points[mod1(k + 1, 10)]
            last_point = star.points[mod1(k + 2, 10)]
            (last_point[1] - first_point[1]) * (y - first_point[2]) -
            (last_point[2] - first_point[2]) * (x - first_point[1]) >= -1e-15
        end
    end

    # As a sector core, a rectangular core flattens to the area of its core boundary, with the
    # resistance of its 61 declared wires.
    row = only(DM.flatten(designs[1], Float64).conductors)
    @test row.r_ex ≈ sqrt(area(rectangle) / pi) rtol=1e-12
    @test DM.tubular_resistance(row.r_in, row.r_ex, row.material.rho) ≈
          copper.rho / (61area(wire)) rtol=1e-12

    # A shell around the rectangular core is the rectangle 2t wider and 2t higher, with
    # sharp corners. Each layer adds 2t(w + h) + 4t² to the rectangle inside it. The jacket
    # wraps the outer rectangle of the insulation.
    pe = Material(kind=:insulator, rho=1e14, eps_r=2.3)
    insulated = build(CableDesign, "insulated-rectangle",
        terminal(:core, stranded(copper; shape=wire, boundary=rectangle)),
        insulation(pe; t=1e-3), jacket(pe; t=0.5e-3))
    layers = [only(filter(placed -> placed.source.tag === tag, insulated.geometry.regions))
              for tag in (:insulation, :jacket)]
    @test layers[1].primitive.inner == rectangle
    @test layers[2].primitive.inner == layers[1].primitive.outer
    @test all(zip(layers, (1e-3, 0.5e-3))) do (layer, t)
        inner = layer.primitive.inner
        area(layer.primitive) ≈ 2t * (inner.w + inner.h) + 4t^2
    end
    @test all(isfinite, compute(CableConstantsProblem(insulated)).R)

    # Around a polygon, every edge moves outward by t, and adjacent moved edges meet at a
    # sharp corner. A layer is the shoelace area of its outer polygon less that of its inner
    # one. A regular hexagon of apothem a grows to apothem a + t. Its layer is
    # 2√3((a + t)² - a²).
    hexagon = Polygon([(4e-3 * cos(k * pi / 3), 4e-3 * sin(k * pi / 3)) for k in 0:5])
    designs = map((hexagon, star, ellipse)) do core_boundary
        build(CableDesign, "insulated-core",
            terminal(:core, stranded(copper; shape=wire, boundary=core_boundary)),
            insulation(pe; t=1e-3))
    end
    layers = [only(filter(placed -> placed.source.tag === :insulation,
                  design.geometry.regions)).primitive for design in designs[1:2]]
    @test all(layers) do layer
        shoelace = [abs(sum(p[1] * q[2] - q[1] * p[2]
                            for (p, q) in zip(points, circshift(points, -1)))) / 2
                    for points in (layer.outer.points, layer.inner.points)]
        area(layer) ≈ shoelace[1] - shoelace[2]
    end
    apothem = 4e-3 * cos(pi / 6)
    @test area(first(layers)) ≈ 2sqrt(3) * ((apothem + 1e-3)^2 - apothem^2)
    # Every shell flattens to the annulus between the radii of equal area, so the insulated
    # hexagon and ellipse cores compute.
    @test all(design -> all(isfinite, compute(CableConstantsProblem(design)).R),
        designs[[1, 3]])
    # The 2 mm wide notch of this outline closes under a 1.5 mm shell.
    notched = Polygon([(0.0, 0.0), (10e-3, 0.0), (10e-3, 10e-3), (6e-3, 10e-3),
        (6e-3, 4e-3), (4e-3, 4e-3), (4e-3, 10e-3), (0.0, 10e-3)])
    @test_throws "too large for this Polygon outline" DM.resolve(notched, DM.Shell(1.5e-3))

    # This polygon is an 8 mm square with a slot cut in from one side. Its centroid lies in
    # the slot, so it does not see every edge.
    slotted = Polygon([(0.0, 0.0), (8e-3, 0.0), (8e-3, 1e-3), (1e-3, 1e-3), (1e-3, 7e-3),
        (8e-3, 7e-3), (8e-3, 8e-3), (0.0, 8e-3)])
    @test_throws "must be star-shaped about its centroid" build(CableDesign, "slotted",
        terminal(:core, stranded(copper; shape=wire, boundary=slotted)))
    for core_boundary in (rectangle, ellipse, star)
        @test_throws "deforms its wires intrinsically; omit compact" stranded(copper;
            shape=wire, compact=true, boundary=core_boundary)
        @test_throws "rectangular strands take a Disk boundary" stranded(copper;
            center=Disk(0.5e-3), shape=Rectangle(0.75e-3, 0.35e-3), boundary=core_boundary)
        @test_throws "Rectangle, Ellipse or Polygon boundary starts at the cable center" build(
            CableDesign, "non-circular-after-part", filler(oil, Disk(1e-3)),
            terminal(:core, stranded(copper; shape=wire, boundary=core_boundary)))
    end
end

@testitem "DataModel / stranding / rectangular strands fill their boundary" tags=[:unit, :parametric] begin
    import LineCableModels.DataModel as DM
    copper = Material(kind=:conductor, rho=1.7241e-8)
    center, strand, boundary = Disk(0.55e-3), Rectangle(0.75e-3, 0.35e-3), Disk(3.5e-3)
    design = build(CableDesign, "rectangular-courses", terminal(:core,
        stranded(copper; center, shape=strand, boundary)))
    strips = filter(placed -> placed.source.primitive == strand, design.geometry.regions)
    courses = [only(placed.placement.patterns).pattern.course for placed in strips]
    @test [count(==(course), courses) for course in 1:9] == [6, 8, 11, 14, 17, 20, 23, 25, 18]
    @test outer_radius(design) == maximum(placed -> placed.primitive.ro, strips) == boundary.r
    # The course heights scale by area from the packed radius to the disk boundary.
    packed = sqrt(center.r^2 + 142area(strand) / pi)
    scale = (boundary.r^2 - center.r^2) / (packed^2 - center.r^2)
    @test all(placed -> area(placed.primitive) ≈ scale * area(strand), strips)
    # The strips of a course close it, the first course starts on the center wire, and each
    # next course starts where the previous one ends.
    radii = [only(unique((placed.primitive.ri, placed.primitive.ro)
                         for placed in strips[courses .== course])) for course in 1:9]
    @test first(first(radii)) == center.r
    @test all(first(radii[course]) == last(radii[course - 1]) for course in 2:9)
    @test all(sum(placed -> placed.primitive.span, strips[courses .== course]) ≈ 2pi
              for course in 1:9)
    # The resistance is that of the declared rectangular strands.
    row = only(DM.flatten(design, Float64).conductors)
    @test DM.tubular_resistance(row.r_in, row.r_ex, row.material.rho) ≈
          copper.rho / (area(center) + 142area(strand)) rtol=1e-12
end

@testitem "DataModel / stranding / a formation after a part fills the ring" tags=[:unit, :parametric] begin
    import LineCableModels.DataModel as DM
    aluminum = Material(kind=:conductor, rho=2.8264e-8)
    copper = Material(kind=:conductor, rho=1.7241e-8)
    gel = Material(kind=:insulator, rho=1e12, eps_r=2.2)
    other_gel = Material(kind=:insulator, rho=1e9, eps_r=4.0)
    oil = Material(kind=:insulator, rho=1e13, eps_r=2.2)
    pe = Material(kind=:insulator, rho=Inf, eps_r=2.3)

    # An OPGW has a gel-filled aluminium tube, then rectangular strands bent around it. The
    # courses follow the rule at the cable center, with the tube in place of the center wire.
    strand = Rectangle(3e-3, 2e-3)
    opgw_designs = map((gel, other_gel)) do filling
        build(CableDesign, "opgw", filler(filling, Disk(2.5e-3)),
            terminal(:core, sheath(aluminum; t=0.5e-3),
                stranded(aluminum; shape=strand, boundary=Disk(7e-3), lay=LayRatio(12))),
            insulation(pe; t=1e-3))
    end
    opgw = first(opgw_designs)
    tube = only(filter(placed -> placed.source.tag === :sheath, opgw.geometry.regions))
    strips = filter(placed -> placed.source.primitive == strand, opgw.geometry.regions)
    courses = [only(entry.pattern.course for entry in placed.placement.patterns
                    if entry.pattern isa DM.BoundedPlacement) for placed in strips]
    @test [count(==(course), courses) for course in 1:2] == [8, 12]
    @test all(placed -> placed.primitive isa DM.BentStrip, strips)
    @test all(placed -> placed.primitive.ri == tube.primitive.ro,
        strips[courses .== 1])
    # The course heights scale so that the outer course ends on the disk boundary.
    packed = sqrt(tube.primitive.ro^2 + 20area(strand) / pi)
    scale = (7e-3^2 - tube.primitive.ro^2) / (packed^2 - tube.primitive.ro^2)
    @test all(placed -> area(placed.primitive) ≈ scale * area(strand), strips)
    @test maximum(placed -> placed.primitive.ro, strips) == 7e-3
    @test all(placed -> placed.primitive.ri == maximum(r -> r.primitive.ro, strips[courses .== 1]),
        strips[courses .== 2])
    @test last(opgw.geometry.regions).primitive.ri == 7e-3

    # The tube and the rectangular strands form one terminal, and the design computes. The
    # gel only sets the inner radius of that conductor. A different gel gives the same Z and Y.
    @test all(isfinite, compute(CableConstantsProblem(opgw)).R)
    opgw_parameters = map(opgw_designs) do design
        system = build(LineCableSystem, design, Pose2(0.0, 10.0);
            connections=Dict(:core => 1))
        compute(LineParametersProblem(system; earth_props=homogeneous(rho=100.0),
            frequencies=[50.0]))
    end
    @test all(isfinite, first(opgw_parameters).Z) && all(isfinite, first(opgw_parameters).Y)
    @test first(opgw_parameters).Z == last(opgw_parameters).Z
    @test first(opgw_parameters).Y == last(opgw_parameters).Y
    # A declared annulus is the same tube.
    declared_tube = build(CableDesign, "opgw", filler(gel, Disk(2.5e-3)),
        terminal(:core, core(aluminum, Annulus(2.5e-3, 3e-3)),
            stranded(aluminum; shape=strand, boundary=Disk(7e-3), lay=LayRatio(12))),
        insulation(pe; t=1e-3))
    for table in (:conductors, :dielectrics)
        @test getfield(DM.flatten(declared_tube, Float64), table) ==
              getfield(DM.flatten(opgw, Float64), table)
    end

    # In a hollow conductor, round wires surround an oil duct. The first course holds the 15
    # wires that fit around the duct, and each next course six more.
    wire = Disk(1.5e-3)
    oil_duct = filler(oil, Disk(6e-3))
    hollow = build(CableDesign, "hollow-conductor", oil_duct,
        terminal(:core, stranded(copper; shape=wire, boundary=Disk(15e-3), lay=LayRatio(14))),
        insulation(pe; t=2e-3))
    wires = filter(placed -> placed.terminal === :core, hollow.geometry.regions)
    courses = [only(entry.pattern.course for entry in placed.placement.patterns
                    if entry.pattern isa DM.BoundedPlacement) for placed in wires]
    @test [count(==(course), courses) for course in 1:3] == [15, 21, 27]
    centers = [centroid(placed.primitive) for placed in wires]
    for (course, radius) in zip(1:3, (7.5e-3, 10.5e-3, 13.5e-3))
        @test all(center -> hypot(center...) ≈ radius, centers[courses .== course])
    end
    @test all(hypot((left .- right)...) >= 2wire.r * (1 - 1e-12)
              for (index, left) in enumerate(centers) for right in centers[(index + 1):end])
    filled = only(filter(placed -> placed.source.tag === :stranded_fill,
        hollow.geometry.regions))
    @test (filled.primitive.outer.ri, filled.primitive.outer.ro) == (6e-3, 15e-3)
    @test area(filled.primitive) ≈ pi * (15e-3^2 - 6e-3^2) - 63area(wire)
    @test all(isfinite, compute(CableConstantsProblem(hollow)).R)

    # A sheath and rectangular armour form one terminal.
    armoured = build(CableDesign, "sheath-armour", terminal(:core, core(copper; r=3e-3)),
        insulation(pe; t=2e-3),
        terminal(:screen, sheath(aluminum; t=0.5e-3),
            armor(aluminum; shape=Rectangle(1e-3, 0.5e-3), n=10)),
        insulation(pe; t=1e-3))
    @test all(isfinite, compute(CableConstantsProblem(armoured)).R)

    # A formation after an insulated core starts at the insulation and computes.
    insulated = build(CableDesign, "stranded-neutral", terminal(:core, core(copper; r=3e-3)),
        insulation(pe; t=2e-3),
        terminal(:neutral, stranded(copper; shape=Disk(0.5e-3), boundary=Disk(7e-3))),
        insulation(pe; t=1e-3))
    neutral = last(DM.flatten(insulated, Float64).conductors)
    @test neutral.r_in == 5e-3
    @test neutral.r_ex ≈ 7e-3
    @test all(isfinite, compute(CableConstantsProblem(insulated)).R)

    @test_throws "has no center of its own" build(CableDesign, "duct-center", oil_duct,
        terminal(:core, stranded(copper; center=Disk(1e-3), shape=wire, boundary=Disk(15e-3))))
    @test_throws "takes natural wires" build(CableDesign, "duct-compact", oil_duct,
        terminal(:core, stranded(copper; shape=wire, compact=true, boundary=Disk(15e-3))))
    @test_throws "takes a Disk boundary around a circular outline" build(CableDesign,
        "duct-sector", oil_duct, terminal(:core, stranded(copper; shape=Disk(0.5e-3),
            boundary=Sector(2pi / 3, 1.5e-3, 8e-3, 0.6e-3))))
    @test_throws "concentric with and outside the preceding part" build(CableDesign,
        "duct-inside", oil_duct, terminal(:core, stranded(copper; shape=wire, boundary=Disk(5e-3))))
    @test_throws "less than one wire diameter around the center wire or the preceding part" build(
        CableDesign, "duct-thin", oil_duct,
        terminal(:core, stranded(copper; shape=wire, boundary=Disk(8e-3))))
    for shape in (wire, strand)
        @test_throws "takes a Disk boundary around a circular outline" build(CableDesign,
            "after-rectangle", filler(oil, Rectangle(4e-3, 4e-3)),
            terminal(:core, stranded(aluminum; shape, boundary=Disk(15e-3))))
    end
    @test_throws "requires center=Disk(...)" build(CableDesign, "rods-at-center",
        terminal(:core, stranded(aluminum; shape=strand, boundary=Disk(7e-3))))
end

@testitem "DataModel / stranding / clipped disks conserve area inside cells" tags=[:unit, :measurements] begin
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

@testitem "DataModel / stranding / sector bundle identity, containment and symmetry" tags=[:unit, :datamodel] begin
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
            members, shapes, courses = DM.mapped_courses(boundary, Disk(a))
            @test courses == 3
            @test length(shapes) == 37
            @test [count(m -> m.course == k, members) for k in 0:3] == [1, 6, 12, 18]
            @test all(isapprox.(first(members).site,
                DM.polygon_centroid(DM.fan_polygon(boundary))))
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
                    center = DM.centroid(shapes[i])
                    any(indices) do j
                        mirror = DM.centroid(shapes[j])
                        hypot(center[1] - mirror[1], center[2] + mirror[2]) < tolerance
                    end
                end
            end
        end

        pose = Pose2(0.012, -0.008, 0.43)
        original, shapes, _ = DM.mapped_courses(boundary, Disk(radius))
        moved, rotated, _ = DM.mapped_courses(DM.resolve(pose, boundary), Disk(radius))
        @test getproperty.(moved, :course) == getproperty.(original, :course)
        @test all(eachindex(shapes)) do i
            expected = DM.centroid(DM.resolve(pose, shapes[i]))
            actual = DM.centroid(rotated[i])
            hypot(expected[1] - actual[1], expected[2] - actual[2]) < 1e-9
        end
    end
end
