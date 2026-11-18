using LineCableModels.DataModel: Ellipse, Polygon, Sector

# A pattern of explicit placements, declared outside the package: a repeat takes any
# pattern.
struct Listed <: AbstractLayoutPattern
    at::Vector{Placement{Float64}}
end
LineCableModels.resolve(_, pattern::Listed, ::Circle) = pattern
function LineCableModels.placements(pattern::Listed, circle::Circle)
    return (placements = pattern.at, shape = circle)
end

@testset "DataModel designs" begin
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.3)
    air = material(kind = :insulator, rho = Inf)

    @testset "a bare wire is one terminal" begin
        wire = build(CableDesign, "bare_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)))
        @test [(row.name, row.scope, row.regions) for row in wire.terminals] ==
              [(:core, [], [1])]
    end

    @testset "an insulation carries no terminal" begin
        wire = build(CableDesign, "insulated_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)),
            Region(:core_insulation, Offset(1.0e-3), polymer))
        @test [(row.name, row.regions) for row in wire.terminals] == [(:core, [1])]
        insulation = wire.regions[2].shape
        @test r_in(insulation) == 0.0425 && r_ex(insulation) ≈ 0.0435
        @test r_ex(wire.outer) ≈ 0.0435
    end

    @testset "the coaxial cables of the references" begin
        semicon = Material(:semicon, 1000.0, 1000.0)
        steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
        for (strand, metal) in ((Circle(1.5e-3), 37π * 1.5e-3^2),
            (Rectangle(3.0e-3, 3.0e-3), π * 1.5e-3^2 + 36 * 9.0e-6))
            # The courses and the armour at their radii, and at no offset from the part inside.
            design,
            contextual = map(((3.0e-3, 6.0e-3, 9.0e-3, 26.5e-3),
                ntuple(_ -> Offset(0), 4))) do radii
                courses = [Lay(Helix(LayRatio(12)),
                               Repeat(CircularPattern(n; r = radii[k]),
                                   Region(Symbol(:core_strands_, k), strand, copper)))
                           for (k, n) in enumerate((6, 12, 18))]
                armour = Lay(Helix(LayRatio(10)),
                    Repeat(CircularPattern(40; r = radii[4]),
                        Region(:armour_wires, Circle(2.0e-3), steel)))
                return build(CableDesign, "coaxial",
                    Terminal(:core,
                        Group(Region(:core_fill, Circle(10.5e-3), air),
                            Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...))),
                    Region(:core_semicon_inner, Offset(1.0e-3), semicon),
                    Region(:core_insulation, Offset(8.0e-3), polymer),
                    Region(:core_semicon_outer, Offset(1.0e-3), semicon),
                    Terminal(:sheath,
                        Region(:sheath_copper, Annulus(20.5e-3, 21.5e-3), copper)),
                    Region(:sheath_bedding, Offset(3.0e-3), polymer),
                    Terminal(:armour,
                        Group(Region(:armour_fill, Offset(4.0e-3), air), armour)),
                    Region(:armour_jacket, Offset(4.0e-3), polymer))
            end
            # The strands of the core and of the armour, each with the overlength and the
            # turns per unit length of its lay.
            core,
            armour = [LineCableModels.resolve(
                          NamedTuple[], Placement(0, 0), 1.0, 0.0, design.tree.items[k])
                      for k in (1, 7)]
            others = [LineCableModels.resolve(
                          NamedTuple[], Placement(0, 0), 1.0, 0.0, contextual.tree.items[k])
                      for k in (1, 7)]
            @test all(hypot(centroid(a.S)...) ≈ hypot(centroid(b.S)...)
            for (strands, declared) in zip((core, armour), others)
            for (a, b) in zip(strands, declared))
            @test all(map((a, b) -> support(a.shape) ≈ support(b.shape), design.regions,
                contextual.regions))
            @test length(design.regions) == 85
            @test [row.name for row in design.terminals] == [:core, :sheath, :armour]
            @test [length(row.regions) for row in design.terminals] == [37, 1, 40]
            @test r_ex(design.outer) ≈ 32.5e-3
            @test [length(row.copies) for row in design.repeats] == [6, 12, 18, 40]
            # The centre wire lies straight. The 36 strands of the courses and the 40 wires
            # of the armour follow helices.
            @test [!iszero(strand.N) for strand in core] == [false; fill(true, 36)]
            @test length(armour) == 40
            @test all(strand -> hypot(centroid(strand.S)...) ≈ 26.5e-3, armour)
            @test all(strand -> inv(strand.N) ≈ 2 * 10 * 26.5e-3, armour)
            @test all(strand -> strand.ℓ ≈ sqrt(1 + (π / 10)^2), armour)
            # The fill of each group is the region after its contents.
            fills = ((38, 1:37), (84, 44:83))
            @test [design.regions[index].source.tag for (index, contents) in fills] ==
                  [:core_fill, :armour_fill]
            @test area(design.regions[38].shape) ≈ π * 10.5e-3^2 - metal
            @test [length(design.regions[index].shape.holes) for (index, contents) in fills] ==
                  [37, 40]
            for (index, contents) in fills
                @test design.regions[index].shape.holes ==
                      [placed.shape for placed in design.regions[contents]]
            end
        end
    end

    @testset "a terminal around terminals, and one name for two parts" begin
        wire = Region(:core_metal, Circle(1.0e-3), copper)
        tube = Region(:core_tube, Offset(1.0e-3), copper)
        nested = build(CableDesign, "nested",
            Terminal(:phase, Stack(Terminal(:inner, wire), Terminal(:outer, tube))))
        @test [(row.name, row.regions) for row in nested.terminals] == [(:phase, [1, 2])]
        shared = build(CableDesign, "shared", Terminal(:core, wire),
            Region(:core_insulation, Offset(1.0e-3), polymer), Terminal(:core, tube))
        @test [(row.name, row.regions) for row in shared.terminals] == [(:core, [1, 3])]
    end

    @testset "a repeat places its copies by any pattern" begin
        pair = Repeat(Listed([Placement(-2.0e-3, 0.0), Placement(2.0e-3, 0.0)]),
            Region(:pair_wires, Circle(1.0e-3), copper))
        design = build(CableDesign, "pair",
            Group(Region(:pair_fill, Circle(3.0e-3), air), Terminal(:pair, pair)))
        @test [centroid(placed.shape) for placed in design.regions[1:2]] ==
              [(-2.0e-3, 0.0), (2.0e-3, 0.0)]
        linear = build(CableDesign, "linear_pair",
            Group(Region(:pair_fill, Circle(3.0e-3), air),
                Terminal(:pair,
                    Repeat(LinearPattern(2, 4.0e-3),
                        Region(:pair_wires, Circle(1.0e-3), copper)))))
        @test [centroid(placed.shape) for placed in linear.regions[1:2]] ==
              [(-2.0e-3, 0.0), (2.0e-3, 0.0)]
    end

    @testset "each copy of a repeat is a scope of its own terminals" begin
        strands = Lay(Helix(LayRatio(12)),
            Repeat(CircularPattern(6; r = 1.0e-3), Region(:core_wires, Circle(0.5e-3), copper)))
        core = Stack(
            Terminal(:core,
                Group(Region(:core_fill, Circle(1.5e-3), air),
                    Stack(Region(:core_centre, Circle(0.5e-3), copper), strands))),
            Region(:core_insulation, Offset(1.0e-3), polymer),
            Terminal(:screen, Region(:core_screen, Offset(0.5e-3), copper)))
        cores = Repeat(CircularPattern(3; r = 3.5e-3), core)
        pipe = Region(:pipe_fill, Circle(6.5e-3), polymer)
        separate = build(CableDesign, "three_cores", Enclosure(pipe, cores))
        @test [(row.name, row.scope) for row in separate.terminals] ==
              [(name, [k + 1]) for k in 1:3 for name in (:core, :screen)]
        @test [length(row.regions) for row in separate.terminals] == [7, 1, 7, 1, 7, 1]
        @test [length(row.copies) for row in separate.repeats] == [3, 6, 6, 6]
        # The resolved repeat keeps one core in its own frame, with the placements of the
        # copies. The wires around the centre wire of the core follow a helix.
        copies = only(separate.tree.items).contents
        @test [!iszero(strand.N)
               for strand in LineCableModels.resolve(
            NamedTuple[], Placement(0, 0), 1.0, 0.0, first(copies.part.items))] ==
              [false; trues(6)]
        @test first(separate.repeats).copies == [1:10, 11:20, 21:30]
        @test all(isapprox(collect((at.x, at.y, at.φ)),
                      [3.5e-3 * cos(φ), 3.5e-3 * sin(φ), φ]; atol = 1.0e-15)
        for (at, φ) in zip(copies.pattern, (k * (2π / 3) for k in 0:2)))
        @test [separate.regions[index].source.tag for index in (8, 18, 28, 31)] ==
              [:core_fill, :core_fill, :core_fill, :pipe_fill]
        fill = separate.regions[31].shape
        @test length(fill.holes) == 3
        @test all(hole -> hole isa Circle, fill.holes)
        @test all(hole -> hole.r ≈ 3.0e-3, fill.holes)
        for (hole, angle) in zip(fill.holes, (0, 2π / 3, 4π / 3))
            @test hole.at.x ≈ 3.5e-3 * cos(angle)
            @test hole.at.y ≈ 3.5e-3 * sin(angle)
        end
        @test area(fill) ≈ π * 6.5e-3^2 - 3π * 3.0e-3^2
        @test [length(separate.regions[index].shape.holes) for index in (8, 18, 28)] ==
              [7, 7, 7]
        bonded = build(CableDesign, "bonded", Terminal(:phase, Enclosure(pipe, cores)))
        @test [(row.name, length(row.regions)) for row in bonded.terminals] ==
              [(:phase, 24)]
    end

    @testset "a laid repeat of a composite part" begin
        insulated = Stack(Terminal(:core, Region(:wire, Circle(1.0e-3), copper)),
            Region(:wire_insulation, Offset(0.5e-3), polymer))
        laid = Lay(Helix(LayRatio(12)), Repeat(CircularPattern(3; r = 6.0e-3), insulated))
        cores = Enclosure(Region(:cores_fill, Circle(7.5e-3), air), laid)
        design = build(CableDesign, "laid_cores", cores)
        @test [placed.source.tag for placed in design.regions] ==
              [repeat([:wire, :wire_insulation], 3); :cores_fill]
        @test [(row.name, row.scope, row.regions) for row in design.terminals] ==
              [(:core, [2], [1]), (:core, [3], [3]), (:core, [4], [5])]
        # The resolved lay keeps its helix over the repeat, with the placements of the copies.
        # Each copy follows the helix at the radius of its placement.
        lay = only(design.tree.items).contents
        @test all(isapprox(collect((at.x, at.y, at.φ)),
                      [6.0e-3 * cos(φ), 6.0e-3 * sin(φ), φ]; atol = 1.0e-15)
        for (at, φ) in zip(lay.part.pattern, (k * (2π / 3) for k in 0:2)))
        @test [row.copies for row in design.repeats] == [[1:2, 3:4, 5:6]]
        radii = [hypot(at.x, at.y) for at in lay.part.pattern]
        @test all(radius -> radius ≈ 6.0e-3, radii)
        @test all(radius -> pitch(lay.helix, radius) ≈ 2 * 12 * 6.0e-3, radii)
        @test all(radius -> overlength(lay.helix, radius) ≈ sqrt(1 + (π / 12)^2), radii)
        @test length(last(design.regions).shape.holes) == 3
    end

    @testset "a fill takes the outlines of its outermost contents" begin
        design = build(CableDesign, "layers",
            Group(Region(:core_fill, Circle(5.0e-3), polymer),
                Terminal(:core,
                    Stack(Region(:core_metal, Circle(1.0e-3), copper),
                        Region(:core_inner, Offset(1.0e-3), polymer),
                        Region(:core_outer, Offset(0.5e-3), polymer)))))
        fill = last(design.regions).shape
        @test length(fill.holes) == 1
        hole = only(fill.holes)
        @test hole isa Circle
        @test hole.r ≈ 2.5e-3
        @test centroid(hole) == (0.0, 0.0)
        @test area(fill) ≈ π * (5.0e-3^2 - 2.5e-3^2)

        hollow = Group(Region(:core_fill, Circle(5.0e-3), polymer),
            Region(:core_tube, Annulus(1.0e-3, 2.0e-3), polymer))
        @test_throws DomainError build(CableDesign, "hollow", hollow)
        @test_throws "a gap or an overlap" build(CableDesign, "hollow", hollow)
    end

    @testset "the checks of a design" begin
        wire = Region(:core_metal, Circle(1.5e-3), copper)
        course = Repeat(CircularPattern(6; r = 3.0e-3), wire)
        @test_throws ArgumentError build(CableDesign, "loose", wire)
        @test_throws DomainError build(CableDesign, "gap", Terminal(:core, wire),
            Region(:core_insulation, Annulus(2.0e-3, 3.0e-3), polymer))
        @test_throws DomainError build(CableDesign, "overlap",
            Terminal(:core, Stack(Region(:core_centre, Circle(2.0e-3), copper), course)))
        @test_throws DomainError build(CableDesign, "protrusion",
            Terminal(:core,
                Group(Region(:core_fill, Circle(4.0e-3), air), Stack(wire, course))))
        tube = Region(:core_tube, Offset(1.0e-3), copper)
        @test_throws "holds the terminals [:a, :b]" build(CableDesign, "two",
            Group(Region(:core_fill, Circle(4.0e-3), air),
                Stack(Terminal(:a, wire), Terminal(:b, tube))))
    end

    @testset "mixed inputs" begin
        design = build(CableDesign, "mixed",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)),
            Region(:core_insulation, Offset(1.0e-3), polymer))
        @test eltype(design) === Float64
        @test all(placed -> eltype(placed.shape) === Float64, design.regions)
        uncertain = build(CableDesign, "mixed_uncertain",
            Terminal(:core,
                Region(:core_metal, Circle(measurement(0.0425, 1.0e-4)), copper)),
            Region(:core_insulation, Offset(1.0e-3), polymer))
        @test eltype(uncertain) === Measurement{Float64}
        filled = build(CableDesign, "filled",
            Group(Region(:core_fill, Circle(5.0e-3), polymer),
                Terminal(:core,
                    Region(:core_metal, Circle(1.5e-3),
                        material(kind = :conductor, rho = 1.7241e-8)))))
        fill = last(filled.regions).shape
        mapped = map(identity, fill)
        @test eltype(only(fill.holes)) === Float64
        @test mapped.outer.r == fill.outer.r
        @test [(hole.r, hole.at.x, hole.at.y, hole.at.φ) for hole in mapped.holes] ==
              [(hole.r, hole.at.x, hole.at.y, hole.at.φ) for hole in fill.holes]
        problem = LineParametersProblem(
            build(CableSystem, [design, design],
                [Placement(0.0, -1.0), Placement(1.0, -1.0)];
                connections = [(core = 1,), (core = 2,)], id = "mixed", line_length = 1.0),
            homogeneous(rho = 100); frequencies = [50.0])
        result = compute(problem, CoaxialFormulation())
        @test result isa LineParameters
    end
    @testset "rectangular polygonal and elliptical containers" begin
        pentagon = Polygon([
            (0.0, 0.0), (0.04, 0.0), (0.05, 0.03), (0.02, 0.05), (-0.01, 0.03)])
        for (container, inside, crossing) in (
            (Rectangle(20e-3, 10e-3), Circle(1.0e-3), Circle(1.0e-3, Placement(0, 4.5e-3))),
            (pentagon, Circle(1.0e-3, Placement(0.02, 0.02)),
                Circle(1.0e-3, Placement(0.02, 0.5e-3))),
            (Ellipse(30e-3, 20e-3), Circle(9.9e-3, Placement(0, 10e-3)),
                Circle(10.1e-3, Placement(0, 10e-3))))
            design = build(CableDesign, "container",
                Enclosure(Region(:container, container, air),
                    Terminal(:core, Region(:core, inside, copper))))
            @test only(design.terminals).regions == [1]
            @test_throws DomainError build(CableDesign, "crossing",
                Enclosure(Region(:container, container, air),
                    Terminal(:core, Region(:core, crossing, copper))))
        end
        sector = Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3),
            1.5e-3, Placement(0, 0))
        for wire in (Circle(1.0e-3, Placement(5.0e-3, 0)),
            Circle(1.0e-3,
                Placement(5.0e-3,
                    (5.0e-3 - 1.0e-3 / sin(π / 3)) * tan(π / 3) -
                    (1.0e-3 + 0.5e-9) / cos(π / 3))))
            design = build(CableDesign, "sector_container",
                Enclosure(Region(:container, sector, air),
                    Terminal(:core, Region(:core, wire, copper))))
            @test only(design.terminals).regions == [1]
        end
        @test_throws DomainError build(CableDesign, "crossing_sector",
            Enclosure(Region(:container, sector, air),
                Terminal(:core,
                    Region(:core, Circle(1.0e-3, Placement(5.0e-3, 5.5e-3)), copper))))
    end

    @testset "contents nest by containment" begin
        sector = Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3),
            1.5e-3, Placement(0, 0))
        design = build(CableDesign, "nested_sector",
            Terminal(:core,
                Group(Region(:outer, Circle(12.0e-3), air),
                    Group(Region(:wedge, sector, air),
                        Region(:wire, Circle(1.0e-3, Placement(5.0e-3, 0)), copper)))))
        fill = last(design.regions).shape
        @test length(fill.holes) == 1
        @test only(fill.holes) == sector
        design = build(CableDesign, "wires_and_tape",
            Region(:inner, Circle(5.0e-3), air),
            Group(Region(:fill, Offset(2.2e-3), air),
                Terminal(:core,
                    Stack(
                        Repeat(CircularPattern(8; r = Offset(0)),
                            Region(:wire, Circle(1.0e-3), copper)),
                        Repeat(CircularPattern(1; r = Offset(0)),
                            Region(:tape, Rectangle(3.0e-3, 0.1e-3), copper))))))
        fill = last(design.regions).shape
        @test length(fill.holes) == 9
    end

    @testset "a one-row duct bank" begin
        concrete = Material(:insulator, Inf, 9)
        design = build(CableDesign, "bank",
            Enclosure(Region(:bank, Rectangle(0.3, 0.2), concrete),
                Repeat(LinearPattern(2, 0.12),
                    Enclosure(Region(:hole, Circle(0.05), air),
                        Terminal(:core, Region(:core, Circle(0.01), copper))))))
        @test [region.source.tag for region in design.regions] ==
              [:core, :hole, :core, :hole, :bank]
        @test [collect(centroid(region.shape)) for region in design.regions] ≈
              [[-0.06, 0.0], [-0.06, 0.0], [0.06, 0.0], [0.06, 0.0], [0.0, 0.0]] rtol = 1e-12
        @test [area(region.shape) for region in design.regions] ≈
              [π * 0.01^2, π * (0.05^2 - 0.01^2), π * 0.01^2,
            π * (0.05^2 - 0.01^2), 0.06 - 2π * 0.05^2] rtol = 1e-12
        @test [(row.copies,) for row in design.repeats] == [([[1, 2], [3, 4]],)]
        @test [length(design.regions[index].shape.holes) for index in (2, 4, 5)] ==
              [1, 1, 2]
        @test [(row.name, row.scope, row.regions) for row in design.terminals] ==
              [(:core, [2], [1]), (:core, [3], [3])]
    end

    @testset "stranded wires inside an ellipse" begin
        design = build(CableDesign, "elliptical_stranded",
            Terminal(:core,
                stranded(copper; center = Circle(1.0e-3), shape = Circle(1.0e-3),
                    n = (6,), lay = LayRatio(12), boundary = Ellipse(4.0e-3, 3.0e-3))))
        @test only(design.terminals).regions == collect(1:7)
        @test_throws DomainError build(CableDesign, "crossing_stranded",
            Terminal(:core,
                stranded(copper; center = Circle(1.0e-3), shape = Circle(1.0e-3),
                    n = (6,), lay = LayRatio(12), boundary = Ellipse(4.0e-3, 2.5e-3))))
    end
    @testset "an elliptical duct with its wall builds" begin
        air = Material(:insulator, Inf)
        pvc = Material(:insulator, Inf, 3)
        design = build(CableDesign, "duct_wall",
            Stack(
                Enclosure(Region(:duct, Ellipse(30e-3, 20e-3), air),
                    Terminal(:core, Region(:core, Circle(0.01), copper))),
                Region(:wall, Offset(3e-3), pvc)))
        wall = last(design.regions).shape
        @test wall isa LineCableModels.PerforatedShape
        @test area(wall) ≈ 0.12 * 1.3221199657742158 * 3e-3 + π * 3e-3^2 rtol = 1e-12
        @test area(design.outer) ≈ area(boundary(wall)) rtol = 1e-12
        @test area(design.outer) ≈
              π * 30e-3 * 20e-3 +
              0.12 * 1.3221199657742158 * 3e-3 + π * 3e-3^2 rtol = 1e-12
    end

    @testset "a duct bank of unequal ducts in two rows" begin
        concrete = Material(:insulator, Inf, 9)
        design = build(CableDesign, "bank",
            Enclosure(Region(:bank, Rectangle(0.3, 0.3), concrete),
                Assembly((
                    Placement(-0.07,
                        -0.07) => Enclosure(Region(:hole, Circle(0.05), air),
                        Terminal(:core, Region(:core, Circle(0.01), copper))),
                    Placement(0.07,
                        0.07) => Enclosure(Region(:hole, Circle(0.04), air),
                        Terminal(:core, Region(:core, Circle(0.01), copper)))))))
        @test [region.source.tag for region in design.regions] ==
              [:core, :hole, :core, :hole, :bank]
        @test [collect(centroid(region.shape)) for region in design.regions[1:4]] ≈
              [[-0.07, -0.07], [-0.07, -0.07], [0.07, 0.07], [0.07, 0.07]] rtol = 1e-12
        @test [area(region.shape) for region in design.regions] ≈
              [π * 0.01^2, π * (0.05^2 - 0.01^2), π * 0.01^2, π * (0.04^2 - 0.01^2),
            0.09 - π * 0.05^2 - π * 0.04^2] rtol = 1e-12
        @test [(at.x, at.y) for (at, part) in only(design.tree.items).contents.members] ==
              [(-0.07, -0.07), (0.07, 0.07)]
        @test isempty(design.repeats)
        @test [length(design.regions[index].shape.holes) for index in (2, 4, 5)] ==
              [1, 1, 2]
        @test [(row.name, row.scope, row.regions) for row in design.terminals] ==
              [(:core, [2], [1]), (:core, [3], [3])]
    end

    # The sector cores touch along their sides, each in its frame at its centroid. A core
    # rotated by 0.02 rad about the center enters its neighbour.
    @testset "the members of an assembly lie apart" begin
        sector = Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(0, 0))
        c = first(centroid(sector))
        core = Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(-c, 0))
        for (angles, apart) in
            (((0, 2π / 3, 4π / 3), true), ((0, 2π / 3 - 0.02, 4π / 3), false))
            part = Stack(
                Enclosure(Region(:fill, Circle(11.0e-3), air),
                    Assembly(Tuple(Placement(c * cos(φ),
                                       c * sin(φ),
                                       φ) => Stack(
                                       Terminal(name, Region(:core, core, copper)),
                                       Region(:insulation, Offset(1.0e-3), polymer))
                    for (φ, name) in zip(angles, (:phase_a, :phase_b, :phase_c))))),
                Region(:jacket, Offset(3.0e-3), polymer))
            if apart
                @test length(build(CableDesign, "three_sector_cores", part).terminals) == 3
            else
                @test_throws DomainError build(CableDesign, "three_sector_cores", part)
            end
        end
    end

    @testset "equal parts built apart are equal and have one hash" begin
        strands() = Lay(Helix(LayRatio(12)),
            Repeat(CircularPattern(6; r = 3.0e-3), Region(:strand, Circle(1.0e-3), copper)))
        cable() = Stack(
            Terminal(:core,
                Group(Region(:core_fill, Circle(4.0e-3), polymer),
                    Stack(Region(:center, Circle(1.0e-3), copper), strands()))),
            Region(:insulation, Offset(2.0e-3), polymer))
        pipe() = Enclosure(Region(:filler, Circle(20.0e-3), polymer),
            Assembly((
                Placement(8.0e-3, 0.0) => cable(), Placement(-8.0e-3, 0.0) => cable())))
        # A lay, a group, a terminal, a stack, an assembly, an enclosure and a resolved repeat.
        nodes() = (strands(), first(cable().items).part, first(cable().items), cable(),
            pipe().contents, pipe(), Repeat([Placement(8.0e-3, 0.0)], cable()))
        for (a, b) in zip(nodes(), nodes())
            @test a == b
            @test hash(a) == hash(b)
        end
        @test Group(Region(:fill, Circle(4.0e-3), polymer), Region(:a, Circle(1.0e-3), copper)) !=
              Group(Region(:fill, Circle(4.0e-3), polymer), Region(:b, Circle(1.0e-3), copper))
        # Resolved parts: a fill has the outlines of its contents as holes, and the
        # insulation of a bar has the arcs of its outline.
        trees() = (
            build(CableDesign, "pair", pipe(), Region(:jacket, Offset(1.0e-3), polymer)).tree,
            build(CableDesign, "bar",
                Terminal(:core, Region(:bar, Rectangle(10.0e-3, 5.0e-3), copper)),
                Region(:insulation, Offset(1.0e-3), polymer)).tree,
            Region(:wedge, Polygon([(0.0, 0.0), (1.0e-3, 0.0), (0.0, 1.0e-3)]), copper))
        for (a, b) in zip(trees(), trees())
            @test a == b
            @test hash(a) == hash(b)
        end
    end
end
