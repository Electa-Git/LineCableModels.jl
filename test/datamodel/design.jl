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
        @test only(wire.assemblies).regions == [1]
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
                        Enclosure(Region(:core_fill, Circle(10.5e-3), air),
                            Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...))),
                    Region(:core_semicon_inner, Offset(1.0e-3), semicon),
                    Region(:core_insulation, Offset(8.0e-3), polymer),
                    Region(:core_semicon_outer, Offset(1.0e-3), semicon),
                    Terminal(:sheath,
                        Region(:sheath_copper, Annulus(20.5e-3, 21.5e-3), copper)),
                    Region(:sheath_bedding, Offset(3.0e-3), polymer),
                    Terminal(:armour,
                        Enclosure(Region(:armour_fill, Offset(4.0e-3), air), armour)),
                    Region(:armour_jacket, Offset(4.0e-3), polymer))
            end
            @test all(map((a, b) -> a.radii ≈ b.radii, design.lays, contextual.lays))
            @test all(map((a, b) -> support(a.shape) ≈ support(b.shape), design.regions,
                contextual.regions))
            @test length(design.regions) == 85
            @test [row.name for row in design.terminals] == [:core, :sheath, :armour]
            @test [length(row.regions) for row in design.terminals] == [37, 1, 40]
            @test r_ex(design.outer) ≈ 32.5e-3
            @test [length(row.copies) for row in design.repeats] == [6, 12, 18, 40]
            @test [row.repeat for row in design.lays] == [1, 2, 3, 4]
            @test all(radius -> radius ≈ 26.5e-3, design.lays[4].radii)
            @test all(pitch -> pitch ≈ 2 * 10 * 26.5e-3, design.lays[4].pitches)
            @test all(factor -> factor ≈ sqrt(1 + (π / 10)^2), design.lays[4].overlengths)
            @test [design.regions[row.region].source.tag for row in design.fills] ==
                  [:core_fill, :armour_fill]
            @test [length(row.contents) for row in design.fills] == [37, 40]
            core = design.regions[first(design.fills).region].shape
            @test area(core) ≈ π * 10.5e-3^2 - metal
            @test [length(design.regions[row.region].shape.holes) for row in design.fills] ==
                  [37, 40]
            for row in design.fills
                shape = design.regions[row.region].shape
                expected = LineCableModels.DifferenceShape(shape.outer,
                    [design.regions[index].shape for index in row.contents])
                @test area(shape) ≈ area(expected)
                @test centroid(shape) == centroid(expected)
                @test distance(shape) == distance(expected)
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
            Enclosure(Region(:pair_fill, Circle(3.0e-3), air), Terminal(:pair, pair)))
        @test [centroid(placed.shape) for placed in design.regions[1:2]] ==
              [(-2.0e-3, 0.0), (2.0e-3, 0.0)]
    end

    @testset "each copy of a repeat is a scope of its own terminals" begin
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
        separate = build(CableDesign, "three_cores", Enclosure(pipe, cores))
        @test [(row.name, row.scope) for row in separate.terminals] ==
              [(name, [(repeat = 1, copy = k)]) for k in 1:3 for name in (:core, :screen)]
        @test [length(row.regions) for row in separate.terminals] == [7, 1, 7, 1, 7, 1]
        @test [length(row.copies) for row in separate.repeats] == [3, 6, 6, 6]
        @test [row.repeat for row in separate.lays] == [2, 3, 4]
        @test all(row -> hypot(row.axis...) ≈ 3.5e-3, separate.lays)
        @test [length(row.contents) for row in separate.fills] == [7, 7, 7, 30]
        fill = separate.regions[last(separate.fills).region].shape
        @test length(fill.holes) == 3
        @test all(hole -> hole isa Circle, fill.holes)
        @test all(hole -> hole.r ≈ 3.0e-3, fill.holes)
        for (hole, angle) in zip(fill.holes, (0, 2π / 3, 4π / 3))
            @test hole.at.x ≈ 3.5e-3 * cos(angle)
            @test hole.at.y ≈ 3.5e-3 * sin(angle)
        end
        @test area(fill) ≈ π * 6.5e-3^2 - 3π * 3.0e-3^2
        @test [length(separate.regions[row.region].shape.holes)
               for row in separate.fills[1:3]] == [7, 7, 7]
        bonded = build(CableDesign, "bonded", Terminal(:phase, Enclosure(pipe, cores)))
        @test [(row.name, length(row.regions)) for row in bonded.terminals] ==
              [(:phase, 24)]
    end

    @testset "a fill takes the outlines of its outermost contents" begin
        design = build(CableDesign, "layers",
            Enclosure(Region(:core_fill, Circle(5.0e-3), polymer),
                Terminal(:core,
                    Stack(Region(:core_metal, Circle(1.0e-3), copper),
                        Region(:core_inner, Offset(1.0e-3), polymer),
                        Region(:core_outer, Offset(0.5e-3), polymer)))))
        fill = design.regions[only(design.fills).region].shape
        @test length(fill.holes) == 1
        hole = only(fill.holes)
        @test hole isa Circle
        @test hole.r ≈ 2.5e-3
        @test centroid(hole) == (0.0, 0.0)
        @test area(fill) ≈ π * (5.0e-3^2 - 2.5e-3^2)

        hollow = Enclosure(Region(:core_fill, Circle(5.0e-3), polymer),
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
                Enclosure(Region(:core_fill, Circle(4.0e-3), air), Stack(wire, course))))
        box = Region(:box_fill, Rectangle(4.0e-3, 4.0e-3), air)
        @test_throws ArgumentError build(CableDesign, "box",
            Enclosure(box, Terminal(:core, wire)))
    end

    @testset "Float32 inputs" begin
        wire = build(CableDesign, "float32",
            Terminal(:core,
                Region(:core_metal, Circle(1.5f-3),
                    material(kind = :conductor, rho = 1.7241f-8))))
        @test r_ex(wire.outer) ≈ 1.5f-3
        filled = build(CableDesign, "float32_fill",
            Enclosure(Region(:core_fill, Circle(5.0f-3), polymer),
                Terminal(:core,
                    Region(:core_metal, Circle(1.5f-3),
                        material(kind = :conductor, rho = 1.7241f-8)))))
        fill = filled.regions[only(filled.fills).region].shape
        @test length(fill.holes) == 1
        @test eltype(only(fill.holes)) === Float32
    end
end
