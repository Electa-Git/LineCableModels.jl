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
            Region(:core_insulation, Shell(1.0e-3), polymer))
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
            courses = [Lay(Helix(LayRatio(12)),
                           Repeat(CircularPattern(n; r = 3.0e-3k),
                               Region(Symbol(:core_strands_, k), strand, copper)))
                       for (k, n) in enumerate((6, 12, 18))]
            armour = Lay(Helix(LayRatio(10)),
                Repeat(CircularPattern(40; r = 26.5e-3),
                    Region(:armour_wires, Circle(2.0e-3), steel)))
            design = build(CableDesign, "coaxial",
                Terminal(:core,
                    Enclosure(Region(:core_fill, Circle(10.5e-3), air),
                        Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...))),
                Region(:core_semicon_inner, Shell(1.0e-3), semicon),
                Region(:core_insulation, Shell(8.0e-3), polymer),
                Region(:core_semicon_outer, Shell(1.0e-3), semicon),
                Terminal(:sheath, Region(:sheath_copper, Annulus(20.5e-3, 21.5e-3), copper)),
                Region(:sheath_bedding, Shell(3.0e-3), polymer),
                Terminal(:armour,
                    Enclosure(Region(:armour_fill, Shell(4.0e-3), air), armour)),
                Region(:armour_jacket, Shell(4.0e-3), polymer))
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
        end
    end

    @testset "a terminal around terminals, and one name for two parts" begin
        wire = Region(:core_metal, Circle(1.0e-3), copper)
        tube = Region(:core_tube, Shell(1.0e-3), copper)
        nested = build(CableDesign, "nested",
            Terminal(:phase, Stack(Terminal(:inner, wire), Terminal(:outer, tube))))
        @test [(row.name, row.regions) for row in nested.terminals] == [(:phase, [1, 2])]
        shared = build(CableDesign, "shared", Terminal(:core, wire),
            Region(:core_insulation, Shell(1.0e-3), polymer), Terminal(:core, tube))
        @test [(row.name, row.regions) for row in shared.terminals] == [(:core, [1, 3])]
    end

    @testset "each copy of a repeat is a scope of its own terminals" begin
        strands = Lay(Helix(LayRatio(12)),
            Repeat(CircularPattern(6; r = 1.0e-3), Region(:core_wires, Circle(0.5e-3), copper)))
        core = Stack(
            Terminal(:core,
                Enclosure(Region(:core_fill, Circle(1.5e-3), air),
                    Stack(Region(:core_centre, Circle(0.5e-3), copper), strands))),
            Region(:core_insulation, Shell(1.0e-3), polymer),
            Terminal(:screen, Region(:core_screen, Shell(0.5e-3), copper)))
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
        bonded = build(CableDesign, "bonded", Terminal(:phase, Enclosure(pipe, cores)))
        @test [(row.name, length(row.regions)) for row in bonded.terminals] ==
              [(:phase, 24)]
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
    end

    @testset "Float32 inputs" begin
        wire = build(CableDesign, "float32",
            Terminal(:core,
                Region(:core_metal, Circle(1.5f-3),
                    material(kind = :conductor, rho = 1.7241f-8))))
        @test r_ex(wire.outer) ≈ 1.5f-3
    end
end
