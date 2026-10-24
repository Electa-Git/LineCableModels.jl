@testset "DataModel systems" begin
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.3)
    air = material(kind = :insulator, rho = Inf)
    bare = build(CableDesign, "bare_wire",
        Terminal(:core, Region(:core_metal, Circle(0.0425), copper)))
    wire = build(CableDesign, "insulated_wire",
        Terminal(:core, Region(:core_metal, Circle(0.0425), copper)),
        Region(:core_insulation, Shell(1.0e-3), polymer))

    @testset "two wires side by side" begin
        for design in (bare, wire)
            system = build(CableSystem, [design, design],
                [Placement(0.0, -1.0), Placement(1.0, -1.0)];
                connections = [(core = 1,), (core = 2,)], id = "two_wires", line_length = 1.0)
            n = length(design.regions)
            @test length(system.regions) == 2n
            @test [(row.cable, row.terminal, row.connection, row.regions)
                   for row in system.terminals] == [(1, 1, 1, [1]), (2, 1, 2, [n + 1])]
            @test centroid(system.regions[n + 1].shape) == (1.0, -1.0)
        end
    end

    @testset "the coaxial cables of the references" begin
        semicon = Material(:semicon, 1000.0, 1000.0)
        steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
        for strand in (Circle(1.5e-3), Rectangle(3.0e-3, 3.0e-3))
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
            system = build(CableSystem, [design], [Placement(0.0, -1.0)];
                connections = [(core = 1, sheath = 2, armour = 3)], id = "coaxial",
                line_length = 1000.0)
            @test length(system.regions) == 85
            @test [(row.terminal, row.connection, length(row.regions))
                   for row in system.terminals] == [(1, 1, 37), (2, 2, 1), (3, 3, 40)]
            @test centroid(system.regions[1].shape) == (0.0, -1.0)
            x, y = centroid(system.regions[first(system.terminals[3].regions)].shape)
            @test hypot(x, y + 1.0) ≈ 26.5e-3
        end
    end

    @testset "touching cables keep their placements" begin
        R = r_ex(wire.outer)
        pair = [Placement(0.0, -1.0), Placement(2R, -1.0)]
        trefoil = [Placement(0.0, -1.0 + 2R / sqrt(3)), Placement(-R, -1.0 - R / sqrt(3)),
            Placement(R, -1.0 - R / sqrt(3))]
        for placements in (pair, trefoil)
            system = build(CableSystem, fill(wire, length(placements)), placements;
                connections = [(core = k,) for k in eachindex(placements)],
                id = "touching", line_length = 1.0)
            @test system.placements == placements
        end
        @test_throws DomainError build(CableSystem, [wire, wire],
            [Placement(0.0, -1.0), Placement(2R - 1.0e-6, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "overlap", line_length = 1.0)
    end

    @testset "a name connects every terminal of that name" begin
        strands = Lay(Helix(LayRatio(12)),
            Repeat(CircularPattern(6; r = 1.0e-3), Region(:core_wires, Circle(0.5e-3), copper)))
        core = Stack(
            Terminal(:core,
                Enclosure(Region(:core_fill, Circle(1.5e-3), air),
                    Stack(Region(:core_centre, Circle(0.5e-3), copper), strands))),
            Region(:core_insulation, Shell(1.0e-3), polymer),
            Terminal(:screen, Region(:core_screen, Shell(0.5e-3), copper)))
        pipe = build(CableDesign, "three_cores",
            Enclosure(Region(:pipe_fill, Circle(6.5e-3), polymer),
                Repeat(CircularPattern(3; r = 3.5e-3), core)))
        system = build(CableSystem, [pipe], [Placement(0.0, -1.0)];
            connections = [(core = 1, screen = 0)], id = "bundle", line_length = 1.0)
        @test [row.connection for row in system.terminals] == [1, 0, 1, 0, 1, 0]
    end

    @testset "the checks of a system" begin
        screened = build(CableDesign, "screened",
            Terminal(:core, Region(:core_metal, Circle(1.0e-3), copper)),
            Region(:core_insulation, Shell(1.0e-3), polymer),
            Terminal(:screen, Region(:core_screen, Shell(0.5e-3), copper)))
        at = [Placement(0.0, -1.0)]
        @test_throws DimensionMismatch build(CableSystem, [wire, wire],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,)], id = "count", line_length = 1.0)
        @test_throws ArgumentError build(CableSystem, [screened], at;
            connections = [(core = 1,)], id = "unnamed", line_length = 1.0)
        @test_throws ArgumentError build(CableSystem, [screened], at;
            connections = [(core = 1, screen = 0, armour = 0)], id = "unknown",
            line_length = 1.0)
        @test_throws DomainError build(CableSystem, [wire], at;
            connections = [(core = -1,)], id = "negative", line_length = 1.0)
        @test_throws DomainError build(CableSystem, [wire], at;
            connections = [(core = 1,)], id = "length", line_length = 0.0)
    end

    @testset "Float32 inputs" begin
        thin = build(CableDesign, "float32",
            Terminal(:core,
                Region(:core_metal, Circle(1.5f-3),
                    material(kind = :conductor, rho = 1.7241f-8))))
        system = build(CableSystem, [thin], [Placement(0.0f0, -1.0f0)];
            connections = [(core = 1,)], id = "float32", line_length = 1.0f0)
        @test system.line_length === 1.0f0
        @test system.placements == [Placement(0.0f0, -1.0f0)]
    end
end
