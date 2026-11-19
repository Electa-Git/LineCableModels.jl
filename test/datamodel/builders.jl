@testset "DataModel builders" begin
    flatten = LineCableModels.flatten
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.3)
    semicon = Material(:semicon, 1000.0, 1000.0)
    steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
    air = material(kind = :insulator, rho = Inf)

    @testset "the wires of the references" begin
        bare = @design "bare_wire" begin
            @terminal :core begin
                core(copper; r = 0.0425)
            end
        end
        insulated = @design "insulated_wire" begin
            @terminal :core begin
                core(copper; r = 0.0425)
            end
            insulation(polymer; t = 1.0e-3)
        end
        @test [(row.name, row.regions) for row in bare.terminals] == [(:core, [1])]
        @test r_ex(bare.outer) == 0.0425
        @test [(placed.source.tag, placed.source.material) for placed in insulated.regions] ==
              [(:core, copper), (:insulation, polymer)]
        @test [(row.name, row.regions) for row in insulated.terminals] == [(:core, [1])]
        @test r_ex(insulated.outer) ≈ 0.0435
    end

    @testset "the coaxial cables of the references" begin
        for strand in (Circle(1.5e-3), Rectangle(3.0e-3, 3.0e-3))
            design = @design "coaxial" begin
                @terminal :core begin
                    core(copper; r = 1.5e-3)
                    stranded(copper; shape = strand, n = (6, 12, 18), lay = LayRatio(12),
                        boundary = Annulus(1.5e-3, 10.5e-3), fill = air)
                end
                screen(semicon; t = 1.0e-3)
                insulation(polymer; t = 8.0e-3)
                screen(semicon; t = 1.0e-3)
                @terminal :sheath begin
                    sheath(copper; t = 1.0e-3)
                end
                bedding(polymer; t = 3.0e-3)
                @terminal :armour begin
                    wires(steel; shape = Circle(2.0e-3), n = 40, lay = LayRatio(10),
                        fill = air, tag = :armor)
                end
                jacket(polymer; t = 4.0e-3)
            end
            # The node tree of the reference, with the radii of its courses and its sheath,
            # and its central wire a part before its courses.
            courses = [Lay(Helix(LayRatio(12)),
                           Repeat(CircularPattern(n; r = 3.0e-3k),
                               Region(Symbol(:core_strands_, k), strand, copper)))
                       for (k, n) in enumerate((6, 12, 18))]
            armour = Lay(Helix(LayRatio(10)),
                Repeat(CircularPattern(40; r = 26.5e-3),
                    Region(:armour_wires, Circle(2.0e-3), steel)))
            reference = build(CableDesign, "coaxial",
                Terminal(:core,
                    Stack(Region(:core_centre, Circle(1.5e-3), copper),
                        Group(Region(:core_fill, Annulus(1.5e-3, 10.5e-3), air),
                            Stack(courses...)))),
                Region(:core_semicon_inner, Offset(1.0e-3), semicon),
                Region(:core_insulation, Offset(8.0e-3), polymer),
                Region(:core_semicon_outer, Offset(1.0e-3), semicon),
                Terminal(:sheath, Region(:sheath_copper, Annulus(20.5e-3, 21.5e-3), copper)),
                Region(:sheath_bedding, Offset(3.0e-3), polymer),
                Terminal(:armour,
                    Group(Region(:armour_fill, Offset(4.0e-3), air), armour)),
                Region(:armour_jacket, Offset(4.0e-3), polymer))

            @test [placed.source.material for placed in design.regions] ==
                  [placed.source.material for placed in reference.regions]
            @test all(isapprox(f(a.shape), f(b.shape); rtol = 1.0e-12, atol = 1.0e-15)
            for (a, b) in zip(design.regions, reference.regions)
            for f in (area, support, distance, shape -> first(centroid(shape)),
                shape -> last(centroid(shape))))
            @test [placed.source.tag
                   for placed in design.regions
                   if placed.shape isa LineCableModels.PerforatedShape] ==
                  [:stranded_fill, :armor_fill]
            @test [(row.name, row.scope, row.regions) for row in design.terminals] ==
                  [(row.name, row.scope, row.regions) for row in reference.terminals]
            @test [row.copies for row in design.repeats] ==
                  [row.copies for row in reference.repeats]

            input, expected = coaxial_input(design), coaxial_input(reference)
            @test input.label == expected.label
            @test all(isapprox(a, b; rtol = 1.0e-12, atol = 1.0e-15)
            for f in (input -> input.r_in, input -> input.r_ex)
            for (a, b) in zip(f(input), f(expected)))
            @test all(isapprox(f(a), f(b); rtol = 1.0e-12, atol = 1.0e-15)
            for (a, b) in zip(input.material, expected.material)
            for f in (tube -> tube.rho, tube -> tube.mu_r, tube -> tube.alpha))
            @test all(isapprox(a, b; rtol = 1.0e-12)
            for f in (input -> input.dielectric_r_in, input -> input.dielectric_r_ex)
            for (a, b) in zip(f(input), f(expected)))
            @test input.dielectric_material == expected.dielectric_material
            @test input.dielectric_ranges == expected.dielectric_ranges
            @test isapprox(input.mu_r_ins, expected.mu_r_ins; rtol = 1.0e-12)
        end
    end

    @testset "a tape, a course of wires and the formations of cables" begin
        @test tape(copper; w = 10.0e-3, t = 0.1e-3, lay = LayRatio(10)) ==
              Group(Region(:tape_fill, Offset(0.1e-3), Material(:insulator, Inf)),
            Lay(Helix(LayRatio(10)),
                Repeat(CircularPattern(1; r = Offset(0)),
                    Region(:tape, Rectangle(10.0e-3, 0.1e-3), copper))))
        @test wires(copper; shape = Circle(1.0e-3), n = 12, lay = LayRatio(8)) ==
              Group(Region(:wires_fill, Offset(2.0e-3), Material(:insulator, Inf)),
            Lay(Helix(LayRatio(8)),
                Repeat(CircularPattern(12; r = Offset(0)),
                    Region(:wires, Circle(1.0e-3), copper))))
        wire = @design "insulated_wire" begin
            @terminal :core begin
                core(copper; r = 0.0425)
            end
            insulation(polymer; t = 1.0e-3)
        end
        R = r_ex(wire.outer)
        formation = trefoil(wire; spacing = 2R, at = Placement(0.0, -1.0))
        touching = [Placement(0.0, -1.0 + 2R / sqrt(3)), Placement(-R, -1.0 - R / sqrt(3)),
            Placement(R, -1.0 - R / sqrt(3))]
        @test isapprox([c for p in formation for c in (p.x, p.y, p.φ)],
            [c for p in touching for c in (p.x, p.y, p.φ)]; rtol = 1.0e-15)
        @test_throws DomainError trefoil(wire; spacing = 2R - 1.0e-6)
        @test isapprox(
            [c
             for p in flat(wire; n = 3, spacing = 0.1, at = Placement(0.0, -1.0))
             for c in (p.x, p.y, p.φ)],
            [-0.1, -1.0, 0.0, 0.0, -1.0, 0.0, 0.1, -1.0, 0.0]; rtol = 1.0e-15)
        @test isapprox(
            [c
             for p in
                 flat(wire; n = 3, spacing = 0.1, φ = π / 2, at = Placement(0.0, -1.0))
             for c in (p.x, p.y, p.φ)],
            [0.0, -1.1, 0.0, 0.0, -1.0, 0.0, 0.0, -0.9, 0.0]; rtol = 1.0e-15)
        @test_throws DomainError flat(wire; n = 3, spacing = 2R - 1.0e-6)
    end

    @testset "the courses of a stranded conductor lie on the part inside it" begin
        oil = Material(:insulator, 1.0e12, 2.2)
        wire = Circle(1.5e-3)
        # A course of one wire with no part inside it is the wire on the center. The same
        # conductor with its central wire declared before its courses has the same coaxial
        # input.
        centered = @design "centered" begin
            @terminal :core begin
                stranded(copper; shape = wire, n = (1, 6, 12, 18), lay = LayRatio(12),
                    boundary = Circle(10.5e-3), fill = air)
                insulation(polymer; t = 8.0e-3)
            end
        end
        @test centroid(first(centered.regions).shape) == (0.0, 0.0)
        composed = @design "composed" begin
            @terminal :core begin
                core(copper; r = 1.5e-3)
                stranded(copper; shape = wire, n = (6, 12, 18), lay = LayRatio(12),
                    boundary = Annulus(1.5e-3, 10.5e-3), fill = air)
                insulation(polymer; t = 8.0e-3)
            end
        end
        input, expected = coaxial_input(centered), coaxial_input(composed)
        @test all(isapprox(only(f(input)), only(f(expected)); rtol = 1.0e-12)
        for f in (input -> input.r_ex, input -> first(input.dielectric_r_ex),
            input -> only(input.material).rho, input -> only(input.material).mu_r,
            input -> input.mu_r_ins))
        # Courses around a filler have no wire on the center.
        hollow = @design "hollow" begin
            Region(:duct, Circle(6.0e-3), oil)
            @terminal :core begin
                stranded(copper; shape = wire, n = (12, 18), lay = LayRatio(12),
                    boundary = Annulus(6.0e-3, 12.0e-3), fill = oil)
            end
            insulation(polymer; t = 8.0e-3)
        end
        input = coaxial_input(hollow)
        @test [only(input.r_in), only(input.r_ex)] ≈ [6.0e-3, 12.0e-3]
        # The wires of a first course of three touch each other about the center.
        r = 1.5e-3 / sinpi(1 / 3)
        three = build(CableDesign, "three",
            Terminal(:core,
                stranded(copper; shape = wire, n = (3, 9), lay = LayRatio(12),
                    boundary = Circle(r + 4.5e-3), fill = air)))
        @test all(placed -> hypot(centroid(placed.shape)...) ≈ r, three.regions[1:3])
        # A rectangular strand bends around the part inside its course.
        @test_throws DomainError build(CableDesign, "rods",
            Terminal(:core,
                stranded(copper; shape = Rectangle(3.0e-3, 3.0e-3), n = (6,),
                    lay = LayRatio(12), boundary = Circle(4.5e-3))))
    end

    @testset "a tape alone and the courses of a stranded core" begin
        taped = @design "taped_wire" begin
            @terminal :core begin
                core(copper; r = 1.0e-3)
            end
            insulation(polymer; t = 1.0e-3)
            @terminal :screen begin
                tape(copper; w = 2.0e-3, t = 0.1e-3, lay = LayRatio(10))
            end
        end
        @test [placed.source.tag
               for placed in taped.regions
               if placed.shape isa LineCableModels.PerforatedShape] == [:tape_fill]
        input = coaxial_input(taped)
        @test [last(input.r_in), last(input.r_ex)] ≈ [2.0e-3, 2.1e-3]
        conductor = stranded(copper; shape = Circle(1.5e-3), n = (1, 6, 12),
            lay = (LayRatio(10), LayRatio(10), LayRatio(12)), boundary = Circle(7.5e-3),
            φ0 = (0, 0, π / 12))
        @test (conductor.fill, conductor.contents.items) ==
              (Region(:stranded_fill, Circle(7.5e-3), Material(:insulator, Inf)),
            [
                Lay(Helix(LayRatio(10)),
                    Repeat(CircularPattern(1; r = Offset(0)),
                        Region(:stranded, Circle(1.5e-3), copper))),
                Lay(Helix(LayRatio(10)),
                    Repeat(CircularPattern(6; r = Offset(0)),
                        Region(:stranded, Circle(1.5e-3), copper))),
                Lay(Helix(LayRatio(12)),
                    Repeat(CircularPattern(12; r = Offset(0), φ0 = π / 12),
                        Region(:stranded, Circle(1.5e-3), copper)))])
    end
end
