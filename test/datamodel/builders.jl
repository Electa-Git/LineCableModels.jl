@testset "DataModel builders" begin
    flatten = LineCableModels.flatten
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.3)
    semicon = Material(:semicon, 1000.0, 1000.0)
    steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
    air = material(kind = :insulator, rho = Inf)

    @testset "the wires of the references" begin
        bare = @cable "bare_wire" begin
            @terminal :core begin
                core(copper; r = 0.0425)
            end
        end
        insulated = @cable "insulated_wire" begin
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
            design = @cable "coaxial" begin
                @terminal :core begin
                    stranded(copper; center = Circle(1.5e-3), shape = strand,
                        n = (6, 12, 18), lay = LayRatio(12), boundary = Circle(10.5e-3),
                        fill = air)
                end
                screen(semicon; t = 1.0e-3)
                insulation(polymer; t = 8.0e-3)
                screen(semicon; t = 1.0e-3)
                @terminal :sheath begin
                    sheath(copper; t = 1.0e-3)
                end
                bedding(polymer; t = 3.0e-3)
                @terminal :armour begin
                    armor(steel; shape = Circle(2.0e-3), n = 40, lay = LayRatio(10),
                        fill = air)
                end
                jacket(polymer; t = 4.0e-3)
            end
            # The node tree of the reference, with the radii of its courses and its sheath.
            courses = [Lay(Helix(LayRatio(12)),
                           Repeat(CircularPattern(n; r = 3.0e-3k),
                               Region(Symbol(:core_strands_, k), strand, copper)))
                       for (k, n) in enumerate((6, 12, 18))]
            armour = Lay(Helix(LayRatio(10)),
                Repeat(CircularPattern(40; r = 26.5e-3),
                    Region(:armour_wires, Circle(2.0e-3), steel)))
            reference = build(CableDesign, "coaxial",
                Terminal(:core,
                    Enclosure(Region(:core_fill, Circle(10.5e-3), air),
                        Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...))),
                Region(:core_semicon_inner, Offset(1.0e-3), semicon),
                Region(:core_insulation, Offset(8.0e-3), polymer),
                Region(:core_semicon_outer, Offset(1.0e-3), semicon),
                Terminal(:sheath, Region(:sheath_copper, Annulus(20.5e-3, 21.5e-3), copper)),
                Region(:sheath_bedding, Offset(3.0e-3), polymer),
                Terminal(:armour,
                    Enclosure(Region(:armour_fill, Offset(4.0e-3), air), armour)),
                Region(:armour_jacket, Offset(4.0e-3), polymer))

            @test [placed.source.material for placed in design.regions] ==
                  [placed.source.material for placed in reference.regions]
            @test all(isapprox(f(a.shape), f(b.shape); rtol = 1.0e-12, atol = 1.0e-15)
            for (a, b) in zip(design.regions, reference.regions)
            for f in (area, support, distance, shape -> first(centroid(shape)),
                shape -> last(centroid(shape))))
            @test [design.regions[row.region].source.tag for row in design.fills] ==
                  [:stranded_fill, :armor_fill]
            @test [(row.region, row.contents) for row in design.fills] ==
                  [(row.region, row.contents) for row in reference.fills]
            @test [(row.name, row.scope, row.regions) for row in design.terminals] ==
                  [(row.name, row.scope, row.regions) for row in reference.terminals]
            @test [row.copies for row in design.repeats] ==
                  [row.copies for row in reference.repeats]
            @test all(isapprox(f(a), f(b))
            for (a, b) in zip(design.lays, reference.lays)
            for f in (row -> row.radii, row -> row.pitches, row -> row.overlengths))

            blueprint, expected = flatten(design), flatten(reference)
            @test [tube.terminal for tube in blueprint.conductors] ==
                  [tube.terminal for tube in expected.conductors]
            @test all(isapprox(f(a), f(b); rtol = 1.0e-12, atol = 1.0e-15)
            for (a, b) in zip(blueprint.conductors, expected.conductors)
            for f in (tube -> tube.r_in, tube -> tube.r_ex, tube -> tube.material.rho,
                tube -> tube.material.mu_r, tube -> tube.material.alpha))
            @test all(isapprox(f(a), f(b); rtol = 1.0e-12)
            for (a, b) in zip(blueprint.dielectrics, expected.dielectrics)
            for f in (layer -> layer.r_in, layer -> layer.r_ex))
            @test [layer.material for layer in blueprint.dielectrics] ==
                  [layer.material for layer in expected.dielectrics]
            @test blueprint.dielectric_ranges == expected.dielectric_ranges
            @test isapprox(blueprint.mu_r_ins, expected.mu_r_ins; rtol = 1.0e-12)
        end
    end

    @testset "Float32 inputs and the default fill" begin
        wire = material(kind = :conductor, rho = 1.7241f-8)
        design = @cable "float32" begin
            @terminal :core begin
                stranded(wire; center = Circle(1.5f-3), shape = Circle(1.5f-3), n = (6,),
                    lay = LayRatio(12.0f0), boundary = Circle(4.5f-3))
            end
        end
        fill = last(design.regions).source
        @test fill.tag == :stranded_fill && fill.material.rho == Inf
        @test r_ex(design.outer) ≈ 4.5f-3
    end
end
