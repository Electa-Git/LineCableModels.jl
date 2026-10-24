@testset "DataModel blueprints" begin
    flatten = LineCableModels.flatten
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.3)
    semicon = Material(:semicon, 1000.0, 1000.0)
    steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
    air = material(kind = :insulator, rho = Inf)

    @testset "a wire is its own tube" begin
        wire = flatten(build(CableDesign, "insulated_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)),
            Region(:core_insulation, Shell(1.0e-3), polymer)))
        tube = only(wire.conductors)
        @test [tube.r_in, tube.r_ex] ≈ [0.0, 0.0425]
        @test tube.material.rho ≈ copper.rho && tube.material.mu_r ≈ 1
        @test [only(wire.dielectrics).r_in, only(wire.dielectrics).r_ex] ≈ [0.0425, 0.0435]
        @test wire.mu_r_ins == [1.0]
    end

    @testset "a wire inside a tube of its copper is a solid conductor" begin
        core = flatten(build(CableDesign, "wire_and_tube",
            Terminal(:core,
                Stack(Region(:core_wire, Circle(1.0e-3), copper),
                    Region(:core_tube, Shell(1.0e-3), copper)))))
        tube = only(core.conductors)
        @test [tube.r_in, tube.r_ex] ≈ [0.0, 2.0e-3]
        @test tube.material.rho ≈ copper.rho && tube.material.mu_r ≈ 1
    end

    @testset "the coaxial cables of the references" begin
        for (strand, core) in (
            (Circle(1.5e-3),
                (2.3580693965940444e-8, 1.0588802203492307, 1.416791518882097)),
            (Rectangle(3.0e-3, 3.0e-3),
                (1.8631839650067386e-8, 1.0072826882814461, 1.4073302667101928)))
            courses = [Lay(Helix(LayRatio(12)),
                           Repeat(CircularPattern(n; r = 3.0e-3k),
                               Region(Symbol(:core_strands_, k), strand, copper)))
                       for (k, n) in enumerate((6, 12, 18))]
            armour = Lay(Helix(LayRatio(10)),
                Repeat(CircularPattern(40; r = 26.5e-3),
                    Region(:armour_wires, Circle(2.0e-3), steel)))
            blueprint = flatten(build(CableDesign, "coaxial",
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
                Region(:armour_jacket, Shell(4.0e-3), polymer)))
            @test [tube.terminal for tube in blueprint.conductors] == [1, 2, 3]
            @test [tube.r_in for tube in blueprint.conductors] ≈ [0.0, 20.5e-3, 24.5e-3]
            @test [tube.r_ex for tube in blueprint.conductors] ≈ [10.5e-3, 21.5e-3, 28.5e-3]
            @test [layer.material for layer in blueprint.dielectrics] ==
                  [semicon, polymer, semicon, polymer, polymer]
            @test blueprint.dielectric_ranges == [1:3, 4:4, 5:5]
            values = [(tube.material.rho, tube.material.mu_r, mu)
                      for (tube, mu) in zip(blueprint.conductors, blueprint.mu_r_ins)]
            expected = [core, (copper.rho, 1.0, 1.0),
                (2.4270115498150639e-7, 5.8715848278702119, 1.1305521523424971)]
            @test all(isapprox(value, target; rtol = 1.0e-12)
            for (row, targets) in zip(values, expected)
            for (value, target) in zip(row, targets))
            @test [tube.material.alpha for tube in blueprint.conductors] ≈
                  [0.00393, 0.00393, 0.0045]
            @test blueprint.assembly_ranges == [1:3] && blueprint.positions == [(0.0, 0.0)]
        end
        # 37 round wires in parallel, the 36 laid ones longer by their overlength.
        wire = copper.rho / (π * 1.5e-3^2)
        resistance = wire / (1 + 36 / sqrt(1 + (π / 12)^2))
        @test resistance * π * 10.5e-3^2 ≈ 2.3580693965940444e-8
    end

    @testset "opposite lays cancel turns" begin
        relative = map((1, -1)) do dir
            courses = [Lay(Helix(LayRatio(12); dir = k == 2 ? dir : 1),
                           Repeat(CircularPattern(n; r = 3.0e-3k),
                               Region(Symbol(:core_strands_, k), Circle(1.5e-3), copper)))
                       for (k, n) in enumerate((6, 12))]
            core = flatten(build(CableDesign, "core",
                Terminal(:core,
                    Enclosure(Region(:core_fill, Circle(7.5e-3), air),
                        Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...))),
                Region(:core_insulation, Shell(8.0e-3), polymer)))
            return only(core.mu_r_ins)
        end
        @test relative[2] < relative[1]
    end

    @testset "designs without coaxial tubes" begin
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
        @test_throws ArgumentError flatten(build(CableDesign, "three_cores",
            Enclosure(pipe, cores)))
        @test_throws ArgumentError flatten(build(CableDesign, "bonded",
            Terminal(:phase, Enclosure(pipe, cores))))
        warm = Material(:conductor, 1.7241e-8, 1, 1, 25, 0.00393)
        @test_throws ArgumentError flatten(build(CableDesign, "two_temperatures",
            Terminal(:core,
                Stack(Region(:core_wire, Circle(1.0e-3), copper),
                    Region(:core_tube, Shell(1.0e-3), warm)))))
    end

    @testset "Float32 inputs" begin
        wire = flatten(build(CableDesign, "float32",
            Terminal(:core,
                Region(:core_metal, Circle(1.5f-3),
                    material(kind = :conductor, rho = 1.7241f-8)))))
        @test only(wire.conductors).r_ex === 1.5f-3
        @test wire.mu_r_ins == Float32[1]
    end
end
