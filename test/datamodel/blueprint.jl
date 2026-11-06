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
            Region(:core_insulation, Offset(1.0e-3), polymer)))
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
                    Region(:core_tube, Offset(1.0e-3), copper)))))
        tube = only(core.conductors)
        @test [tube.r_in, tube.r_ex] ≈ [0.0, 2.0e-3]
        @test tube.material.rho ≈ copper.rho && tube.material.mu_r ≈ 1
    end

    @testset "a lay over a strand at the radius 0 keeps it straight" begin
        strand = Repeat(CircularPattern(1; r = 0),
            Region(:core_metal, Circle(0.0425), copper))
        laid = build(CableDesign, "laid_wire",
            Terminal(:core, Lay(Helix(LayRatio(12)), strand)),
            Region(:core_insulation, Offset(1.0e-3), polymer))
        straight = build(CableDesign, "straight_wire", Terminal(:core, strand),
            Region(:core_insulation, Offset(1.0e-3), polymer))
        @test only(laid.lays).pitches == [Inf]
        blueprint, expected = flatten(laid), flatten(straight)
        @test only(blueprint.conductors).material.rho ==
              only(expected.conductors).material.rho
        @test blueprint.mu_r_ins == expected.mu_r_ins == [1.0]
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
                Region(:core_semicon_inner, Offset(1.0e-3), semicon),
                Region(:core_insulation, Offset(8.0e-3), polymer),
                Region(:core_semicon_outer, Offset(1.0e-3), semicon),
                Terminal(:sheath, Region(:sheath_copper, Annulus(20.5e-3, 21.5e-3), copper)),
                Region(:sheath_bedding, Offset(3.0e-3), polymer),
                Terminal(:armour,
                    Enclosure(Region(:armour_fill, Offset(4.0e-3), air), armour)),
                Region(:armour_jacket, Offset(4.0e-3), polymer)))
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
            @test isempty(only(blueprint.interior_ranges))
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
                Region(:core_insulation, Offset(8.0e-3), polymer)))
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
            Region(:core_insulation, Offset(1.0e-3), polymer),
            Terminal(:screen, Region(:core_screen, Offset(0.5e-3), copper)))
        cores = Repeat(CircularPattern(3; r = 3.5e-3), core)
        pipe = Region(:pipe_fill, Circle(6.5e-3), polymer)
        three = flatten(build(CableDesign, "three_cores", Enclosure(pipe, cores)))
        @test [tube.assembly for tube in three.conductors] == [2, 2, 3, 3, 4, 4]
        @test length(three.positions) == 4 && three.positions[1] == (0.0, 0.0)
        @test [layer.material for layer in three.dielectrics[three.interior_ranges[1]]] ==
              [polymer]
        @test_throws ArgumentError flatten(build(CableDesign, "bonded",
            Terminal(:phase, Enclosure(pipe, cores))))
        warm = Material(:conductor, 1.7241e-8, 1, 1, 25, 0.00393)
        @test_throws ArgumentError flatten(build(CableDesign, "two_temperatures",
            Terminal(:core,
                Stack(Region(:core_wire, Circle(1.0e-3), copper),
                    Region(:core_tube, Offset(1.0e-3), warm)))))
    end

    # A block that tiles no annulus about the center is the tube of its area and perimeter
    # (Ametani). The layer around it takes the thickness of its band.
    @testset "a sector core is the tube of its area and perimeter" begin
        sector = LineCableModels.Sector(
            0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(0, 0))
        blueprint = flatten(build(CableDesign, "sector_core",
            Terminal(:core, Region(:core, sector, copper)),
            Region(:insulation, Offset(1.0e-3), polymer)))
        tube = only(blueprint.conductors)
        @test [tube.r_in, tube.r_ex] ≈ [2.3105861599593243e-3, 5.652135443739411e-3] rtol = 1e-12
        @test π * (tube.r_ex^2 - tube.r_in^2) ≈ 83.59095280005376e-6 rtol = 1e-12
        @test 2π * tube.r_ex ≈ 35.51341437429244e-3 rtol = 1e-12
        @test tube.material.rho ≈ copper.rho rtol = 1e-12
        layer = only(blueprint.dielectrics)
        @test [layer.r_in, layer.r_ex] ≈ [tube.r_ex, tube.r_ex + 1.0e-3] rtol = 1e-12
    end

    # An open tape and its fill tile one layer. Under a radial field the two carry their
    # currents in parallel, and the admittivity of the layer is the area-weighted sum of theirs.
    @testset "an open dielectric tape is one layer of the composite material" begin
        polyester = material(kind = :insulator, rho = 1.0e15, eps_r = 3.3, tan_delta = 2.0e-3)
        compound = material(kind = :insulator, rho = 1.0e13, eps_r = 2.2, tan_delta = 5.0e-4)
        blueprint = flatten(build(CableDesign, "open_tape",
            Terminal(:core, Region(:core, Circle(10.0e-3), copper)),
            tape(polyester; w = 10.0e-3, t = 0.2e-3, lay = LayRatio(10), fill = compound),
            Region(:insulation, Offset(1.0e-3), polymer)))
        @test length(blueprint.dielectrics) == 2
        layer = first(blueprint.dielectrics)
        f = 10.0e-3 * 0.2e-3 / (π * (10.2e-3^2 - 10.0e-3^2))
        @test layer.material isa Material{Float64, :insulator}
        @test [layer.r_in, layer.r_ex] ≈ [10.0e-3, 10.2e-3] rtol = 1e-12
        @test 1 / layer.material.rho ≈ f / 1.0e15 + (1 - f) / 1.0e13 rtol = 1e-12
        @test layer.material.eps_r ≈ f * 3.3 + (1 - f) * 2.2 rtol = 1e-12
        @test layer.material.eps_r * layer.material.tan_delta ≈
              f * 3.3 * 2.0e-3 + (1 - f) * 2.2 * 5.0e-4 rtol = 1e-12
        for ω in 2π .* (50.0, 1.0e3, 1.0e6)
            @test 1 / layer.material.rho +
                  ω * LineCableModels.ε₀ * layer.material.eps_r *
                  (layer.material.tan_delta + im) ≈
                  f * (1 / 1.0e15 + ω * LineCableModels.ε₀ * 3.3 * (2.0e-3 + im)) +
                  (1 - f) * (1 / 1.0e13 + ω * LineCableModels.ε₀ * 2.2 * (5.0e-4 + im)) rtol = 1e-12
        end
    end

    # Construction 10. The cores lie at their centroids. The root has no conductor, so its fill
    # and its common jacket form its interior.
    @testset "the cores of a three-core sector cable sit at their centroids" begin
        sector = LineCableModels.Sector(
            0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(0, 0))
        c = first(LineCableModels.centroid(sector))
        core = LineCableModels.Sector(
            0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(-c, 0))
        angles = (0, 2π / 3, 4π / 3)
        blueprint = flatten(build(CableDesign, "three_sector_cores",
            Stack(
                Enclosure(Region(:fill, Circle(11.0e-3), air),
                    Assembly(Tuple(Placement(c * cos(φ),
                                       c * sin(φ),
                                       φ) => Stack(
                                       Terminal(name, Region(:core, core, copper)),
                                       Region(:insulation, Offset(1.0e-3), polymer))
                    for (φ, name) in zip(angles, (:phase_a, :phase_b, :phase_c))))),
                Region(:jacket, Offset(3.0e-3), polymer))))
        @test [tube.assembly for tube in blueprint.conductors] == [2, 3, 4]
        @test all(isapprox([tube.r_in, tube.r_ex],
                      [2.3105861599593243e-3, 5.652135443739411e-3]; rtol = 1e-12)
        for tube in blueprint.conductors)
        @test [collect(p) for p in blueprint.positions] ≈
              vcat([[0.0, 0.0]], [[c * cos(φ), c * sin(φ)] for φ in angles]) rtol = 1e-12
        @test isempty(blueprint.assembly_ranges[1])
        @test blueprint.assembly_ranges[2:4] == [1:1, 2:2, 3:3]
        @test blueprint.dielectric_ranges == [1:1, 2:2, 3:3]
        @test blueprint.interior_ranges[1] == 4:5 &&
              all(isempty, blueprint.interior_ranges[2:4])
        @test [layer.material for layer in blueprint.dielectrics] ==
              [polymer, polymer, polymer, air, polymer]
        @test all(isapprox(layer.r_ex, layer.r_in + 1.0e-3; rtol = 1e-12)
        for layer in blueprint.dielectrics[1:3])
    end

    # Construction 15: three coaxials in a trefoil, each on the wall of a steel pipe. The
    # members' conductors precede the pipe's, and the fill is the interior of the root.
    @testset "three coaxials in a pipe flatten to their rows" begin
        d = 2 * 22.0e-3 / sqrt(3)
        θs = (π / 2, 7π / 6, 11π / 6)
        coax = Stack(Terminal(:core, Region(:core, Circle(10.0e-3), copper)),
            Region(:insulation, Offset(8.0e-3), polymer),
            Terminal(:screen, Region(:screen, Offset(1.0e-3), copper)),
            Region(:jacket, Offset(3.0e-3), polymer))
        design = build(CableDesign, "pipe",
            Stack(
                Enclosure(Region(:fill, Circle(22.0e-3 * (1 + 2 / sqrt(3))), air),
                    Assembly(Tuple(Placement(d * cos(θ), d * sin(θ)) => coax for θ in θs))),
                Terminal(:pipe, Region(:pipe, Offset(8.0e-3), steel)),
                Region(:coating, Offset(3.0e-3), polymer)))
        blueprint = flatten(design)
        R = 0.04740341184434353
        @test [(tube.terminal, tube.assembly) for tube in blueprint.conductors] ==
              [(1, 2), (2, 2), (3, 3), (4, 3), (5, 4), (6, 4), (7, 1)]
        @test [[tube.r_in, tube.r_ex] for tube in blueprint.conductors] ≈
              [[0, 0.010], [0.018, 0.019], [0, 0.010], [0.018, 0.019], [0, 0.010],
            [0.018, 0.019], [R, R + 0.008]] rtol = 1e-12
        @test [[layer.r_in, layer.r_ex] for layer in blueprint.dielectrics] ≈
              [[0.010, 0.018], [0.019, 0.022], [0.010, 0.018], [0.019, 0.022],
            [0.010, 0.018], [0.019, 0.022], [0, R], [R + 0.008, R + 0.011]] rtol = 1e-12
        @test [layer.material for layer in blueprint.dielectrics] ==
              [polymer, polymer, polymer, polymer, polymer, polymer, air, polymer]
        @test blueprint.dielectric_ranges == [1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 8:8]
        @test blueprint.assembly_ranges == [7:7, 1:2, 3:4, 5:6]
        @test blueprint.interior_ranges == [7:7, 1:0, 3:2, 5:4]
        @test [collect(p) for p in blueprint.positions] ≈
              vcat([[0.0, 0.0]], [[d * cos(θ), d * sin(θ)] for θ in θs]) rtol = 1e-12
        @test blueprint.mu_r_ins ≈ ones(7) rtol = 1e-12
        system = build(CableSystem, [design], [Placement(0, -1)];
            connections = [(core = 1, screen = 0, pipe = 0)], id = "pipe", line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
            frequencies = [50.0])
        @test_throws ArgumentError(
            "the coaxial formulation computes one concentric assembly per cable; the cable pipe holds 4 assemblies") compute(
            problem, CoaxialFormulation())
    end
end
