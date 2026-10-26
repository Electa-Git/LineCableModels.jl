@testset "Engine numerical input" begin
    CoaxialInput = LineCableModels.Engine.CoaxialInput
    FormulationOptions = LineCableModels.FormulationOptions
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.3)
    semicon = Material(:semicon, 1000.0, 1000.0)
    steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
    air = material(kind = :insulator, rho = Inf)
    reductions = FormulationOptions(CoaxialFormulation)

    @testset "the wires of the references" begin
        metal = Terminal(:core, Region(:core_metal, Circle(0.0425), copper))
        for (design, outer) in (
            (build(CableDesign, "bare_wire", metal), 0.0425),
            (
                build(CableDesign, "insulated_wire", metal,
                    Region(:core_insulation, Offset(1.0e-3), polymer)),
                0.0435))
            system = build(CableSystem, [design, design],
                [Placement(0.0, -1.0), Placement(1.0, -1.0)];
                connections = [(core = 1,), (core = 2,)], id = design.id, line_length = 1.0)
            problem = LineParametersProblem(system, homogeneous(rho = 0.1); frequencies = [50.0])
            input = CoaxialInput(problem, reductions)
            @test input.cable == [1, 2] && input.label == ["cable:1:core", "cable:2:core"]
            @test input.r_in == [0.0, 0.0] && input.r_ex == [0.0425, 0.0425]
            @test [tube.rho for tube in input.material] ≈ [copper.rho, copper.rho]
            @test [tube.mu_r for tube in input.material] ≈ [1, 1] &&
                  input.mu_r_ins == [1, 1]
            @test input.r_ins_ex ≈ [outer, outer] && input.radius ≈ [outer, outer]
            @test input.x == [0.0, 1.0] && input.y == [-1.0, -1.0]
            @test input.assemblies == [1:1, 2:2]
            @test input.permutation == [1, 2] && input.kept == [1, 2]
            @test isempty(input.bundles) && isempty(input.eliminated)
            @test input.terminals == input.label
        end
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
            design = build(CableDesign, "coaxial",
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
            system = build(CableSystem, [design], [Placement(0.0, -1.0)];
                connections = [(core = 1, sheath = 2, armour = 3)], id = "coaxial",
                line_length = 1.0)
            problem = LineParametersProblem(system, homogeneous(rho = 100, eps_r = 10);
                frequencies = [50.0])
            input = CoaxialInput(problem, reductions)
            @test input.r_in ≈ [0.0, 20.5e-3, 24.5e-3]
            @test input.r_ex ≈ [10.5e-3, 21.5e-3, 28.5e-3]
            values = [(tube.rho, tube.mu_r, mu)
                      for (tube, mu) in zip(input.material, input.mu_r_ins)]
            expected = [core, (copper.rho, 1.0, 1.0),
                (2.4270115498150639e-7, 5.8715848278702119, 1.1305521523424971)]
            @test all(isapprox(value, target; rtol = 1.0e-12)
            for (row, targets) in zip(values, expected)
            for (value, target) in zip(row, targets))
            @test input.r_ins_ex ≈ [20.5e-3, 24.5e-3, 32.5e-3]
            @test input.layers == [1:3, 4:4, 5:5]
            @test input.layer_material == [semicon, polymer, semicon, polymer, polymer]
            @test input.layer_r_in ≈ [10.5e-3, 11.5e-3, 19.5e-3, 21.5e-3, 28.5e-3]
            @test input.layer_r_ex ≈ [11.5e-3, 19.5e-3, 20.5e-3, 24.5e-3, 32.5e-3]
            @test input.x == [0.0] && input.y == [-1.0]
            @test input.radius == [support(design.outer)] && input.assemblies == [1:3]
            @test input.connection == [1, 2, 3] && input.permutation == [1, 2, 3]
            @test input.kept == [1, 2, 3] && isempty(input.eliminated)
            @test input.terminals == ["cable:1:core", "cable:1:sheath", "cable:1:armour"]
        end
    end

    @testset "the reduction plan" begin
        design = build(CableDesign, "three_tubes",
            Terminal(:core, Region(:core_metal, Circle(10.0e-3), copper)),
            Region(:core_insulation, Offset(10.0e-3), polymer),
            Terminal(:sheath, Region(:sheath_copper, Offset(1.0e-3), copper)),
            Region(:sheath_bedding, Offset(3.0e-3), polymer),
            Terminal(:armour, Region(:armour_steel, Offset(4.0e-3), steel)))
        earth = homogeneous(rho = 100, eps_r = 10)
        bonded = build(CableSystem, [design], [Placement(0.0, -1.0)];
            connections = [(core = 1, sheath = 0, armour = 1)], id = "bonded",
            line_length = 1.0)
        problem = LineParametersProblem(bonded, earth; frequencies = [50.0])
        joined = "cable:1:core+cable:1:armour"
        for (bundle, kron, bundles, kept, terminals) in (
            (true, true, [(1, 2)], [1], [joined]),
            (true, false, [(1, 2)], [1, 3], [joined, "cable:1:sheath"]),
            (false, true, Tuple{Int, Int}[], [1, 2], ["cable:1:core", "cable:1:armour"]),
            (false, false, Tuple{Int, Int}[], [1, 2, 3],
                ["cable:1:core", "cable:1:armour", "cable:1:sheath"]))
            options = FormulationOptions(CoaxialFormulation; bundle, kron)
            input = CoaxialInput(problem, options)
            @test input.permutation == [1, 3, 2]
            @test input.bundles == bundles && input.kept == kept
            @test input.eliminated == setdiff(1:3, kept) && input.terminals == terminals
        end
        two = build(
            CableSystem, [design, design], [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [
                (core = 1, sheath = 0, armour = 0), (core = 2, sheath = 0, armour = 0)],
            id = "two_cables", line_length = 1.0)
        input = CoaxialInput(LineParametersProblem(two, earth; frequencies = [50.0]), reductions)
        @test input.permutation == [1, 4, 2, 3, 5, 6] && input.kept == [1, 2]
        @test input.terminals == ["cable:1:core", "cable:2:core"]
        @test input.layers == [1:1, 2:2, 3:2, 3:3, 4:4, 5:4]
        @test input.assemblies == [1:3, 4:6]
    end

    @testset "a Float32 design arrives in Float64" begin
        copper32 = Material(:conductor, 1.7241f-8, 1, 1, 20, 0.00393f0)
        polymer32 = Material(:insulator, 1.0f14, 2.3f0)
        wire = build(CableDesign, "insulated_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425f0), copper32)),
            Region(:core_insulation, Offset(1.0f-3), polymer32))
        system = build(CableSystem, [wire], [Placement(0.0f0, -1.0f0)];
            connections = [(core = 1,)], id = "float32", line_length = 1.0f0)
        problem = LineParametersProblem(system, homogeneous(rho = 100.0f0);
            temperature = 20.0f0, frequencies = Float32[50])
        input = CoaxialInput(problem, reductions)
        @test input isa CoaxialInput{Float64}
        @test input.material isa Vector{Material{Float64, :conductor}}
        @test only(input.layer_material) isa Material{Float64, :insulator}
        @test only(input.layer_material).eps_r == 2.3f0 && input.r_ex == [0.0425f0]
    end
end
