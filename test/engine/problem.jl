@testset "Engine problems" begin
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    wire = build(CableDesign, "bare_wire",
        Terminal(:core, Region(:core_metal, Circle(0.0425), copper)))
    earth = homogeneous(rho = 100, eps_r = 10)
    system = build(CableSystem, [wire, wire], [Placement(0.0, -1.0), Placement(1.0, -1.0)];
        connections = [(core = 1,), (core = 2,)], id = "two_wires", line_length = 1.0)

    @testset "the values take the type the inputs promote to" begin
        problem = LineParametersProblem(system, earth; frequencies = [50, 1000])
        @test problem isa LineParametersProblem{Float64}
        @test problem.temperature === 20.0
        copper32 = Material(:conductor, 1.7241f-8, 1, 1, 20, 0.00393f0)
        wire32 = build(CableDesign, "bare_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425f0), copper32)))
        system32 = build(CableSystem, [wire32, wire32],
            [Placement(0.0f0, -1.0f0), Placement(1.0f0, -1.0f0)];
            connections = [(core = 1,), (core = 2,)], id = "two_wires", line_length = 1.0f0)
        earth32 = homogeneous(rho = 100.0f0, eps_r = 10.0f0)
        frequencies = Float32[50, 1000]
        @test LineParametersProblem(system32, earth32; temperature = 20.0f0, frequencies) isa
              LineParametersProblem{Float32}
        @test LineParametersProblem(system32, earth32; temperature = 20, frequencies) isa
              LineParametersProblem{Float64}
    end

    @testset "the frequencies increase, positive and finite, at a finite temperature" begin
        for frequencies in
            (Float64[], [0.0, 50.0], [50.0, Inf], [1000.0, 50.0], [50.0, 50.0])
            @test_throws DomainError LineParametersProblem(system, earth; frequencies)
        end
        @test_throws DomainError LineParametersProblem(
            system, earth; temperature = NaN, frequencies = [50.0])
    end

    @testset "a problem needs a live conductor" begin
        grounded = build(CableSystem, [wire, wire],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 0,), (core = 0,)], id = "grounded", line_length = 1.0)
        @test_throws ArgumentError LineParametersProblem(grounded, earth; frequencies = [50.0])
    end

    @testset "a cable crosses the surface of the earth, or touches it" begin
        crossing = build(CableSystem, [wire], [Placement(0.0, -0.02)];
            connections = [(core = 1,)], id = "crossing", line_length = 1.0)
        @test_throws ArgumentError LineParametersProblem(crossing, earth; frequencies = [50.0])
        touching = build(CableSystem, [wire], [Placement(0.0, -0.0425)];
            connections = [(core = 1,)], id = "touching", line_length = 1.0)
        @test LineParametersProblem(touching, earth; frequencies = [50.0]).system ===
              touching
    end
end
