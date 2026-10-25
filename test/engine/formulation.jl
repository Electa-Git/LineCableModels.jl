@testset "Engine formulations" begin
    Engine = LineCableModels.Engine
    backend = LineCableModels.backend

    @testset "a backend runs each formulation type" begin
        @test backend(LineParametersFormulation) === LineCableModels.CoaxialBackend()
        @test backend(typeof(FEM())) === LineCableModels.GmshBackend()
    end

    @testset "the finite-element formulation and its options" begin
        fem = FEM()
        @test (fem.insulation_admittance, fem.semicon_admittance, fem.resistivity,
            fem.properties) === (Engine.InsulationAdmittance.Formula(Val(:lossless)),
            Engine.SemiconAdmittance.Formula(Val(:lossless)),
            Engine.Resistivity.Formula(Val(:linear)),
            Engine.EarthProperties.Formula(Val(:constant)))
        @test FEM(; Γ = [complex(0.0, 1.0)]).options.Γ == [complex(0.0, 1.0)]
        @test_throws DomainError FEM(; Γ = complex(NaN))
        @test_throws DomainError FEM(; Γ = ComplexF64[])
    end

    @testset "the finite-element computation needs the Gmsh extension" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        wire = build(CableDesign, "bare_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)))
        system = build(
            CableSystem, [wire, wire], [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_wires", line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 100); frequencies = [50.0])
        @test_throws "using Gmsh" compute(problem, FEM())
    end
end
