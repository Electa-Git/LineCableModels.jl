# A formulation and a backend that the test declares. `compute` runs a formulation on the
# backend that its type names, and the method of the backend computes.
struct Routed <: LineCableModels.AbstractFormulation end
struct Recorder <: LineCableModels.AbstractBackend end
LineCableModels.backend(::Type{Routed}) = Recorder()
function LineCableModels.compute(
        ::Recorder, problem::LineParametersProblem, formulation::Routed; kw...)
    return (problem, formulation, values(kw))
end

@testset "Commons" begin
    @testset "validate returns what it checks" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        @test LineCableModels.validate(copper) === copper
        earth = homogeneous(rho = 100, eps_r = 10)
        @test LineCableModels.validate(earth) === earth
    end

    @testset "resolve gives a layer its place around a part" begin
        @test LineCableModels.resolve(Circle(1.0e-3), Offset(1.0e-3)) ==
              Annulus(1.0e-3, 2.0e-3)
    end

    @testset "build gives a cable its regions and terminals" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        wire = build(CableDesign, "wire",
            Terminal(:core, Region(:core_metal, Circle(1.0e-3), copper)))
        @test [(row.name, row.regions) for row in wire.terminals] == [(:core, [1])]
    end

    @testset "compute runs a formulation on the backend that its type names" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        wire = build(CableDesign, "bare_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)))
        system = build(
            CableSystem, [wire, wire], [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_wires", line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 100); frequencies = [50.0])
        @test compute(problem, Routed(); trace = true) ==
              (problem, Routed(), (trace = true,))
    end

    @testset "options hold every declared key in the type of its default" begin
        FormulationOptions = LineCableModels.FormulationOptions
        ComputationOptions = LineCableModels.ComputationOptions
        defaults = LineCableModels.defaults
        reductions = FormulationOptions(LineParametersFormulation; kron = false)
        @test (reductions.bundle, reductions.kron, reductions.transposition) ==
              (true, false, false)
        @test typeof(reductions) === typeof(FormulationOptions(LineParametersFormulation))
        @test getfield(FormulationOptions(FEM), :values) == defaults(FEM)
        @test FormulationOptions(FEM; Γ = 1).Γ === complex(1.0)
        @test @inferred((options -> options.Γ)(FormulationOptions(FEM))) === complex(0.0)
        @test FormulationOptions(FEM; Γ = 0.5im) == FormulationOptions(FEM; Γ = 0.5im)
        @test FormulationOptions(FEM; Γ = 0.5im) != FormulationOptions(FEM)
        @test_throws ArgumentError FormulationOptions(LineParametersFormulation; bundled = false)
        @test_throws DomainError FormulationOptions(FEM; Γ = complex(NaN))
        @test ComputationOptions(LineCableModels.CoaxialBackend; trace = true).trace
        @test_throws ArgumentError ComputationOptions(
            LineCableModels.CoaxialBackend; traced = true)
    end

    @testset "select fits a selection to its declared default" begin
        select = LineCableModels.select
        @test select((rtol = 1.0e-8, maxevals = 10), (maxevals = 20.0,)) ===
              (rtol = 1.0e-8, maxevals = 20)
        @test eltype(select(complex(0.0), [1, 2])) === ComplexF64
        @test_throws ArgumentError select((rtol = 1.0e-8,), (atol = 0.0,))
    end
end
