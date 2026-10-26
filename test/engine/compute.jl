using LinearAlgebra: I, lu!, norm
using JET

@testset "Engine coaxial computation" begin
    Engine = LineCableModels.Engine
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.3)
    design = build(CableDesign, "coaxial",
        Terminal(:core, Region(:core_metal, Circle(10e-3), copper)),
        Region(:core_insulation, Offset(10e-3), polymer),
        Terminal(:sheath, Region(:sheath_copper, Annulus(20e-3, 21e-3), copper)),
        Region(:sheath_jacket, Offset(3e-3), polymer))
    system = build(CableSystem, [design, design],
        [Placement(0.0, -1.0), Placement(1.0, -1.0)];
        connections = [(core = 1, sheath = 0), (core = 1, sheath = 2)],
        id = "two_cables", line_length = 1.0)
    earth = homogeneous(rho = 100)
    problem = LineParametersProblem(system, earth; frequencies = [50.0, 1e4, 1e6])
    formulation = CoaxialFormulation()
    result = compute(problem, formulation; trace = true)
    unreduced = compute(problem, CoaxialFormulation(; bundle = false, kron = false))
    @test size(result.Z) == (2, 2, 3) && size(result.Y) == (2, 2, 3) &&
          result.terminals == ["cable:1:core+cable:2:core", "cable:2:sheath"]
    # The rows of C follow the conductor order of the unreduced result.
    C = [label in ("cable:1:core", "cable:2:core") ? 1.0 : 0.0
         for label in unreduced.terminals]
    C = hcat(C, Float64.(unreduced.terminals .== "cable:2:sheath"))
    @test all(norm(result.Y[:, :, k] - transpose(C) * unreduced.Y[:, :, k] * C) <=
              1e-12 * norm(result.Y[:, :, k]) for k in eachindex(problem.frequencies))
    @test all(norm(result.Z[:, :, k] -
                   inv(transpose(C) * inv(unreduced.Z[:, :, k]) * C)) <=
              1e-12 * norm(result.Z[:, :, k]) for k in eachindex(problem.frequencies))

    transposed = compute(problem,
        CoaxialFormulation(; bundle = false, kron = false, transposition = true))
    for (computed, original) in ((transposed.Z, unreduced.Z), (transposed.Y, unreduced.Y))
        @test all(norm(computed[:, :, k] -
                       sum(original[circshift(1:4, j), circshift(1:4, j), k]
                  for j in 0:3) / 4) <= 1e-14 * norm(original[:, :, k])
        for k in eachindex(problem.frequencies))
    end

    untraced = compute(problem, formulation; trace = false)
    @test all(size(matrix) == (0, 0, 0) for matrix in untraced.trace) &&
          reinterpret(UInt64, vec(untraced.Z)) == reinterpret(UInt64, vec(result.Z)) &&
          reinterpret(UInt64, vec(untraced.Y)) == reinterpret(UInt64, vec(result.Y))
    @test (@inferred compute(problem, CoaxialFormulation())).Z == result.Z

    input = Engine.CoaxialInput(problem, formulation.options)
    buffers = Engine.initialize_buffers(formulation, problem, input)
    n = length(input.cable)
    c = length(input.radius)
    m = length(input.kept)
    nf = length(problem.frequencies)
    Z = Array{ComplexF64, 3}(undef, m, m, nf)
    Y = Array{ComplexF64, 3}(undef, m, m, nf)
    trace = (Zin = Array{ComplexF64, 3}(undef, n, n, nf),
        Pin = Array{ComplexF64, 3}(undef, n, n, nf),
        Zg = Array{ComplexF64, 3}(undef, c, c, nf),
        Pg = Array{ComplexF64, 3}(undef, c, c, nf))
    options = LineCableModels.ComputationOptions(Engine.LCMBackend; trace = true)
    JET.@test_opt target_modules=(LineCableModels,) Engine.solve!(
        Z, Y, trace, problem, input, formulation, buffers, options)

    # The loop allocates only the pivot vectors of its four LU decompositions.
    for nf in (11, 101)
        problem = LineParametersProblem(system, earth;
            frequencies = 10.0 .^ range(1, 6; length = nf))
        input = Engine.CoaxialInput(problem, formulation.options)
        buffers = Engine.initialize_buffers(formulation, problem, input)
        m = length(input.kept)
        c = length(input.radius)
        Z = Array{ComplexF64, 3}(undef, m, m, nf)
        Y = Array{ComplexF64, 3}(undef, m, m, nf)
        trace = (Zin = Array{ComplexF64, 3}(undef, 0, 0, 0),
            Pin = Array{ComplexF64, 3}(undef, 0, 0, 0),
            Zg = Array{ComplexF64, 3}(undef, 0, 0, 0),
            Pg = Array{ComplexF64, 3}(undef, 0, 0, 0))
        options = LineCableModels.ComputationOptions(Engine.LCMBackend; trace = false)
        Engine.solve!(Z, Y, trace, problem, input, formulation, buffers, options)
        factor = Matrix{ComplexF64}(I, size(buffers.factor)...)
        reduced = Matrix{ComplexF64}(I, size(buffers.reduced)...)
        earth_factor = Matrix{ComplexF64}(I, c, c)
        lu!(factor)
        lu!(reduced)
        lu!(earth_factor)
        pivots = (@allocated lu!(factor)) + (@allocated lu!(factor)) +
                 (@allocated lu!(reduced)) + (@allocated lu!(earth_factor))
        # The arguments pass through an anonymous function, where their types are concrete.
        @test only(map(
            (z, y, t, p, i, f, b,
                o) -> @allocated(LineCableModels.Engine.solve!(
                z, y, t, p, i, f, b, o)),
            (Z,), (Y,), (trace,), (problem,), (input,), (formulation,), (buffers,),
            (options,))) == nf * pivots
    end

    copper32 = Material(:conductor, 1.7241f-8, 1.0f0, 1.0f0, 20.0f0, 0.00393f0)
    polymer32 = Material(:insulator, 1.0f14, 2.3f0)
    design32 = build(CableDesign, "coaxial",
        Terminal(:core, Region(:core_metal, Circle(10.0f-3), copper32)),
        Region(:core_insulation, Offset(10.0f-3), polymer32),
        Terminal(:sheath, Region(:sheath_copper, Annulus(20.0f-3, 21.0f-3), copper32)),
        Region(:sheath_jacket, Offset(3.0f-3), polymer32))
    system32 = build(CableSystem, [design32, design32],
        [Placement(0.0f0, -1.0f0), Placement(1.0f0, -1.0f0)];
        connections = [(core = 1, sheath = 0), (core = 1, sheath = 2)],
        id = "two_cables", line_length = 1.0f0)
    result32 = compute(
        LineParametersProblem(system32, homogeneous(rho = 100.0f0);
            frequencies = Float32[50, 1e4, 1e6]),
        formulation)
    @test eltype(result32.Z) == eltype(result32.Y) == ComplexF64
    @test_throws ArgumentError("LCMBackend does not implement FEM") compute(
        Engine.LCMBackend(), problem, FEM())
end
