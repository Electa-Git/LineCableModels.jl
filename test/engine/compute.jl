using LinearAlgebra: I, lu!, norm
using JET

# An earth-return identifier that the test declares for the mixed pairs, so that one slot holds
# three distinct formulas.
function LineCableModels.defaults(::Type{LineCableModels.Engine.EarthReturn.Formula{:crossing}})
    (;)
end
function LineCableModels.validate(
        formula::LineCableModels.Engine.EarthReturn.Formula{:crossing},
        ::typeof(LineCableModels.Engine.potential!))
    return formula
end
function (::LineCableModels.Engine.EarthReturn.Formula{:crossing})(
        ::Val{:P}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    return complex(1.0e10)
end

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

    @test_throws ArgumentError("LCMBackend does not implement FEM") compute(
        Engine.LCMBackend(), problem, FEM())

    @testset "a layered earth under the bottommost model" begin
        reference = Dict(REFERENCES)["two_bare_wires"]
        air, bottom = reference.earth.layers
        layered = LineParametersProblem(reference.system,
            EarthModel((air, EarthLayer(100.0, 1.0, 1.0, 5.0), bottom));
            temperature = reference.temperature, frequencies = reference.frequencies)
        result = compute(layered, CoaxialFormulation())
        expected = compute(reference, CoaxialFormulation())
        @test reinterpret(UInt64, result.Z) == reinterpret(UInt64, expected.Z)
        @test reinterpret(UInt64, result.Y) == reinterpret(UInt64, expected.Y)
    end

    @testset "a layered earth under the equivalent earth of Xue, against its FEM row" begin
        reference = Dict(REFERENCES)["two_bare_wires"]
        earth = EarthModel((EarthLayer(Inf, 1.0, 1.0), EarthLayer(100.0, 10.0, 1.0, 5.0),
            EarthLayer(500.0, 20.0, 1.0, 10.0), EarthLayer(50.0, 5.0, 1.0)))
        problem = LineParametersProblem(reference.system, earth;
            temperature = reference.temperature, frequencies = [1000.0])
        result = compute(problem, CoaxialFormulation(equivalent_earth = :xue2021))
        rows = Dict((source, quantity) => zeros(ComplexF64, 2, 2)
        for source in ("analytical", "fem"), quantity in ("Z", "Y"))
        path = joinpath(@__DIR__, "..", "references", "two_bare_wires", "xue2021.csv")
        for line in eachline(path)
            fields = split(line, ',')
            haskey(rows, (fields[1], get(fields, 3, ""))) || continue
            row = parse(Int, fields[4])
            column = parse(Int, fields[5])
            value = complex(parse(Float64, fields[6]), parse(Float64, fields[7]))
            rows[(fields[1], fields[3])][row, column] = value
        end
        # The FEM deviations are the README's status of the file, at its printed precision.
        for (quantity, computed, status) in
            (("Z", result.Z, 0.00503), ("Y", result.Y, 0.00735))
            analytical = rows[("analytical", quantity)]
            scale = maximum(abs(analytical[i, i]) for i in 1:2)
            @test maximum(abs, computed[:, :, 1] - analytical) <= 1e-9 * scale
            @test round(maximum(abs, rows[("fem", quantity)] - analytical) / scale;
                sigdigits = 3) == status
        end
    end

    @testset "the layer map and the media terms" begin
        mixed = Dict(REFERENCES)["mixed_bare_wires"]
        @test Engine.CoaxialInput(mixed, CoaxialFormulation().options).earth_layer == [1, 2]
        reference = Dict(REFERENCES)["two_bare_wires"]
        air, bottom = reference.earth.layers
        design = first(reference.system.designs)
        system = build(CableSystem, [design, design],
            [Placement(0.0, -1.0), Placement(1.0, -6.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_bare_wires", line_length = 1.0)
        problem = LineParametersProblem(system,
            EarthModel((air, EarthLayer(100.0, 1.0, 1.0, 5.0), bottom));
            frequencies = [50.0, 1.0e4])
        formulation = CoaxialFormulation()
        @test Engine.CoaxialInput(problem, formulation.options).earth_layer == [2, 3]
        result = compute(problem, formulation)
        expected = compute(
            LineParametersProblem(system, EarthModel((air, bottom));
                frequencies = [50.0, 1.0e4]),
            formulation)
        @test isequal(result.Z, expected.Z)
        @test isequal(result.Y, expected.Y)
        problem = LineParametersProblem(
            reference.system, homogeneous(rho = 0.1, mu_r = 2.0);
            frequencies = [50.0])
        input = Engine.CoaxialInput(problem, formulation.options)
        buffers = Engine.initialize_buffers(formulation, problem, input)
        media = problem.earth.layers
        s = im * (2π * 50.0)
        Engine.to_buffer!(CoaxialFormulation(), buffers, input, media, s, 1)
        ε₀ = 8.8541878128e-12
        μ₀ = 4π * 1e-7
        @test isapprox(buffers.σ̃, [s * ε₀, 1 / 0.1 + s * ε₀]; rtol = 1.0e-15) &&
              isapprox(buffers.μ, [μ₀, 2μ₀]; rtol = 1.0e-15)
    end

    @testset "the three arms of a slot" begin
        mixed = Dict(REFERENCES)["mixed_bare_wires"]
        arms = CoaxialFormulation(;
            earth_admittance = (air = :ideal, earth = :martins2026, mixed = :crossing))
        @test length(unique(map(typeof, Tuple(arms.earth_admittance)))) == 3
        input = Engine.CoaxialInput(mixed, arms.options)
        buffers = Engine.initialize_buffers(arms, mixed, input)
        media = mixed.earth.layers
        s = complex(0.0, 2π * first(mixed.frequencies))
        Engine.to_buffer!(arms, buffers, input, media, s, 1)
        map(formula -> Engine.to_buffer!(formula, buffers, input, media, s, 1), arms.formulas)
        Pg = zeros(ComplexF64, 2, 2)
        Engine.potential!(arms, Pg, buffers, input, media, s, 1)
        @test Pg[1, 1] ==
              arms.earth_admittance.air(Val(:P), Val(:self), 1, buffers, input, media, s, 1)
        @test Pg[2, 2] == Engine.from_buffer(arms.earth_admittance.earth, buffers).Pg[2, 2]
        @test Pg[1, 2] == Pg[2, 1] == complex(1.0e10)
        JET.@test_opt target_modules=(LineCableModels,) Engine.potential!(
            arms, Pg, buffers, input, media, s, 1)

        @test_throws "treats no pair of a receiver in the air and a source in the earth" compute(
            mixed, arms)
        closed = CoaxialFormulation(earth_impedance = (air = :gary1976, earth = :saad1996))
        places = CoaxialFormulation(earth_impedance = (
            air = :gary1976, earth = :saad1996, mixed = :lucca1994))
        @test compute(mixed, closed) isa Engine.LineParameters
        bytes = Int[]
        for formulation in (CoaxialFormulation(), closed, places)
            input = Engine.CoaxialInput(mixed, formulation.options)
            buffers = Engine.initialize_buffers(formulation, mixed, input)
            m = length(input.kept)
            nf = length(mixed.frequencies)
            Z = Array{ComplexF64, 3}(undef, m, m, nf)
            Y = similar(Z)
            trace = (Zin = Array{ComplexF64, 3}(undef, 0, 0, 0),
                Pin = Array{ComplexF64, 3}(undef, 0, 0, 0),
                Zg = Array{ComplexF64, 3}(undef, 0, 0, 0),
                Pg = Array{ComplexF64, 3}(undef, 0, 0, 0))
            options = LineCableModels.ComputationOptions(Engine.LCMBackend; trace = false)
            Engine.solve!(Z, Y, trace, mixed, input, formulation, buffers, options)
            push!(bytes,
                only(map(
                    (z,
                        y,
                        t,
                        p,
                        i,
                        f,
                        b,
                        o) -> @allocated(LineCableModels.Engine.solve!(
                        z, y, t, p, i, f, b, o)),
                    (Z,), (Y,), (trace,), (mixed,), (input,), (formulation,), (buffers,),
                    (options,))))
        end
        @test bytes[1] == bytes[2]
        @test bytes[1] == bytes[3]
    end
end
