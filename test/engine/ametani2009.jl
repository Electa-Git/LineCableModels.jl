@testset "Engine ametani2009 earth" begin
    Engine = LineCableModels.Engine
    μ₀ = 4π * 1e-7
    formulation = CoaxialFormulation(earth_impedance = (
        air = :gary1976, earth = :saad1996, mixed = :ametani2009))
    reference = Dict(REFERENCES)["mixed_bare_wires"]
    input = Engine.CoaxialInput(reference, formulation.options)
    n = length(input.radius)
    target = zeros(ComplexF64, n, n, length(reference.frequencies))
    for line in eachline(joinpath(
        @__DIR__, "..", "references", "mixed_bare_wires", "ametani2009.csv"))
        fields = split(line, ',')
        startswith(line, "analytical") || continue
        k = searchsortedfirst(reference.frequencies, parse(Float64, fields[2]))
        target[parse(Int, fields[4]), parse(Int, fields[5]), k] = complex(
            parse(Float64, fields[6]), parse(Float64, fields[7]))
    end
    result = compute(reference, formulation; trace = true)
    @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target[:, :, k]) /
                  maximum(abs(target[i, i, k]) for i in 1:n)
    for k in eachindex(reference.frequencies)) <= 1.0e-12
    @test (@inferred compute(reference, formulation)) isa Engine.LineParameters

    design = reference.system.designs[1]
    for placements in ([Placement(0.0, 5.0), Placement(1.0, -0.5)],
        [Placement(1.0, -0.5), Placement(0.0, 5.0)])
        system = build(CableSystem, [design, design], placements;
            connections = [(core = 1,), (core = 2,)], id = "uneven", line_length = 1.0)
        problem = LineParametersProblem(system, reference.earth;
            frequencies = [50.0, 1.0e4, 1.0e6])
        result = compute(problem, formulation; trace = true)
        for (k, frequency) in enumerate(problem.frequencies)
            s = im * (2π * frequency)
            h₁ = 5.0
            h₂ = 0.5
            y = 1.0
            σ = 1 / 0.1
            m = sqrt(s * μ₀ * σ)
            h_e = inv(m)
            H = h₁ + h₂ + 2h_e
            S = sqrt(H^2 + y^2)
            D = sqrt((h₁ + h₂)^2 + y^2)
            Z_m = s * μ₀ / (2π) * exp(-h₂ / h_e) * log(S / D)
            @test isapprox(result.trace.Zg[1, 2, k], Z_m; rtol = 1.0e-12)
            @test isapprox(result.trace.Zg[2, 1, k], Z_m; rtol = 1.0e-12)
        end
    end

    @test_throws "treats no pair of a receiver in the air and a source in the air" compute(
        reference, CoaxialFormulation(earth_impedance = :ametani2009))
    @test_throws "treats no pair of a receiver in the earth and a source in the earth" compute(
        Dict(REFERENCES)["two_bare_wires"], CoaxialFormulation(earth_impedance = :ametani2009))
    @test_throws "gives no potential!" CoaxialFormulation(earth_admittance = :ametani2009)
    magnetic = LineParametersProblem(reference.system, homogeneous(rho = 0.1, mu_r = 2.0);
        frequencies = [10.0, 100.0])
    @test_throws r":ametani2009.*earth layer 1 has mu_r = 2.0" compute(
        magnetic, CoaxialFormulation(earth_impedance = (mixed = :ametani2009,)))
end
