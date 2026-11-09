@testset "Engine gary1976 earth" begin
    Engine = LineCableModels.Engine
    μ₀ = 4π * 1e-7
    formulation = CoaxialFormulation(earth_impedance = :gary1976)
    @testset "$case" for case in ("overhead_bare_wires",)
        problem = Dict(REFERENCES)[case]
        input = Engine.CoaxialInput(problem, formulation.options)
        n = length(input.radius)
        target = zeros(ComplexF64, n, n, length(problem.frequencies))
        for line in eachline(joinpath(@__DIR__, "..", "references", case, "gary1976.csv"))
            fields = split(line, ',')
            startswith(line, "analytical") || continue
            k = searchsortedfirst(problem.frequencies, parse(Float64, fields[2]))
            target[parse(Int, fields[4]), parse(Int, fields[5]), k] = complex(
                parse(Float64, fields[6]), parse(Float64, fields[7]))
        end
        result = compute(problem, formulation; trace = true)
        @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target[:, :, k]) /
                      maximum(abs(target[i, i, k]) for i in 1:n)
        for k in eachindex(problem.frequencies)) <= 1.0e-12
        @test (@inferred compute(problem, formulation)) isa Engine.LineParameters
    end

    for case in ("overhead_bare_wires",)
        reference = Dict(REFERENCES)[case]
        problem = LineParametersProblem(reference.system, reference.earth;
            frequencies = [50.0, 1.0e4, 1.0e6])
        input = Engine.CoaxialInput(problem, formulation.options)
        result = compute(problem, formulation; trace = true)
        for (k, frequency) in enumerate(problem.frequencies)
            s = im * (2π * frequency)
            σ = inv(last(problem.earth.layers).rho)
            p = 1 / sqrt(s * μ₀ * σ)
            for i in eachindex(input.radius), j in eachindex(input.radius)

                r = input.radius[i]
                h = abs(input.y[i])
                D = hypot(input.x[i] - input.x[j], input.y[i] - input.y[j])
                S = sqrt((input.x[i] - input.x[j])^2 +
                         (abs(input.y[i]) + abs(input.y[j]) + 2p)^2)
                expected = s * μ₀ / (2π) * (i == j ? log(2(h + p) / r) : log(S / D))
                @test isapprox(result.trace.Zg[i, j, k], expected; rtol = 1.0e-12)
            end
        end
    end

    @test_throws "treats no pair of a receiver in the earth and a source in the earth" compute(
        Dict(REFERENCES)["two_bare_wires"], formulation)
    @test_throws "gives no potential!" CoaxialFormulation(earth_admittance = :gary1976)
    @test_throws DomainError compute(Dict(REFERENCES)["mixed_bare_wires"], formulation)
    reference = Dict(REFERENCES)["overhead_bare_wires"]
    magnetic = LineParametersProblem(reference.system, homogeneous(rho = 0.1, mu_r = 2.0);
        frequencies = [10.0, 100.0])
    @test_logs (:warn,
        "the formula :gary1976 holds for mu_r_earth in 1.0 .. 1.0; " *
        "the computation has mu_r_earth in 2.0 .. 2.0") compute(magnetic, formulation)
end
