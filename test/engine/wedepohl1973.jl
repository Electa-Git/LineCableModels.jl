@testset "Engine wedepohl1973 earth" begin
    Engine = LineCableModels.Engine
    μ₀ = 4π * 1e-7
    formulation = CoaxialFormulation(earth_impedance = :wedepohl1973)
    @testset "$case" for case in ("two_bare_wires", "two_insulated_wires",
        "coaxial_circular_strands", "coaxial_rectangular_strands")
        problem = Dict(REFERENCES)[case]
        input = Engine.CoaxialInput(problem, formulation.options)
        n = length(input.radius)
        target = zeros(ComplexF64, n, n, length(problem.frequencies))
        for line in
            eachline(joinpath(@__DIR__, "..", "references", case, "wedepohl1973.csv"))
            fields = split(line, ',')
            startswith(line, "analytical") || continue
            k = searchsortedfirst(problem.frequencies, parse(Float64, fields[2]))
            target[parse(Int, fields[4]), parse(Int, fields[5]), k] = complex(
                parse(Float64, fields[6]), parse(Float64, fields[7]))
        end
        result = if startswith(case, "two_")
            @test_logs (:warn, r":wedepohl1973") compute(problem, formulation; trace = true)
        else
            @test_logs compute(problem, formulation; trace = true)
        end
        @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target[:, :, k]) /
                      maximum(abs(target[i, i, k]) for i in 1:n)
        for k in eachindex(problem.frequencies)) <= 1.0e-12
        inferred = if startswith(case, "two_")
            @test_logs (:warn, r":wedepohl1973") @inferred compute(problem, formulation)
        else
            @test_logs @inferred compute(problem, formulation)
        end
        @test inferred isa Engine.LineParameters
    end

    for case in ("two_bare_wires", "coaxial_circular_strands")
        reference = Dict(REFERENCES)[case]
        problem = LineParametersProblem(reference.system, reference.earth;
            frequencies = [50.0, 1.0e4, 1.0e6])
        input = Engine.CoaxialInput(problem, formulation.options)
        result = if case == "two_bare_wires"
            @test_logs (:warn, r":wedepohl1973") compute(problem, formulation; trace = true)
        else
            @test_logs compute(problem, formulation; trace = true)
        end
        for (k, frequency) in enumerate(problem.frequencies)
            s = im * (2π * frequency)
            m = sqrt(s * μ₀ / last(problem.earth.layers).rho)
            for i in eachindex(input.radius), j in eachindex(input.radius)

                r = input.radius[i]
                h = abs(input.y[i])
                H = h + abs(input.y[j])
                y = abs(input.x[i] - input.x[j])
                D = hypot(y, input.y[i] - input.y[j])
                expected = s * μ₀ / (2π) * (i == j ?
                            -log(exp(MathConstants.γ) * m * r / 2) + 1 / 2 - 4m * h / 3 :
                            -log(exp(MathConstants.γ) * m * D / 2) + 1 / 2 - 2m * H / 3)
                @test isapprox(result.trace.Zg[i, j, k], expected; rtol = 1.0e-12)
            end
        end
    end

    @test_throws "treats no pair of a receiver in the air and a source in the air" compute(
        Dict(REFERENCES)["overhead_bare_wires"], formulation)
    @test_throws "gives no potential!" CoaxialFormulation(earth_admittance = :wedepohl1973)
    reference = Dict(REFERENCES)["two_bare_wires"]
    magnetic = LineParametersProblem(reference.system, homogeneous(rho = 0.1, mu_r = 2.0);
        frequencies = [10.0, 100.0])
    result = @test_logs compute(magnetic, formulation; trace = true)
    input = Engine.CoaxialInput(magnetic, formulation.options)
    μ = 2μ₀
    for (k, frequency) in enumerate(magnetic.frequencies)
        s = im * (2π * frequency)
        m = sqrt(s * μ / 0.1)
        for i in eachindex(input.radius), j in eachindex(input.radius)

            r = input.radius[i]
            h = abs(input.y[i])
            H = h + abs(input.y[j])
            D = hypot(input.x[i] - input.x[j], input.y[i] - input.y[j])
            expected = s * μ / (2π) * (i == j ?
                        -log(exp(MathConstants.γ) * m * r / 2) + 1 / 2 - 4m * h / 3 :
                        -log(exp(MathConstants.γ) * m * D / 2) + 1 / 2 - 2m * H / 3)
            @test isapprox(result.trace.Zg[i, j, k], expected; rtol = 1.0e-12)
        end
    end
    @test_logs compute(
        LineParametersProblem(reference.system, reference.earth;
            frequencies = [10.0, 100.0]),
        formulation)
    @test_logs (:warn, r":wedepohl1973") compute(
        LineParametersProblem(reference.system, reference.earth; frequencies = [1.0e4]),
        formulation)
    system = build(CableSystem, [first(reference.system.designs)],
        [first(reference.system.placements)]; connections = [(core = 1,)], id = "one_wire",
        line_length = 1.0)
    @test_logs compute(
        LineParametersProblem(system, reference.earth;
            frequencies = [1.0e4]), formulation)
    @test_logs (:warn, r":wedepohl1973") compute(
        LineParametersProblem(system, reference.earth; frequencies = [1.0e6]), formulation)
end
