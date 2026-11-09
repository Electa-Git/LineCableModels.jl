using SpecialFunctions: besselk

@testset "Engine saad1996 earth" begin
    Engine = LineCableModels.Engine
    μ₀ = 4π * 1e-7
    formulation = CoaxialFormulation(earth_impedance = :saad1996)
    @testset "$case" for case in ("two_bare_wires", "two_insulated_wires",
        "coaxial_circular_strands", "coaxial_rectangular_strands")
        problem = Dict(REFERENCES)[case]
        input = Engine.CoaxialInput(problem, formulation.options)
        n = length(input.radius)
        target = zeros(ComplexF64, n, n, length(problem.frequencies))
        for line in eachline(joinpath(@__DIR__, "..", "references", case, "saad1996.csv"))
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

    for case in ("two_bare_wires", "coaxial_circular_strands")
        reference = Dict(REFERENCES)[case]
        problem = LineParametersProblem(reference.system, reference.earth;
            frequencies = [50.0, 1.0e4, 1.0e6])
        input = Engine.CoaxialInput(problem, formulation.options)
        result = compute(problem, formulation; trace = true)
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
                            besselk(0, m * r) + 2exp(-2m * h) / (4 + m^2 * r^2) :
                            besselk(0, m * D) + 2exp(-m * H) / (4 + m^2 * y^2))
                @test isapprox(result.trace.Zg[i, j, k], expected; rtol = 1.0e-12)
            end
        end
    end

    @test_throws "treats no pair of a receiver in the air and a source in the air" compute(
        Dict(REFERENCES)["overhead_bare_wires"], formulation)
    @test_throws "gives no potential!" CoaxialFormulation(earth_admittance = :saad1996)
    reference = Dict(REFERENCES)["two_bare_wires"]
    magnetic = LineParametersProblem(reference.system, homogeneous(rho = 0.1, mu_r = 2.0);
        frequencies = [10.0, 100.0])
    @test_logs (:warn,
        "the formula :saad1996 holds for mu_r_earth in 1.0 .. 1.0; " *
        "the computation has mu_r_earth in 2.0 .. 2.0") compute(magnetic, formulation)
end
