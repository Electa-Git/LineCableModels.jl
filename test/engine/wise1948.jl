using QuadGK: quadgk
using SpecialFunctions: expint

function wise_potential(H, y, β², a)
    return first(quadgk(
        λ -> exp(-H * λ) * cos(y * λ) / (a * λ + sqrt(λ^2 + β²)),
        sort([0.0,
            [v for q in (abs(sqrt(β²)), abs(sqrt(β²) / a), 1 / H)
             for v in (q / 4, q, 4q)]...,
            Inf])...;
        rtol = 1.0e-13, atol = 0.0, maxevals = 10^8, order = 11))
end

function wise_approximation(g, a, u)
    r₁ = u * (1 - sqrt(1 - 4 / (a + 1))) / 2
    r₂ = u * (1 + sqrt(1 - 4 / (a + 1))) / 2
    return (-r₂ * exp(g * r₁) * (-expint(g * r₁)) +
            r₁ * exp(g * r₂) * (-expint(g * r₂))) / ((a + 1) * (r₂ - r₁))
end

@testset "Engine wise1948 earth" begin
    Engine = LineCableModels.Engine
    μ₀ = 4π * 1e-7
    ε₀ = 8.8541878128e-12
    formulation = CoaxialFormulation(earth_admittance = :wise1948)
    reference = Dict(REFERENCES)["overhead_bare_wires"]
    wet = LineParametersProblem(reference.system, homogeneous(rho = 100.0, eps_r = 10.0);
        frequencies = reference.frequencies)
    input = Engine.CoaxialInput(reference, formulation.options)
    n = length(input.radius)
    target = zeros(ComplexF64, n, n, length(reference.frequencies))
    for line in eachline(joinpath(
        @__DIR__, "..", "references", "overhead_bare_wires", "wise1948.csv"))
        fields = split(line, ',')
        startswith(line, "analytical") && fields[3] == "Pg" || continue
        k = searchsortedfirst(reference.frequencies, parse(Float64, fields[2]))
        target[parse(Int, fields[4]), parse(Int, fields[5]), k] = complex(
            parse(Float64, fields[6]), parse(Float64, fields[7]))
    end
    result = compute(reference, formulation; trace = true)
    @test maximum(maximum(abs, result.trace.Pg[:, :, k] - target[:, :, k]) /
                  maximum(abs(target[i, i, k]) for i in 1:n)
    for k in eachindex(reference.frequencies)) <= 1.0e-9

    for problem in (reference, wet)
        result = @inferred compute(problem, formulation; trace = true)
        @test result isa Engine.LineParameters
        for (k, frequency) in enumerate(problem.frequencies)
            s = im * (2π * frequency)
            ε₁ = ε₀ * problem.earth.layers[1].eps_r
            σ̃₁ = s * ε₁
            σ̃₂ = inv(problem.earth.layers[2].rho) + s * ε₀ * problem.earth.layers[2].eps_r
            a = σ̃₂ / σ̃₁
            β² = s * μ₀ *
                 (inv(problem.earth.layers[2].rho) +
                  s * ε₀ * (problem.earth.layers[2].eps_r - problem.earth.layers[1].eps_r))
            for i in 1:n, j in 1:n

                h_i = abs(input.y[i])
                h_j = abs(input.y[j])
                H = h_i + h_j
                y = abs(input.x[i] - input.x[j])
                D = hypot(y, h_i - h_j)
                D′ = hypot(y, H)
                target[i, j, k] = ((i == j ? log(2h_i / input.radius[i]) : log(D′ / D)) +
                                   2wise_potential(H, y, β², a)) / (2π * ε₁)
            end
        end
        @test maximum(maximum(abs, result.trace.Pg[:, :, k] - target[:, :, k]) /
                      maximum(abs(target[i, i, k]) for i in 1:n)
        for k in eachindex(problem.frequencies)) <= 1.0e-8
    end

    deviations = Float64[]
    for g in (0.5, 2.0, 5.0), x in (1.0, 3.0, 6.0)

        σ = 0.01
        ω = σ / (18x * ε₀)
        s = im * ω
        a = 15 - 18x * im
        β² = s * μ₀ * (σ + s * ε₀ * 14)
        u = sqrt(β²) / abs(sqrt(β²))
        h = g / (2abs(sqrt(β²)))
        system = build(CableSystem, [reference.system.designs[1]], [Placement(0.0, h)];
            connections = [(core = 1,)], id = "wise", line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 15.0);
            frequencies = [ω / (2π)])
        Pg = compute(problem, formulation; trace = true).trace.Pg
        MN = (2π * ε₀ * Pg[1, 1, 1] - log(2h / 0.0425)) / 2
        A = wise_approximation(g, a, u)
        push!(deviations,
            max(abs(real(MN) - real(A)) / abs(real(A)),
                abs(imag(MN) - imag(A)) / abs(imag(A))))
    end
    @test maximum(deviations) < 0.01

    @test_throws "treats no pair of a receiver in the earth and a source in the earth" compute(
        Dict(REFERENCES)["two_bare_wires"], formulation)
    @test_throws "treats no pair of a receiver in the air and a source in the earth" compute(
        Dict(REFERENCES)["mixed_bare_wires"],
        CoaxialFormulation(earth_admittance =
        (air = :wise1948, earth = :ideal, mixed = :wise1948)))
    @test_throws "gives no impedance!" CoaxialFormulation(earth_impedance = :wise1948)
    magnetic = LineParametersProblem(reference.system, homogeneous(rho = 0.1, mu_r = 2.0);
        frequencies = [10.0, 100.0])
    @test_logs (:warn,
        "the formula :wise1948 holds for mu_r_earth in 1.0 .. 1.0; " *
        "the computation has mu_r_earth in 2.0 .. 2.0") compute(magnetic, formulation)
    gamma = CoaxialFormulation(earth_admittance = :wise1948, Γ = 1e-3im)
    @test_logs (:warn,
        "the formula :wise1948 holds for Γ in 0.0 .. 0.0; " *
        "the computation has Γ in 0.001 .. 0.001") LineCableModels.validate(
        LineCableModels.domain(Dict(REFERENCES)["two_bare_wires"], gamma.options), gamma)
    @test_throws ArgumentError CoaxialFormulation(
        earth_admittance = :wise1948 => (integral = (rtol = -1.0,),))
    buffers = Engine.initialize_buffers(formulation, reference, input)
    media = reference.earth.layers
    s = complex(0.0, 2π * 50.0)
    Engine.to_buffer!(formulation, buffers, input, media, s, 1)
    formula = formulation.earth_admittance.air
    formula(Val(:P), Val(:self), 1, buffers, input, media, s, 1)
    formula(Val(:P), Val(:mutual), 1, 2, buffers, input, media, s, 1)
    @test only(map((f, b, c, m, z) -> @allocated(f(Val(:P), Val(:self), 1, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
    @test only(map(
        (f, b, c, m, z) -> @allocated(f(Val(:P), Val(:mutual), 1, 2, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
end
