using QuadGK: quadgk

function wise_integral(H, y, β²)
    return first(quadgk(
        λ -> exp(-H * λ) * cos(y * λ) / (λ + sqrt(λ^2 + β²)),
        sort([0.0, abs(sqrt(β²)) / 4, abs(sqrt(β²)), 4abs(sqrt(β²)),
            1 / (4H), 1 / H, 4 / H, Inf])...;
        rtol = 1.0e-13, atol = 0.0, maxevals = 10^8, order = 11))
end

function wise_series(r, θ)
    s₂ = s₂′ = s₄ = s₄′ = σ₁ = σ₂ = σ₃ = σ₄ = zero(complex(r))
    for k in 0:11
        s₂ += (-1)^k * (r / 2)^(4k + 2) * cos((4k + 2) * θ) /
              (Float64(factorial(big(2k + 1))) * Float64(factorial(big(2k + 2))))
        s₂′ += (-1)^k * (r / 2)^(4k + 2) * sin((4k + 2) * θ) /
               (Float64(factorial(big(2k + 1))) * Float64(factorial(big(2k + 2))))
        s₄ += (-1)^k * (r / 2)^(4k + 4) * cos((4k + 4) * θ) /
              (Float64(factorial(big(2k + 2))) * Float64(factorial(big(2k + 3))))
        s₄′ += (-1)^k * (r / 2)^(4k + 4) * sin((4k + 4) * θ) /
               (Float64(factorial(big(2k + 2))) * Float64(factorial(big(2k + 3))))
        σ₁ += (-1)^k * r^(4k + 1) * cos((4k + 1) * θ) /
              (prod(Float64, 3:2:(4k + 1); init = 1.0)^2 * (4k + 3))
        σ₃ += (-1)^k * r^(4k + 3) * cos((4k + 3) * θ) /
              (prod(Float64, 3:2:(4k + 3))^2 * (4k + 5))
        σ₂ += (-1)^k * (sum(inv, 1:(2k + 2)) - 1 / (4k + 4)) *
              (r / 2)^(4k + 2) * cos((4k + 2) * θ) /
              (Float64(factorial(big(2k + 1))) * Float64(factorial(big(2k + 2))))
        σ₄ += (-1)^k * (sum(inv, 1:(2k + 3)) - 1 / (4k + 6)) *
              (r / 2)^(4k + 4) * cos((4k + 4) * θ) /
              (Float64(factorial(big(2k + 2))) * Float64(factorial(big(2k + 3))))
    end
    γ = exp(MathConstants.γ)
    P = (π / 8) * (1 - s₄) + 0.5 * log(2 / (γ * r)) * s₂ + 0.5 * θ * s₂′ -
        σ₁ / sqrt(2) + 0.5 * σ₂ + σ₃ / sqrt(2)
    Q = 0.25 + 0.5 * log(2 / (γ * r)) * (1 - s₄) - 0.5 * θ * s₄′ +
        σ₁ / sqrt(2) - (π / 8) * s₂ + σ₃ / sqrt(2) - 0.5 * σ₄
    return P + im * Q
end

@testset "Engine wise1934 earth" begin
    Engine = LineCableModels.Engine
    μ₀ = 4π * 1e-7
    ε₀ = 8.8541878128e-12
    formulation = CoaxialFormulation(earth_impedance = :wise1934)
    reference = Dict(REFERENCES)["overhead_bare_wires"]
    wet = LineParametersProblem(reference.system, homogeneous(rho = 100.0, eps_r = 10.0);
        frequencies = reference.frequencies)
    input = Engine.CoaxialInput(reference, formulation.options)
    n = length(input.radius)
    target = zeros(ComplexF64, n, n, length(reference.frequencies))
    for line in eachline(joinpath(
        @__DIR__, "..", "references", "overhead_bare_wires", "wise1934.csv"))
        fields = split(line, ',')
        startswith(line, "analytical") || continue
        k = searchsortedfirst(reference.frequencies, parse(Float64, fields[2]))
        target[parse(Int, fields[4]), parse(Int, fields[5]), k] = complex(
            parse(Float64, fields[6]), parse(Float64, fields[7]))
    end
    result = compute(reference, formulation; trace = true)
    @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target[:, :, k]) /
                  maximum(abs(target[i, i, k]) for i in 1:n)
    for k in eachindex(reference.frequencies)) <= 1.0e-9

    carson = compute(reference, CoaxialFormulation(earth_impedance = :carson1926); trace = true)
    @test maximum(maximum(abs, result.trace.Zg[:, :, k] - carson.trace.Zg[:, :, k]) /
                  maximum(abs(carson.trace.Zg[i, i, k]) for i in 1:n)
    for k in eachindex(reference.frequencies)) <= 1.0e-8

    for problem in (reference, wet)
        result = @inferred compute(problem, formulation; trace = true)
        @test result isa Engine.LineParameters
        for (k, frequency) in enumerate(problem.frequencies)
            s = im * (2π * frequency)
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
                target[i, j, k] = s * μ₀ / (2π) *
                                  ((i == j ? log(2h_i / input.radius[i]) : log(D′ / D)) +
                                   2wise_integral(H, y, β²))
            end
        end
        @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target[:, :, k]) /
                      maximum(abs(target[i, i, k]) for i in 1:n)
        for k in eachindex(problem.frequencies)) <= 1.0e-8
    end

    σ = 0.01
    for r in (0.25, 1.0, 2.0, 5.0)
        α = r^2 / 5
        frequency = α / (2π * μ₀ * σ)
        ω = 2π * frequency
        s = im * ω
        S = sqrt(1 +
                 im * (wet.earth.layers[2].eps_r - 1) * ω * ε₀ * wet.earth.layers[2].rho)
        problem = LineParametersProblem(wet.system, wet.earth; frequencies = [frequency])
        Zg = compute(problem, formulation; trace = true).trace.Zg[:, :, 1]
        mutual = s * μ₀ / (2π) * log(sqrt(5)) +
                 (ω * μ₀ / π) * wise_series(r * S, atan(1 / 2))
        self = s * μ₀ / (2π) * log(2 / 0.0425) +
               (ω * μ₀ / π) * wise_series(2sqrt(α) * S, 0.0)
        expected = [self mutual; mutual self]
        @test maximum(abs, Zg - expected) / maximum(abs, diag(expected)) <= 1.0e-8
    end

    @test_throws "treats no pair of a receiver in the earth and a source in the earth" compute(
        Dict(REFERENCES)["two_bare_wires"], formulation)
    @test_throws "treats no pair of a receiver in the air and a source in the earth" compute(
        Dict(REFERENCES)["mixed_bare_wires"],
        CoaxialFormulation(earth_impedance =
        (air = :wise1934, earth = :saad1996, mixed = :wise1934)))
    @test_throws "gives no potential!" CoaxialFormulation(earth_admittance = :wise1934)
    magnetic = LineParametersProblem(reference.system, homogeneous(rho = 0.1, mu_r = 2.0);
        frequencies = [10.0, 100.0])
    @test_logs (:warn,
        "the formula :wise1934 holds for mu_r_earth in 1.0 .. 1.0; " *
        "the computation has mu_r_earth in 2.0 .. 2.0") compute(magnetic, formulation)
    gamma = CoaxialFormulation(earth_impedance = :wise1934, Γ = 1e-3im)
    @test_logs (:warn,
        "the formula :wise1934 holds for Γ in 0.0 .. 0.0; " *
        "the computation has Γ in 0.001 .. 0.001") LineCableModels.validate(
        LineCableModels.domain(Dict(REFERENCES)["two_bare_wires"], gamma.options), gamma)
    @test_throws ArgumentError CoaxialFormulation(
        earth_impedance = :wise1934 => (integral = (rtol = -1.0,),))
    buffers = Engine.initialize_buffers(formulation, reference, input)
    media = reference.earth.layers
    s = complex(0.0, 2π * 50.0)
    Engine.to_buffer!(formulation, buffers, input, media, s, 1)
    formula = formulation.earth_impedance.air
    formula(Val(:Z), Val(:self), 1, buffers, input, media, s, 1)
    formula(Val(:Z), Val(:mutual), 1, 2, buffers, input, media, s, 1)
    @test only(map((f, b, c, m, z) -> @allocated(f(Val(:Z), Val(:self), 1, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
    @test only(map(
        (f, b, c, m, z) -> @allocated(f(Val(:Z), Val(:mutual), 1, 2, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
end
