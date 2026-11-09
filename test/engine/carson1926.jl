using QuadGK: quadgk

function carson_series(r, θ)
    s₂ = s₂′ = s₄ = s₄′ = σ₁ = σ₂ = σ₃ = σ₄ = 0.0
    for k in 0:7
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
    return complex(P, Q)
end

@testset "Engine carson1926 earth" begin
    Engine = LineCableModels.Engine
    μ₀ = 4π * 1e-7
    formulation = CoaxialFormulation(earth_impedance = :carson1926)
    reference = Dict(REFERENCES)["overhead_bare_wires"]
    input = Engine.CoaxialInput(reference, formulation.options)
    n = length(input.radius)
    target = zeros(ComplexF64, n, n, length(reference.frequencies))
    for line in eachline(joinpath(
        @__DIR__, "..", "references", "overhead_bare_wires", "carson1926.csv"))
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

    h = 1.0
    y = 1.0
    r = 0.0425
    σ = 10.0
    for (k, frequency) in enumerate(reference.frequencies)
        s = im * (2π * frequency)
        γ² = s * μ₀ * σ
        H = 2h
        for i in 1:n, j in 1:n

            J = first(quadgk(
                λ -> exp(-H * λ) * cos((i == j ? 0.0 : y) * λ) /
                     (λ + sqrt(λ^2 + γ²)),
                sort([0.0, abs(sqrt(γ²)) / 4, abs(sqrt(γ²)), 4abs(sqrt(γ²)),
                    1 / (4H), 1 / H, 4 / H, Inf])...;
                rtol = 1.0e-13, atol = 0.0, maxevals = 10^8, order = 11))
            target[i, j, k] = s * μ₀ / (2π) *
                              ((i == j ? log(2h / r) : log(hypot(y, H) / y)) + 2J)
        end
    end
    @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target[:, :, k]) /
                  maximum(abs(target[i, i, k]) for i in 1:n)
    for k in eachindex(reference.frequencies)) <= 1.0e-8

    for r in (0.25, 1.0, 2.0, 5.0)
        α = r^2 / 5
        frequency = α / (2π * μ₀ * σ)
        ω = 2π * frequency
        s = im * ω
        problem = LineParametersProblem(reference.system, reference.earth;
            frequencies = [frequency])
        Zg = compute(problem, formulation; trace = true).trace.Zg[:, :, 1]
        mutual = s * μ₀ / (2π) * log(sqrt(5)) +
                 (ω * μ₀ / π) * carson_series(r, atan(1 / 2))
        self = s * μ₀ / (2π) * log(2 / 0.0425) +
               (ω * μ₀ / π) * carson_series(2sqrt(α), 0.0)
        expected = [self mutual; mutual self]
        @test maximum(abs, Zg - expected) / maximum(abs, diag(expected)) <= 1.0e-8
    end

    @test_throws "treats no pair of a receiver in the earth and a source in the earth" compute(
        Dict(REFERENCES)["two_bare_wires"], formulation)
    @test_throws "treats no pair of a receiver in the air and a source in the earth" compute(
        Dict(REFERENCES)["mixed_bare_wires"],
        CoaxialFormulation(earth_impedance =
        (air = :carson1926, earth = :saad1996, mixed = :carson1926)))
    @test_throws "gives no potential!" CoaxialFormulation(earth_admittance = :carson1926)
    magnetic = LineParametersProblem(reference.system, homogeneous(rho = 0.1, mu_r = 2.0);
        frequencies = [10.0, 100.0])
    @test_logs (:warn,
        "the formula :carson1926 holds for mu_r_earth in 1.0 .. 1.0; " *
        "the computation has mu_r_earth in 2.0 .. 2.0") compute(magnetic, formulation)
    gamma = CoaxialFormulation(earth_impedance = :carson1926, earth_admittance = :ideal,
        Γ = 1e-3im)
    @test_logs (:warn,
        "the formula :carson1926 holds for Γ in 0.0 .. 0.0; " *
        "the computation has Γ in 0.001 .. 0.001") (:warn,
        "the formula :ideal holds for Γ in 0.0 .. 0.0; " *
        "the computation has Γ in 0.001 .. 0.001") compute(reference, gamma)
    @test_throws ArgumentError CoaxialFormulation(
        earth_impedance = :carson1926 => (integral = (rtol = -1.0,),))
    @test (@inferred compute(
        reference, CoaxialFormulation(earth_impedance = :carson1926,
            earth_admittance = :ideal))) isa Engine.LineParameters

    buffers = Engine.initialize_buffers(formulation, reference, input)
    media = reference.earth.layers
    s = complex(0.0, 2π * 50.0)
    formula = formulation.earth_impedance.air
    formula(Val(:Z), Val(:self), 1, buffers, input, media, s, 1)
    formula(Val(:Z), Val(:mutual), 1, 2, buffers, input, media, s, 1)
    @test only(map((f, b, c, m, z) -> @allocated(f(Val(:Z), Val(:self), 1, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
    @test only(map(
        (f, b, c, m, z) -> @allocated(f(Val(:Z), Val(:mutual), 1, 2, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
end
