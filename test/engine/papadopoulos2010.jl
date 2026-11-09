using QuadGK: quadgk
using SpecialFunctions: besselk

function papadopoulos_axis(integrand, H, κ₁², κ₂², c)
    b = [v for q in (abs(sqrt(κ₁²)), abs(sqrt(κ₂²)), abs(c), 1 / H)
         for v in (q / 4, q, 4q) if !iszero(q)]
    return first(quadgk(integrand, sort(unique([0.0, b..., Inf]))...;
        rtol = 1.0e-13, atol = 0.0, maxevals = 10^8, order = 11))
end

function papadopoulos_entries(input, earth, s)
    μ₀ = 4π * 1e-7
    ε₀ = 8.8541878128e-12
    air, layer = earth.layers
    μ₂ = μ₀ * layer.mu_r
    σ̃₁ = s * ε₀ * air.eps_r
    σ̃₂ = 1 / layer.rho + s * ε₀ * layer.eps_r
    γ₁² = s * μ₀ * σ̃₁
    γ₂² = s * μ₂ * σ̃₂
    kx² = -s^2 * μ₂ * ε₀ * layer.eps_r
    κ₁² = -s^2 * μ₀ * ε₀ * (layer.mu_r * layer.eps_r - air.eps_r)
    κ₂² = γ₂² + kx²
    c = sqrt(κ₂²) * σ̃₁ / σ̃₂
    n = length(input.radius)
    Zg = zeros(ComplexF64, n, n)
    Pg = zeros(ComplexF64, n, n)
    for i in 1:n, j in 1:n

        r_i = input.radius[i]
        h_i = abs(input.y[i])
        h_j = abs(input.y[j])
        H = h_i + h_j
        y = i == j ? r_i : abs(input.x[i] - input.x[j])
        D = i == j ? r_i : hypot(y, h_i - h_j)
        D′ = hypot(y, H)
        Λ = besselk(0, sqrt(κ₂²) * D) - besselk(0, sqrt(κ₂²) * D′)
        Zg[i, j] = s * μ₂ / (2π) * (Λ +
                    2 * papadopoulos_axis(
            λ -> begin
                a₁ = sqrt(λ^2 + κ₁²)
                a₂ = sqrt(λ^2 + κ₂²)
                Δ = a₂ * μ₀ + a₁ * μ₂
                F = μ₀ * exp(-a₂ * H) / Δ
                F * cos(y * λ)
            end, H, κ₁², κ₂², c))
        Pg[i, j] = s / (2π * σ̃₂) * (Λ +
                    2 * papadopoulos_axis(
            λ -> begin
                a₁ = sqrt(λ^2 + κ₁²)
                a₂ = sqrt(λ^2 + κ₂²)
                Δ = a₂ * μ₀ + a₁ * μ₂
                F = μ₀ * exp(-a₂ * H) / Δ
                G = μ₀ * μ₂ * a₂ * (γ₂² - γ₁²) * exp(-a₂ * H) /
                    (Δ * (a₂ * γ₁² * μ₂ + a₁ * γ₂² * μ₀))
                (F + G) * cos(y * λ)
            end, H, κ₁², κ₂², c))
    end
    return (Zg, Pg)
end

@testset "Engine papadopoulos2010 earth" begin
    Engine = LineCableModels.Engine
    formulation = CoaxialFormulation(earth_impedance = :papadopoulos2010,
        earth_admittance = :papadopoulos2010)
    for case in ("two_bare_wires", "two_insulated_wires",
        "coaxial_circular_strands", "coaxial_rectangular_strands")
        @testset "$case" begin
            reference = Dict(REFERENCES)[case]
            input = Engine.CoaxialInput(reference, formulation.options)
            n = length(input.radius)
            target_Z = zeros(ComplexF64, n, n, length(reference.frequencies))
            target_P = zeros(ComplexF64, n, n, length(reference.frequencies))
            for line in eachline(joinpath(
                @__DIR__, "..", "references", case, "papadopoulos2010.csv"))
                fields = split(line, ',')
                startswith(line, "analytical") || continue
                k = searchsortedfirst(reference.frequencies, parse(Float64, fields[2]))
                target = fields[3] == "Zg" ? target_Z : target_P
                target[parse(Int, fields[4]), parse(Int, fields[5]), k] = complex(
                    parse(Float64, fields[6]), parse(Float64, fields[7]))
            end
            result = compute(reference, formulation; trace = true)
            @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target_Z[:, :, k]) /
                          maximum(abs(target_Z[i, i, k]) for i in 1:n)
            for k in eachindex(reference.frequencies)) <= 1.0e-9
            @test maximum(maximum(abs, result.trace.Pg[:, :, k] - target_P[:, :, k]) /
                          maximum(abs(target_P[i, i, k]) for i in 1:n)
            for k in eachindex(reference.frequencies)) <= 1.0e-9
            @test (@inferred compute(reference, formulation)) isa Engine.LineParameters
            for (k, frequency) in enumerate(reference.frequencies)
                target_Z[:, :, k],
                target_P[:, :, k] = papadopoulos_entries(
                    input, reference.earth, im * (2π * frequency))
            end
            @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target_Z[:, :, k]) /
                          maximum(abs(target_Z[i, i, k]) for i in 1:n)
            for k in eachindex(reference.frequencies)) <= 1.0e-8
            @test maximum(maximum(abs, result.trace.Pg[:, :, k] - target_P[:, :, k]) /
                          maximum(abs(target_P[i, i, k]) for i in 1:n)
            for k in eachindex(reference.frequencies)) <= 1.0e-8
        end
    end

    for (earth, frequencies) in (
        (homogeneous(rho = 100.0, eps_r = 10.0),
            Dict(REFERENCES)["two_bare_wires"].frequencies),
        (homogeneous(rho = 100.0, eps_r = 10.0, mu_r = 2.0), [50.0, 1.0e4, 1.0e6]))
        reference = LineParametersProblem(Dict(REFERENCES)["two_bare_wires"].system,
            earth; frequencies)
        input = Engine.CoaxialInput(reference, formulation.options)
        n = length(input.radius)
        target_Z = zeros(ComplexF64, n, n, length(frequencies))
        target_P = zeros(ComplexF64, n, n, length(frequencies))
        result = compute(reference, formulation; trace = true)
        for (k, frequency) in enumerate(frequencies)
            target_Z[:, :, k],
            target_P[:, :, k] = papadopoulos_entries(
                input, earth, im * (2π * frequency))
        end
        @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target_Z[:, :, k]) /
                      maximum(abs(target_Z[i, i, k]) for i in 1:n)
        for k in eachindex(frequencies)) <= 1.0e-8
        @test maximum(maximum(abs, result.trace.Pg[:, :, k] - target_P[:, :, k]) /
                      maximum(abs(target_P[i, i, k]) for i in 1:n)
        for k in eachindex(frequencies)) <= 1.0e-8
    end

    @test_throws "treats no pair of a receiver in the air and a source in the air" compute(
        Dict(REFERENCES)["overhead_bare_wires"], formulation)
    @test_throws "treats no pair of a receiver in the air and a source in the earth" compute(
        Dict(REFERENCES)["mixed_bare_wires"], CoaxialFormulation(earth_impedance = (
            air = :carson1926, earth = :saad1996, mixed = :papadopoulos2010)))
    gamma = CoaxialFormulation(earth_impedance = :papadopoulos2010,
        earth_admittance = :papadopoulos2010, Γ = 1e-3im)
    @test_logs (:warn,
        "the formula :papadopoulos2010 holds for Γ in 0.0 .. 0.0; " *
        "the computation has Γ in 0.001 .. 0.001") LineCableModels.validate(
        LineCableModels.domain(Dict(REFERENCES)["two_bare_wires"], gamma.options), gamma)
    @test_throws ArgumentError CoaxialFormulation(
        earth_impedance = :papadopoulos2010 => (integral = (rtol = -1.0,),))

    reference = Dict(REFERENCES)["two_bare_wires"]
    input = Engine.CoaxialInput(reference, formulation.options)
    buffers = Engine.initialize_buffers(formulation, reference, input)
    media = reference.earth.layers
    s = complex(0.0, 2π * 50.0)
    Engine.to_buffer!(formulation, buffers, input, media, s, 1)
    formula = formulation.earth_impedance.earth
    formula(Val(:Z), Val(:self), 1, buffers, input, media, s, 1)
    formula(Val(:Z), Val(:mutual), 1, 2, buffers, input, media, s, 1)
    @test only(map((f, b, c, m, z) -> @allocated(f(Val(:Z), Val(:self), 1, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
    @test only(map(
        (f, b, c, m, z) -> @allocated(f(Val(:Z), Val(:mutual), 1, 2, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
    formula = formulation.earth_admittance.earth
    formula(Val(:P), Val(:self), 1, buffers, input, media, s, 1)
    formula(Val(:P), Val(:mutual), 1, 2, buffers, input, media, s, 1)
    @test only(map((f, b, c, m, z) -> @allocated(f(Val(:P), Val(:self), 1, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
    @test only(map(
        (f, b, c, m, z) -> @allocated(f(Val(:P), Val(:mutual), 1, 2, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
end
