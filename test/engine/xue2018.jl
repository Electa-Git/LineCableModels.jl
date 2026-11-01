using QuadGK: quadgk
using SpecialFunctions: besselk

function xue_axis(integrand, H, κ₁², κ₂², c)
    b = [v for q in (abs(sqrt(κ₁²)), abs(sqrt(κ₂²)), abs(c), 1 / H)
         for v in (q / 4, q, 4q) if !iszero(q)]
    return first(quadgk(integrand, sort(unique([0.0, b..., Inf]))...;
        rtol = 1.0e-13, atol = 0.0, maxevals = 10^8, order = 11))
end

function xue_entries(input, earth, s, Γ)
    μ₀ = 4π * 1e-7
    ε₀ = 8.8541878128e-12
    air, layer = earth.layers
    σ̃₁ = s * ε₀ * air.eps_r
    σ̃₂ = 1 / layer.rho + s * ε₀ * layer.eps_r
    k_a² = -s * μ₀ * σ̃₁
    k_e² = -s * μ₀ * σ̃₂
    k_v² = -Γ^2
    τ = k_a² / k_e²
    κ₁² = k_v² - k_a²
    κ₂² = k_v² - k_e²
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
        Zg[i, j] = s * μ₀ / (2π) * (Λ +
                    2 * xue_axis(
            λ -> begin
                u₀ = sqrt(λ^2 + k_v² - k_a²)
                u₁ = sqrt(λ^2 + k_v² - k_e²)
                e = exp(-u₁ * H) * cos(y * λ)
                Δ₄ = e * λ^2 / ((u₀ + u₁) * u₁^2)
                Δ₆ = e / ((u₀ + u₁) * u₁^2)
                Δ₇ = e / ((u₀ + τ * u₁) * u₁^2)
                Δ₄ - k_e² * Δ₆ + k_v² * Δ₇
            end, H, κ₁², κ₂², c))
        Pg[i, j] = s / (2π * σ̃₂) * (Λ +
                    2 * xue_axis(
            λ -> begin
                u₀ = sqrt(λ^2 + k_v² - k_a²)
                u₁ = sqrt(λ^2 + k_v² - k_e²)
                e = exp(-u₁ * H) * cos(y * λ)
                Δ₅ = e * λ^2 / ((u₀ + τ * u₁) * u₁^2)
                Δ₆ = e / ((u₀ + u₁) * u₁^2)
                Δ₇ = e / ((u₀ + τ * u₁) * u₁^2)
                Δ₅ - k_e² * Δ₆ + k_v² * Δ₇
            end, H, κ₁², κ₂², c))
    end
    return (Zg, Pg)
end

@testset "Engine xue2018 earth" begin
    Engine = LineCableModels.Engine
    formulation = CoaxialFormulation(earth_impedance = :xue2018,
        earth_admittance = :xue2018)
    for case in ("two_bare_wires", "two_insulated_wires",
        "coaxial_circular_strands", "coaxial_rectangular_strands")
        @testset "$case" begin
            reference = Dict(REFERENCES)[case]
            input = Engine.CoaxialInput(reference, formulation.options)
            n = length(input.radius)
            target_Z = zeros(ComplexF64, n, n, length(reference.frequencies))
            target_P = zeros(ComplexF64, n, n, length(reference.frequencies))
            for line in eachline(joinpath(
                @__DIR__, "..", "references", case, "xue2018.csv"))
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
                target_P[:, :, k] = xue_entries(
                    input, reference.earth, im * (2π * frequency), 0.0im)
            end
            @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target_Z[:, :, k]) /
                          maximum(abs(target_Z[i, i, k]) for i in 1:n)
            for k in eachindex(reference.frequencies)) <= 1.0e-8
            @test maximum(maximum(abs, result.trace.Pg[:, :, k] - target_P[:, :, k]) /
                          maximum(abs(target_P[i, i, k]) for i in 1:n)
            for k in eachindex(reference.frequencies)) <= 1.0e-8
        end
    end

    @testset "two_insulated_wires gamma" begin
        frequencies = Dict(REFERENCES)["two_insulated_wires"].frequencies
        Γ_k = im .* 1.01 .* (2π .* frequencies) ./ 299792458
        gamma = CoaxialFormulation(earth_impedance = :xue2018,
            earth_admittance = :xue2018, Γ = Γ_k)
        reference = Dict(REFERENCES)["two_insulated_wires"]
        input = Engine.CoaxialInput(reference, gamma.options)
        n = length(input.radius)
        target_Z = zeros(ComplexF64, n, n, length(reference.frequencies))
        target_P = zeros(ComplexF64, n, n, length(reference.frequencies))
        for line in eachline(joinpath(
            @__DIR__, "..", "references", "two_insulated_wires", "xue2018-gamma.csv"))
            fields = split(line, ',')
            startswith(line, "analytical") || continue
            k = searchsortedfirst(reference.frequencies, parse(Float64, fields[2]))
            target = fields[3] == "Zg" ? target_Z : target_P
            target[parse(Int, fields[4]), parse(Int, fields[5]), k] = complex(
                parse(Float64, fields[6]), parse(Float64, fields[7]))
        end
        result = compute(reference, gamma; trace = true)
        @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target_Z[:, :, k]) /
                      maximum(abs(target_Z[i, i, k]) for i in 1:n)
        for k in eachindex(reference.frequencies)) <= 1.0e-9
        @test maximum(maximum(abs, result.trace.Pg[:, :, k] - target_P[:, :, k]) /
                      maximum(abs(target_P[i, i, k]) for i in 1:n)
        for k in eachindex(reference.frequencies)) <= 1.0e-9
        @test (@inferred compute(reference, gamma)) isa Engine.LineParameters
        for (k, frequency) in enumerate(reference.frequencies)
            target_Z[:, :, k],
            target_P[:, :, k] = xue_entries(
                input, reference.earth, im * (2π * frequency), Γ_k[k])
        end
        @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target_Z[:, :, k]) /
                      maximum(abs(target_Z[i, i, k]) for i in 1:n)
        for k in eachindex(reference.frequencies)) <= 1.0e-8
        @test maximum(maximum(abs, result.trace.Pg[:, :, k] - target_P[:, :, k]) /
                      maximum(abs(target_P[i, i, k]) for i in 1:n)
        for k in eachindex(reference.frequencies)) <= 1.0e-8
    end

    result = compute(Dict(REFERENCES)["two_bare_wires"],
        CoaxialFormulation(earth_impedance = :xue2018, earth_admittance = :xue2018, Γ = 1e-3im))
    @test all(isfinite, result.Z) && all(isfinite, result.Y)
    @test_throws "treats no pair of a receiver in the air and a source in the air" compute(
        Dict(REFERENCES)["overhead_bare_wires"], formulation)
    @test_throws "treats no pair of a receiver in the air and a source in the earth" compute(
        Dict(REFERENCES)["mixed_bare_wires"], CoaxialFormulation(earth_impedance = (
            air = :carson1926, earth = :saad1996, mixed = :xue2018)))
    @test_throws r":xue2018.*earth layer 1 has mu_r = 2.0" compute(
        LineParametersProblem(Dict(REFERENCES)["two_bare_wires"].system,
            homogeneous(rho = 0.1, mu_r = 2.0);
            frequencies = Dict(REFERENCES)["two_bare_wires"].frequencies), formulation)
    @test_throws ArgumentError CoaxialFormulation(
        earth_impedance = :xue2018 => (integral = (rtol = -1.0,),))

    frequencies = [0.1, 10.0, 100.0]
    Z_delta = Array{ComplexF64, 3}[]
    P_delta = Array{ComplexF64, 3}[]
    for radius in (0.0425, 0.00425)
        copper = Material(:conductor, 1.7241e-8, 1, 0.999994, 20, 0.00393)
        design = build(CableDesign, "bare_wire",
            Terminal(:core, Region(:core_metal, Circle(radius), copper)))
        system = build(CableSystem, [design, design],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_bare_wires", line_length = 1.0)
        reference = LineParametersProblem(system, homogeneous(rho = 0.1); frequencies)
        candidate = compute(reference, formulation; trace = true)
        baseline = compute(reference, CoaxialFormulation(); trace = true)
        push!(Z_delta, candidate.trace.Zg - baseline.trace.Zg)
        push!(P_delta, candidate.trace.Pg - baseline.trace.Pg)
    end
    @test all(50 <= maximum(abs, Z_delta[1][:, :, k]) /
                    maximum(abs, Z_delta[2][:, :, k]) <= 110
    for k in eachindex(frequencies))
    @test all(50 <= maximum(abs, P_delta[1][:, :, k]) /
                    maximum(abs, P_delta[2][:, :, k]) <= 110
    for k in eachindex(frequencies))

    reference = Dict(REFERENCES)["two_bare_wires"]
    input = Engine.CoaxialInput(reference, formulation.options)
    buffers = Engine.initialize_buffers(formulation, reference, input)
    media = reference.earth.layers
    s = complex(0.0, 2π * 50.0)
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
    buffers.gamma[1] = 1.0e-3im
    formula = formulation.earth_impedance.earth
    formula(Val(:Z), Val(:self), 1, buffers, input, media, s, 1)
    formula(Val(:Z), Val(:mutual), 1, 2, buffers, input, media, s, 1)
    @test only(map((f, b, c, m, z) -> @allocated(f(Val(:Z), Val(:self), 1, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
    @test only(map(
        (f, b, c, m, z) -> @allocated(f(Val(:Z), Val(:mutual), 1, 2, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
end
