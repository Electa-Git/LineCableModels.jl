using QuadGK: quadgk
using SpecialFunctions: besselk

function pollaczek_integral(H, y, γ²)
    return first(quadgk(
        λ -> exp(-sqrt(λ^2 + γ²) * H) * cos(y * λ) / (λ + sqrt(λ^2 + γ²)),
        sort([0.0, abs(sqrt(γ²)) / 4, abs(sqrt(γ²)), 4abs(sqrt(γ²)),
            1 / (4H), 1 / H, 4 / H, Inf])...;
        rtol = 1.0e-13, atol = 0.0, maxevals = 10^8, order = 11))
end

function pollaczek_mixed(h_a, h_g, y, γ²)
    H = h_a + h_g
    return first(quadgk(
        λ -> exp(-h_a * λ - h_g * sqrt(λ^2 + γ²)) * cos(y * λ) /
             (λ + sqrt(λ^2 + γ²)),
        sort([0.0, abs(sqrt(γ²)) / 4, abs(sqrt(γ²)), 4abs(sqrt(γ²)),
            1 / (4H), 1 / H, 4 / H, Inf])...;
        rtol = 1.0e-13, atol = 0.0, maxevals = 10^8, order = 11))
end

@testset "Engine pollaczek1926 earth" begin
    Engine = LineCableModels.Engine
    μ₀ = 4π * 1e-7
    formulation = CoaxialFormulation(earth_impedance = :pollaczek1926)
    classical = CoaxialFormulation(earth_impedance = (
        air = :carson1926, earth = :pollaczek1926, mixed = :pollaczek1926))
    for case in ("two_bare_wires", "two_insulated_wires",
        "coaxial_circular_strands", "coaxial_rectangular_strands")
        @testset "$case" begin
            reference = Dict(REFERENCES)[case]
            input = Engine.CoaxialInput(reference, formulation.options)
            n = length(input.radius)
            target = zeros(ComplexF64, n, n, length(reference.frequencies))
            for line in eachline(joinpath(
                @__DIR__, "..", "references", case, "pollaczek1926.csv"))
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
            @test (@inferred compute(reference, formulation)) isa Engine.LineParameters
            for (k, frequency) in enumerate(reference.frequencies)
                s = im * (2π * frequency)
                γ² = s * μ₀ / reference.earth.layers[2].rho
                γ = sqrt(γ²)
                for i in 1:n, j in 1:n

                    h_i = abs(input.y[i])
                    h_j = abs(input.y[j])
                    H = h_i + h_j
                    y = abs(input.x[i] - input.x[j])
                    D = hypot(y, h_i - h_j)
                    D′ = hypot(y, H)
                    target[i, j, k] = s * μ₀ / (2π) *
                                      (i == j ?
                                       besselk(0, γ * input.radius[i]) -
                                       besselk(0, 2γ * h_i) +
                                       2pollaczek_integral(2h_i, 0.0, γ²) :
                                       besselk(0, γ * D) - besselk(0, γ * D′) +
                                       2pollaczek_integral(H, y, γ²))
                end
            end
            @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target[:, :, k]) /
                          maximum(abs(target[i, i, k]) for i in 1:n)
            for k in eachindex(reference.frequencies)) <= 1.0e-8
        end
    end

    reference = Dict(REFERENCES)["mixed_bare_wires"]
    input = Engine.CoaxialInput(reference, classical.options)
    n = length(input.radius)
    target = zeros(ComplexF64, n, n, length(reference.frequencies))
    for line in eachline(joinpath(
        @__DIR__, "..", "references", "mixed_bare_wires", "pollaczek1926.csv"))
        fields = split(line, ',')
        startswith(line, "analytical") || continue
        k = searchsortedfirst(reference.frequencies, parse(Float64, fields[2]))
        target[parse(Int, fields[4]), parse(Int, fields[5]), k] = complex(
            parse(Float64, fields[6]), parse(Float64, fields[7]))
    end
    result = compute(reference, classical; trace = true)
    @test maximum(maximum(abs, result.trace.Zg[:, :, k] - target[:, :, k]) /
                  maximum(abs(target[i, i, k]) for i in 1:n)
    for k in eachindex(reference.frequencies)) <= 1.0e-9
    @test (@inferred compute(reference, classical)) isa Engine.LineParameters
    for (k, frequency) in enumerate(reference.frequencies)
        s = im * (2π * frequency)
        γ² = s * μ₀ / reference.earth.layers[2].rho
        target[1, 2, k] = target[2, 1, k] = s * μ₀ / π *
                                            pollaczek_mixed(1.0, 1.0, 1.0, γ²)
    end
    @test maximum(max(abs(result.trace.Zg[1, 2, k] - target[1, 2, k]),
                      abs(result.trace.Zg[2, 1, k] - target[2, 1, k])) /
                  maximum(abs(result.trace.Zg[i, i, k]) for i in 1:n)
    for k in eachindex(reference.frequencies)) <= 1.0e-8

    design = reference.system.designs[1]
    for placements in ([Placement(0.0, 5.0), Placement(1.0, -0.5)],
        [Placement(1.0, -0.5), Placement(0.0, 5.0)])
        system = build(CableSystem, [design, design], placements;
            connections = [(core = 1,), (core = 2,)], id = "uneven", line_length = 1.0)
        problem = LineParametersProblem(system, reference.earth;
            frequencies = [50.0, 1.0e4, 1.0e6])
        result = compute(problem, classical; trace = true)
        target = [s * μ₀ / π * pollaczek_mixed(5.0, 0.5, 1.0, s * μ₀ / 0.1)
                  for s in im .* (2π .* problem.frequencies)]
        @test maximum(max(abs(result.trace.Zg[1, 2, k] - target[k]),
                          abs(result.trace.Zg[2, 1, k] - target[k])) / abs(target[k])
        for k in eachindex(problem.frequencies)) <= 1.0e-8
    end

    @test_throws "treats no pair of a receiver in the air and a source in the air" compute(
        Dict(REFERENCES)["overhead_bare_wires"], formulation)
    @test_throws "gives no potential!" CoaxialFormulation(earth_admittance = :pollaczek1926)
    magnetic = LineParametersProblem(Dict(REFERENCES)["two_bare_wires"].system,
        homogeneous(rho = 0.1, mu_r = 2.0); frequencies = [10.0, 100.0])
    @test_logs (:warn,
        "the formula :pollaczek1926 holds for mu_r_earth in 1.0 .. 1.0; " *
        "the computation has mu_r_earth in 2.0 .. 2.0") compute(magnetic, formulation)
    gamma = CoaxialFormulation(earth_impedance = :pollaczek1926, Γ = 1e-3im)
    @test_logs (:warn,
        "the formula :pollaczek1926 holds for Γ in 0.0 .. 0.0; " *
        "the computation has Γ in 0.001 .. 0.001") LineCableModels.validate(
        LineCableModels.domain(Dict(REFERENCES)["two_bare_wires"], gamma.options), gamma)
    @test_throws ArgumentError CoaxialFormulation(
        earth_impedance = :pollaczek1926 => (integral = (rtol = -1.0,),))

    buffers = Engine.initialize_buffers(classical, reference, input)
    media = reference.earth.layers
    s = complex(0.0, 2π * 50.0)
    formula = classical.earth_impedance.mixed
    formula(Val(:Z), Val(:mutual), 1, 2, buffers, input, media, s, 1)
    @test only(map(
        (f, b, c, m, z) -> @allocated(f(Val(:Z), Val(:mutual), 1, 2, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0

    reference = Dict(REFERENCES)["two_bare_wires"]
    input = Engine.CoaxialInput(reference, formulation.options)
    buffers = Engine.initialize_buffers(formulation, reference, input)
    media = reference.earth.layers
    formula = formulation.earth_impedance.earth
    formula(Val(:Z), Val(:self), 1, buffers, input, media, s, 1)
    formula(Val(:Z), Val(:mutual), 1, 2, buffers, input, media, s, 1)
    @test only(map((f, b, c, m, z) -> @allocated(f(Val(:Z), Val(:self), 1, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
    @test only(map(
        (f, b, c, m, z) -> @allocated(f(Val(:Z), Val(:mutual), 1, 2, b, c, m, z, 1)),
        (formula,), (buffers,), (input,), (media,), (s,))) == 0
end
