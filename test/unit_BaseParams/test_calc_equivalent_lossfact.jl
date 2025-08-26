@testitem "BaseParams: calc_equivalent_lossfact unit tests" setup = [defaults] begin
    @testset "Basic Functionality" begin
        # Example from docstring: G_eq=1e-8 S·m, C_eq=1e-10 F/m, ω=2π*50
        G_eq = 1e-8
        C_eq = 1e-10
        ω = 2 * pi * 50
        result = calc_equivalent_lossfact(G_eq, C_eq, ω)
        expected = G_eq / (ω * C_eq)
        @test isapprox(result, expected; atol=TEST_TOL)
        @test result > 0
    end

    @testset "Edge Cases" begin
        # Zero conductance
        result = calc_equivalent_lossfact(0.0, 1e-10, 2 * pi * 50)
        @test isapprox(result, 0.0; atol=TEST_TOL)
        # Zero capacitance (should be Inf)
        result = calc_equivalent_lossfact(1e-8, 0.0, 2 * pi * 50)
        @test isinf(result)
        # Zero frequency (should be Inf)
        result = calc_equivalent_lossfact(1e-8, 1e-10, 0.0)
        @test isinf(result)
        # Inf/NaN
        @test isnan(calc_equivalent_lossfact(NaN, 1e-10, 2 * pi * 50))
        @test isnan(calc_equivalent_lossfact(1e-8, NaN, 2 * pi * 50))
        @test isnan(calc_equivalent_lossfact(1e-8, 1e-10, NaN))
        @test isinf(calc_equivalent_lossfact(Inf, 1e-10, 2 * pi * 50))
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        r = calc_equivalent_lossfact(Float32(1e-8), Float32(1e-10), Float32(2 * pi * 50))
        d = calc_equivalent_lossfact(1e-8, 1e-10, 2 * pi * 50)
        @test isapprox(r, d; atol=1e-6)
    end

    @testset "Physical Behavior" begin
        # Increases with G_eq
        r1 = calc_equivalent_lossfact(1e-8, 1e-10, 2 * pi * 50)
        r2 = calc_equivalent_lossfact(2e-8, 1e-10, 2 * pi * 50)
        @test r2 > r1
        # Decreases with C_eq
        r3 = calc_equivalent_lossfact(1e-8, 2e-10, 2 * pi * 50)
        @test r3 < r1
        # Decreases with ω
        r4 = calc_equivalent_lossfact(1e-8, 1e-10, 2 * pi * 100)
        @test r4 < r1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        # All Float64
        r1 = calc_equivalent_lossfact(1e-8, 1e-10, 2 * pi * 50)
        @test typeof(r1) == Float64
        # All Measurement
        r2 = calc_equivalent_lossfact(measurement(1e-8, 1e-10), measurement(1e-10, 1e-12), measurement(2 * pi * 50, 0.1))
        @test r2 isa Measurement{Float64}
        # Mixed: G_eq as Measurement
        r3 = calc_equivalent_lossfact(measurement(1e-8, 1e-10), 1e-10, 2 * pi * 50)
        @test r3 isa Measurement{Float64}
        # Mixed: C_eq as Measurement
        r4 = calc_equivalent_lossfact(1e-8, measurement(1e-10, 1e-12), 2 * pi * 50)
        @test r4 isa Measurement{Float64}
        # Mixed: ω as Measurement
        r5 = calc_equivalent_lossfact(1e-8, 1e-10, measurement(2 * pi * 50, 0.1))
        @test r5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        G_eq = measurement(1e-8, 1e-10)
        C_eq = measurement(1e-10, 1e-12)
        ω = measurement(2 * pi * 50, 0.1)
        result = calc_equivalent_lossfact(G_eq, C_eq, ω)
        @test result isa Measurement{Float64}
        @test uncertainty(result) > 0
    end
end
