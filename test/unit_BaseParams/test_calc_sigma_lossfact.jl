@testitem "BaseParams: calc_sigma_lossfact unit tests" setup = [defaults] begin
    @testset "Basic Functionality" begin
        # Example from docstring: G_eq=2.7169e-9 S·m, r_in=0.01 m, r_ext=0.02 m
        G_eq = 2.7169e-9
        r_in = 0.01
        r_ext = 0.02
        result = calc_sigma_lossfact(G_eq, r_in, r_ext)
        expected = G_eq * log(r_ext / r_in) / (2 * pi)
        @test isapprox(result, expected; atol=TEST_TOL)
        @test result > 0
    end

    @testset "Edge Cases" begin
        # Zero conductance
        result = calc_sigma_lossfact(0.0, 0.01, 0.02)
        @test isapprox(result, 0.0; atol=TEST_TOL)
        # Collapsing geometry: r_ext == r_in
        result = calc_sigma_lossfact(1e-9, 0.01, 0.01)
        @test isapprox(result, 0.0; atol=TEST_TOL)
        # Very large radii
        result = calc_sigma_lossfact(1e-9, 1e3, 1e6)
        expected = 1e-9 * log(1e6 / 1e3) / (2 * pi)
        @test isapprox(result, expected; atol=TEST_TOL)
        # Inf/NaN
        @test isnan(calc_sigma_lossfact(NaN, 0.01, 0.02))
        @test isnan(calc_sigma_lossfact(1e-9, NaN, 0.02))
        @test isnan(calc_sigma_lossfact(1e-9, 0.01, NaN))
        @test isinf(calc_sigma_lossfact(Inf, 0.01, 0.02))
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        r = calc_sigma_lossfact(Float32(1e-9), Float32(0.01), Float32(0.02))
        d = calc_sigma_lossfact(1e-9, 0.01, 0.02)
        @test isapprox(r, d; atol=TEST_TOL)
    end

    @testset "Physical Behavior" begin
        # Increases with G_eq
        r1 = calc_sigma_lossfact(1e-9, 0.01, 0.02)
        r2 = calc_sigma_lossfact(2e-9, 0.01, 0.02)
        @test r2 > r1
        # Increases with log(r_ext/r_in)
        r3 = calc_sigma_lossfact(1e-9, 0.01, 0.04)
        @test r3 > r1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        # All Float64
        r1 = calc_sigma_lossfact(1e-9, 0.01, 0.02)
        @test typeof(r1) == Float64
        # All Measurement
        r2 = calc_sigma_lossfact(measurement(1e-9, 1e-11), measurement(0.01, 1e-5), measurement(0.02, 1e-5))
        @test r2 isa Measurement{Float64}
        # Mixed: G_eq as Measurement
        r3 = calc_sigma_lossfact(measurement(1e-9, 1e-11), 0.01, 0.02)
        @test r3 isa Measurement{Float64}
        # Mixed: radius_in as Measurement
        r4 = calc_sigma_lossfact(1e-9, measurement(0.01, 1e-5), 0.02)
        @test r4 isa Measurement{Float64}
        # Mixed: radius_ext as Measurement
        r5 = calc_sigma_lossfact(1e-9, 0.01, measurement(0.02, 1e-5))
        @test r5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        G_eq = measurement(1e-9, 1e-11)
        r_in = measurement(0.01, 1e-5)
        r_ext = measurement(0.02, 1e-5)
        result = calc_sigma_lossfact(G_eq, r_in, r_ext)
        @test result isa Measurement{Float64}
        @test uncertainty(result) > 0
    end
end
