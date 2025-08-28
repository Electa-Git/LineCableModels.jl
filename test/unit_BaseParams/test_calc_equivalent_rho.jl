@testitem "BaseParams: calc_equivalent_rho unit tests" setup = [defaults] begin
    @testset "Basic Functionality" begin
        # Example from docstring: R=0.01 Ω, r_ext=0.02 m, r_in=0.01 m
        result = calc_equivalent_rho(0.01, 0.02, 0.01)
        expected = 0.01 * π * (0.02^2 - 0.01^2)
        @test isapprox(result, expected; atol=TEST_TOL)
        @test result > 0
    end

    @testset "Edge Cases" begin
        # Zero resistance
        result = calc_equivalent_rho(0.0, 0.02, 0.01)
        @test isapprox(result, 0.0; atol=TEST_TOL)
        # Zero thickness (r_ext == r_in)
        result = calc_equivalent_rho(0.01, 0.01, 0.01)
        @test isapprox(result, 0.0; atol=TEST_TOL)
        # Very large radii
        result = calc_equivalent_rho(0.01, 1e6, 1e3)
        expected = 0.01 * π * (1e6^2 - 1e3^2)
        @test isapprox(result, expected; atol=TEST_TOL)
        # Inf/NaN
        @test isnan(calc_equivalent_rho(NaN, 0.02, 0.01))
        @test isnan(calc_equivalent_rho(0.01, NaN, 0.01))
        @test isnan(calc_equivalent_rho(0.01, 0.02, NaN))
        @test isinf(calc_equivalent_rho(Inf, 0.02, 0.01))
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        r = calc_equivalent_rho(Float32(0.01), Float32(0.02), Float32(0.01))
        d = calc_equivalent_rho(0.01, 0.02, 0.01)
        @test isapprox(r, d; atol=TEST_TOL)
    end

    @testset "Physical Behavior" begin
        # Increases with resistance
        r1 = calc_equivalent_rho(0.01, 0.02, 0.01)
        r2 = calc_equivalent_rho(0.02, 0.02, 0.01)
        @test r2 > r1
        # Increases with area
        r3 = calc_equivalent_rho(0.01, 0.03, 0.01)
        @test r3 > r1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        # All Float64
        r1 = calc_equivalent_rho(0.01, 0.02, 0.01)
        @test typeof(r1) == Float64
        # All Measurement
        r2 = calc_equivalent_rho(measurement(0.01, 1e-4), measurement(0.02, 1e-5), measurement(0.01, 1e-5))
        @test r2 isa Measurement{Float64}
        # Mixed: R as Measurement
        r3 = calc_equivalent_rho(measurement(0.01, 1e-4), 0.02, 0.01)
        @test r3 isa Measurement{Float64}
        # Mixed: radius_ext_con as Measurement
        r4 = calc_equivalent_rho(0.01, measurement(0.02, 1e-5), 0.01)
        @test r4 isa Measurement{Float64}
        # Mixed: radius_in_con as Measurement
        r5 = calc_equivalent_rho(0.01, 0.02, measurement(0.01, 1e-5))
        @test r5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        R = measurement(0.01, 1e-4)
        r_ext = measurement(0.02, 1e-5)
        r_in = measurement(0.01, 1e-5)
        result = calc_equivalent_rho(R, r_ext, r_in)
        @test result isa Measurement{Float64}
        @test uncertainty(result) > 0
    end
end
