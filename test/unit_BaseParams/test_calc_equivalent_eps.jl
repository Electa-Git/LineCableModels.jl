@testitem "BaseParams: calc_equivalent_eps unit tests" setup = [defaults] begin
    @testset "Basic Functionality" begin
        # Example from docstring: C_eq=1e-10 F/m, r_ext=0.01 m, r_in=0.005 m
        result = calc_equivalent_eps(1e-10, 0.01, 0.005)
        expected = (1e-10 * log(0.01 / 0.005)) / (2 * pi) / ε₀
        @test isapprox(result, expected; atol=TEST_TOL)
        @test result > 0
    end

    @testset "Edge Cases" begin
        # Zero capacitance
        result = calc_equivalent_eps(0.0, 0.01, 0.005)
        @test isapprox(result, 0.0; atol=TEST_TOL)
        # Collapsing geometry: r_ext == r_in
        result = calc_equivalent_eps(1e-10, 0.01, 0.01)
        @test isapprox(result, 0.0; atol=TEST_TOL)
        # Very large radii
        result = calc_equivalent_eps(1e-10, 1e6, 1e3)
        expected = (1e-10 * log(1e6 / 1e3)) / (2 * pi) / ε₀
        @test isapprox(result, expected; atol=TEST_TOL)
        # Inf/NaN
        @test isnan(calc_equivalent_eps(NaN, 0.01, 0.005))
        @test isnan(calc_equivalent_eps(1e-10, NaN, 0.005))
        @test isnan(calc_equivalent_eps(1e-10, 0.01, NaN))
        @test isinf(calc_equivalent_eps(Inf, 0.01, 0.005))
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        r = calc_equivalent_eps(Float32(1e-10), Float32(0.01), Float32(0.005))
        d = calc_equivalent_eps(1e-10, 0.01, 0.005)
        @test isapprox(r, d; atol=1e-6)
    end

    @testset "Physical Behavior" begin
        # Increases with capacitance
        r1 = calc_equivalent_eps(1e-10, 0.01, 0.005)
        r2 = calc_equivalent_eps(2e-10, 0.01, 0.005)
        @test r2 > r1
        # Increases with log(r_ext/r_in)
        r3 = calc_equivalent_eps(1e-10, 0.02, 0.005)
        @test r3 > r1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        # All Float64
        r1 = calc_equivalent_eps(1e-10, 0.01, 0.005)
        @test typeof(r1) == Float64
        # All Measurement
        r2 = calc_equivalent_eps(measurement(1e-10, 1e-12), measurement(0.01, 1e-5), measurement(0.005, 1e-5))
        @test r2 isa Measurement{Float64}
        # Mixed: C_eq as Measurement
        r3 = calc_equivalent_eps(measurement(1e-10, 1e-12), 0.01, 0.005)
        @test r3 isa Measurement{Float64}
        # Mixed: radius_ext as Measurement
        r4 = calc_equivalent_eps(1e-10, measurement(0.01, 1e-5), 0.005)
        @test r4 isa Measurement{Float64}
        # Mixed: radius_in as Measurement
        r5 = calc_equivalent_eps(1e-10, 0.01, measurement(0.005, 1e-5))
        @test r5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        C_eq = measurement(1e-10, 1e-12)
        r_ext = measurement(0.01, 1e-5)
        r_in = measurement(0.005, 1e-5)
        result = calc_equivalent_eps(C_eq, r_ext, r_in)
        @test result isa Measurement{Float64}
        @test uncertainty(result) > 0
    end
end
