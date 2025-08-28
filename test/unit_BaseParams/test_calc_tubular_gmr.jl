@testitem "BaseParams: calc_tubular_gmr unit tests" setup = [defaults] begin
    using Measurements: measurement, value, uncertainty

    @testset "Basic Functionality" begin
        # Example from docstring
        radius_ext = 0.02
        radius_in = 0.01
        mu_r = 1.0
        gmr = calc_tubular_gmr(radius_ext, radius_in, mu_r)
        # Manual calculation for expected value
        term1 = (radius_in^4 / (radius_ext^2 - radius_in^2)^2) * log(radius_ext / radius_in)
        term2 = (3 * radius_in^2 - radius_ext^2) / (4 * (radius_ext^2 - radius_in^2))
        Lin = (μ₀ * mu_r / (2 * π)) * (term1 - term2)
        expected = exp(log(radius_ext) - (2 * π / μ₀) * Lin)
        @test isapprox(gmr, expected; atol=TEST_TOL)
        @test gmr > 0
        @test_throws ArgumentError calc_tubular_gmr(radius_in, radius_ext, mu_r)
        @test_throws ArgumentError calc_tubular_gmr(0.0, radius_in, mu_r)
    end

    @testset "Edge Cases" begin
        # Thin shell: radius_ext ≈ radius_in
        radius_ext = 0.01
        radius_in = 0.01
        mu_r = 1.0
        gmr = calc_tubular_gmr(radius_ext, radius_in, mu_r)
        @test isapprox(gmr, radius_ext; atol=TEST_TOL)

        # Infinitely thick tube: radius_in ≫ 0, radius_in / radius_ext ≈ 0
        radius_ext = 1.0
        radius_in = 1e-12
        mu_r = 1.0
        gmr = calc_tubular_gmr(radius_ext, radius_in, mu_r)
        @test isapprox(gmr, 0.7788; atol=1e-4)

        # radius_in = 0 (solid cylinder)
        radius_ext = 0.02
        radius_in = 0.0
        mu_r = 1.0
        gmr = calc_tubular_gmr(radius_ext, radius_in, mu_r)
        @test isapprox(gmr, 0.7788 * radius_ext; atol=1e-4)

        # radius_ext < radius_in (should throw)
        radius_ext = 0.01
        radius_in = 0.02
        mu_r = 1.0
        @test_throws ArgumentError calc_tubular_gmr(radius_ext, radius_in, mu_r)
    end

    @testset "Numerical Consistency" begin
        # Float64
        gmr1 = calc_tubular_gmr(0.02, 0.01, 1.0)
        # Measurement{Float64}
        gmr2 = calc_tubular_gmr(measurement(0.02, 1e-4), measurement(0.01, 1e-4), measurement(1.0, 0.01))
        @test isapprox(value(gmr2), gmr1; atol=TEST_TOL)
        @test uncertainty(gmr2) > 0
    end

    @testset "Physical Behavior" begin
        # GMR increases with radius_ext
        gmr1 = calc_tubular_gmr(0.01, 0.005, 1.0)
        gmr2 = calc_tubular_gmr(0.02, 0.005, 1.0)
        @test gmr2 > gmr1
        # GMR decreases with mu_r
        gmr1 = calc_tubular_gmr(0.02, 0.01, 0.5)
        gmr2 = calc_tubular_gmr(0.02, 0.01, 2.0)
        @test gmr2 < gmr1
    end

    @testset "Type Stability & Promotion" begin
        # All Float64
        gmr = calc_tubular_gmr(0.02, 0.01, 1.0)
        @test typeof(gmr) == Float64
        # All Measurement
        gmr = calc_tubular_gmr(measurement(0.02, 1e-4), measurement(0.01, 1e-4), measurement(1.0, 0.01))
        @test gmr isa Measurement{Float64}
        # Mixed: radius_ext as Measurement
        gmr = calc_tubular_gmr(measurement(0.02, 1e-4), 0.01, 1.0)
        @test gmr isa Measurement{Float64}
        # Mixed: radius_in as Measurement
        gmr = calc_tubular_gmr(0.02, measurement(0.01, 1e-4), 1.0)
        @test gmr isa Measurement{Float64}
        # Mixed: mu_r as Measurement
        gmr = calc_tubular_gmr(0.02, 0.01, measurement(1.0, 0.01))
        @test gmr isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        radius_ext = measurement(0.02, 1e-4)
        radius_in = measurement(0.01, 1e-4)
        mu_r = measurement(1.0, 0.01)
        gmr = calc_tubular_gmr(radius_ext, radius_in, mu_r)
        @test gmr isa Measurement{Float64}
        @test uncertainty(gmr) > 0
    end
end
