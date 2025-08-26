@testitem "BaseParams: calc_tubular_inductance unit tests" setup = [defaults] begin
    using Measurements
    # Basic Functionality
    @testset "Basic Functionality" begin
        # Example from docstring: radius_in = 0.01, radius_ext = 0.02, mu_r = 1.0
        L = calc_tubular_inductance(0.01, 0.02, 1.0)
        expected = 1.0 * μ₀ / (2 * π) * log(0.02 / 0.01)
        @test isapprox(L, expected, atol=TEST_TOL)
    end

    # Edge Cases
    @testset "Edge Cases" begin
        # Very thin tube (radius_ext ≈ radius_in)
        r_in = 0.01
        r_ext = 0.010001
        L_thin = calc_tubular_inductance(r_in, r_ext, 1.0)
        expected_thin = 1.0 * μ₀ / (2 * π) * log(r_ext / r_in)
        @test isapprox(L_thin, expected_thin, atol=TEST_TOL)
        # Large radii
        L_large = calc_tubular_inductance(1e3, 2e3, 1.0)
        expected_large = 1.0 * μ₀ / (2 * π) * log(2e3 / 1e3)
        @test isapprox(L_large, expected_large, atol=TEST_TOL)
        # mu_r = 0 (non-magnetic)
        @test isapprox(calc_tubular_inductance(0.01, 0.02, 0.0), 0.0, atol=TEST_TOL)
    end

    # Numerical Consistency
    @testset "Numerical Consistency" begin
        # Float32
        Lf = calc_tubular_inductance(Float32(0.01), Float32(0.02), Float32(1.0))
        expectedf = Float32(μ₀) / (2f0 * Float32(π)) * log(Float32(0.02) / Float32(0.01))
        @test isapprox(Lf, expectedf, atol=Float32(TEST_TOL))
        # Rational
        Lr = calc_tubular_inductance(1 // 100, 1 // 50, 1 // 1)
        expectedr = (1 // 1) * μ₀ / (2 * π) * log((1 // 50) / (1 // 100))
        @test isapprox(Lr, expectedr, atol=TEST_TOL)
    end

    # Physical Behavior
    @testset "Physical Behavior" begin
        # L increases with mu_r
        L1 = calc_tubular_inductance(0.01, 0.02, 1.0)
        L2 = calc_tubular_inductance(0.01, 0.02, 2.0)
        @test L2 > L1
        # L increases with radius_ext
        L3 = calc_tubular_inductance(0.01, 0.03, 1.0)
        @test L3 > L1
        # L decreases with radius_in
        L4 = calc_tubular_inductance(0.02, 0.03, 1.0)
        @test L4 < L3
    end

    # Type Stability & Promotion
    @testset "Type Stability & Promotion" begin
        # All Float64
        Lf = calc_tubular_inductance(0.01, 0.02, 1.0)
        @test typeof(Lf) == Float64
        # All Measurement
        rinm = measurement(0.01, 1e-5)
        rextm = measurement(0.02, 1e-5)
        murm = measurement(1.0, 1e-3)
        Lm = calc_tubular_inductance(rinm, rextm, murm)
        @test Lm isa Measurement{Float64}
        # Mixed: radius_in as Measurement
        Lmix1 = calc_tubular_inductance(rinm, 0.02, 1.0)
        @test Lmix1 isa Measurement{Float64}
        # Mixed: radius_ext as Measurement
        Lmix2 = calc_tubular_inductance(0.01, rextm, 1.0)
        @test Lmix2 isa Measurement{Float64}
        # Mixed: mu_r as Measurement
        Lmix3 = calc_tubular_inductance(0.01, 0.02, murm)
        @test Lmix3 isa Measurement{Float64}
    end

    # Uncertainty Quantification
    @testset "Uncertainty Quantification" begin
        rinm = measurement(0.01, 1e-5)
        rextm = measurement(0.02, 1e-5)
        murm = measurement(1.0, 1e-3)
        Lm = calc_tubular_inductance(rinm, rextm, murm)
        # Analytical propagation: L = mu_r * μ₀ / (2π) * log(r_ext/r_in)
        μ = 1.0 * μ₀ / (2 * π) * log(0.02 / 0.01)
        # Partial derivatives
        dL_drin = -murm * μ₀ / (2 * π) * (1 / rinm) / (rextm / rinm)
        dL_drext = murm * μ₀ / (2 * π) * (1 / rextm) / (rextm / rinm)
        dL_dmurm = μ₀ / (2 * π) * log(0.02 / 0.01)
        σ2 = (value(dL_drin) * uncertainty(rinm))^2 + (value(dL_drext) * uncertainty(rextm))^2 + (value(dL_dmurm) * uncertainty(murm))^2
        @test isapprox(value(Lm), μ, atol=TEST_TOL)
        @test isapprox(uncertainty(Lm), sqrt(σ2), atol=TEST_TOL)
    end
end
