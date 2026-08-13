@testitem "BaseParams: calc_temperature_correction unit tests" setup = [defaults] begin
    using Measurements
    # Basic Functionality
    @testset "Basic Functionality" begin
        # Example from docstring: alpha = 0.00393, Top = 75.0, T0 = 20.0
        k = calc_temperature_correction(0.00393, 75.0, 20.0)
        @test isapprox(k, 1.2161, atol=1e-4)

        # Default T0 (should use T₀ constant)
        k2 = calc_temperature_correction(0.00393, 75.0)
        k2_ref = calc_temperature_correction(0.00393, 75.0, T₀)
        @test isapprox(k2, k2_ref, atol=TEST_TOL)
    end

    # Edge Cases
    @testset "Edge Cases" begin
        # Zero temperature difference
        @test isapprox(calc_temperature_correction(0.00393, 20.0, 20.0), 1.0, atol=TEST_TOL)
        # Negative alpha (unusual, but mathematically valid)
        @test isapprox(calc_temperature_correction(-0.001, 30.0, 20.0), 0.99, atol=TEST_TOL)
        # Large temperature difference within ΔTmax
        @test isapprox(calc_temperature_correction(0.00393, 20.0 + (ΔTmax - 1), 20.0), 1 + 0.00393 * (ΔTmax - 1), atol=TEST_TOL)
    end

    # Numerical Consistency
    @testset "Numerical Consistency" begin
        # Float32
        kf = calc_temperature_correction(Float32(0.00393), Float32(75.0), Float32(20.0))
        @test isapprox(kf, 1.2161f0, atol=Float32(1e-4))
    end

    # Physical Behavior
    @testset "Physical Behavior" begin
        # Correction increases with temperature for positive alpha
        k1 = calc_temperature_correction(0.00393, 50.0, 20.0)
        k2 = calc_temperature_correction(0.00393, 80.0, 20.0)
        @test k2 > k1
        # Correction decreases with temperature for negative alpha
        k3 = calc_temperature_correction(-0.001, 50.0, 20.0)
        k4 = calc_temperature_correction(-0.001, 80.0, 20.0)
        @test k4 < k3
    end

    # Type Stability & Promotion
    @testset "Type Stability & Promotion" begin
        # All Float64
        kf = calc_temperature_correction(0.00393, 75.0, 20.0)
        @test typeof(kf) == Float64
        # All Measurement
        αm = measurement(0.00393, 1e-5)
        Topm = measurement(75.0, 0.1)
        T0m = measurement(20.0, 0.1)
        km = calc_temperature_correction(αm, Topm, T0m)
        @test km isa Measurement{Float64}
        # Mixed: alpha as Measurement
        kmix1 = calc_temperature_correction(αm, 75.0, 20.0)
        @test kmix1 isa Measurement{Float64}
        # Mixed: Top as Measurement
        kmix2 = calc_temperature_correction(0.00393, Topm, 20.0)
        @test kmix2 isa Measurement{Float64}
        # Mixed: T0 as Measurement
        kmix3 = calc_temperature_correction(0.00393, 75.0, T0m)
        @test kmix3 isa Measurement{Float64}
    end

    # Uncertainty Quantification
    @testset "Uncertainty Quantification" begin
        αm = measurement(0.00393, 1e-5)
        Topm = measurement(75.0, 0.1)
        T0m = measurement(20.0, 0.1)
        km = calc_temperature_correction(αm, Topm, T0m)
        # Analytical propagation: k = 1 + α*(Top-T0)
        # σ² = (Top-T0)²*σ_α² + α²*σ_Top² + α²*σ_T0²
        μ = 1 + 0.00393 * (75.0 - 20.0)
        σ2 = (75.0 - 20.0)^2 * 1e-5^2 + 0.00393^2 * 0.1^2 + 0.00393^2 * 0.1^2
        @test isapprox(value(km), μ, atol=TEST_TOL)
        @test isapprox(uncertainty(km), sqrt(σ2), atol=TEST_TOL)
    end
end
