@testitem "BaseParams: calc_tubular_resistance unit tests" setup = [defaults] begin
    # Basic Functionality
    @testset "Basic Functionality" begin
        # Example from docstring
        radius_in = 0.01
        radius_ext = 0.02
        rho = 1.7241e-8
        alpha = 0.00393
        T0 = 20.0
        Top = 25.0
        expected = calc_temperature_correction(alpha, Top, T0) * rho / (π * (radius_ext^2 - radius_in^2))
        R = calc_tubular_resistance(radius_in, radius_ext, rho, alpha, T0, Top)
        @test isapprox(R, expected, atol=TEST_TOL)
    end

    # Edge Cases
    @testset "Edge Cases" begin
        # Zero thickness (radius_in == radius_ext): cross-section = 0, expect Inf or error
        r_in = 0.01
        r_ext = 0.01
        rho = 1.7241e-8
        alpha = 0.00393
        T0 = 20.0
        Top = 25.0
        # Should return Inf (division by zero)
        R = calc_tubular_resistance(r_in, r_ext, rho, alpha, T0, Top)
        @test isinf(R)
        # Very thin tube (radius_ext - radius_in ≈ eps)
        r_in2 = 0.01
        r_ext2 = 0.01 + eps()
        R2 = calc_tubular_resistance(r_in2, r_ext2, rho, alpha, T0, Top)
        @test R2 > 0
        # Large radii
        R3 = calc_tubular_resistance(1.0, 2.0, rho, alpha, T0, Top)
        @test R3 < 1e-8
        # Negative temperature coefficient (mathematically valid)
        R4 = calc_tubular_resistance(0.01, 0.02, rho, -0.001, T0, Top)
        expected4 = calc_temperature_correction(-0.001, Top, T0) * rho / (π * (0.02^2 - 0.01^2))
        @test isapprox(R4, expected4, atol=TEST_TOL)
    end

    # Numerical Consistency
    @testset "Numerical Consistency" begin
        # Float32
        Rf = calc_tubular_resistance(Float32(0.01), Float32(0.02), Float32(1.7241e-8), Float32(0.00393), Float32(20.0), Float32(25.0))
        expectedf = calc_temperature_correction(Float32(0.00393), Float32(25.0), Float32(20.0)) * Float32(1.7241e-8) / (π * (Float32(0.02)^2 - Float32(0.01)^2))
        @test isapprox(Rf, expectedf, atol=Float32(TEST_TOL))
    end

    # Physical Behavior
    @testset "Physical Behavior" begin
        rho = 1.7241e-8
        alpha = 0.00393
        T0 = 20.0
        Top = 25.0
        # Resistance decreases as cross-section increases
        R_small = calc_tubular_resistance(0.01, 0.015, rho, alpha, T0, Top)
        R_large = calc_tubular_resistance(0.01, 0.03, rho, alpha, T0, Top)
        @test R_large < R_small
        # Resistance increases with increasing resistivity
        R_lowrho = calc_tubular_resistance(0.01, 0.02, 1e-8, alpha, T0, Top)
        R_highrho = calc_tubular_resistance(0.01, 0.02, 1e-7, alpha, T0, Top)
        @test R_highrho > R_lowrho
        # Resistance increases with increasing temperature (for positive alpha)
        R_T1 = calc_tubular_resistance(0.01, 0.02, rho, alpha, T0, 25.0)
        R_T2 = calc_tubular_resistance(0.01, 0.02, rho, alpha, T0, 75.0)
        @test R_T2 > R_T1
    end

    # Type Stability & Promotion
    @testset "Type Stability & Promotion" begin
        # All Float64
        Rf = calc_tubular_resistance(0.01, 0.02, 1.7241e-8, 0.00393, 20.0, 25.0)
        @test typeof(Rf) == Float64
        # All Measurement
        using Measurements
        rin_m = measurement(0.01, 1e-6)
        rext_m = measurement(0.02, 1e-6)
        rho_m = measurement(1.7241e-8, 1e-10)
        alpha_m = measurement(0.00393, 1e-5)
        T0_m = measurement(20.0, 0.1)
        Top_m = measurement(25.0, 0.1)
        Rm = calc_tubular_resistance(rin_m, rext_m, rho_m, alpha_m, T0_m, Top_m)
        @test Rm isa Measurement{Float64}
        # Mixed: first argument as Measurement
        Rmix1 = calc_tubular_resistance(rin_m, 0.02, 1.7241e-8, 0.00393, 20.0, 25.0)
        @test Rmix1 isa Measurement{Float64}
        # Mixed: middle argument as Measurement
        Rmix2 = calc_tubular_resistance(0.01, 0.02, rho_m, 0.00393, 20.0, 25.0)
        @test Rmix2 isa Measurement{Float64}
        # Mixed: last argument as Measurement
        Rmix3 = calc_tubular_resistance(0.01, 0.02, 1.7241e-8, 0.00393, 20.0, Top_m)
        @test Rmix3 isa Measurement{Float64}
    end

    # Uncertainty Quantification
    @testset "Uncertainty Quantification" begin
        rin_m = measurement(0.01, 1e-6)
        rext_m = measurement(0.02, 1e-6)
        rho_m = measurement(1.7241e-8, 1e-10)
        alpha_m = measurement(0.00393, 1e-5)
        T0_m = measurement(20.0, 0.1)
        Top_m = measurement(25.0, 0.1)
        Rm = calc_tubular_resistance(rin_m, rext_m, rho_m, alpha_m, T0_m, Top_m)
        # Analytical propagation (approximate, neglecting correlations):
        ΔA = π * (value(rext_m)^2 - value(rin_m)^2)
        k = value(calc_temperature_correction(alpha_m, Top_m, T0_m))
        μ = k * value(rho_m) / ΔA
        @test isapprox(value(Rm), μ, atol=TEST_TOL)
        # Uncertainty should be nonzero and scale with input uncertainties
        @test uncertainty(Rm) > 0
    end
end