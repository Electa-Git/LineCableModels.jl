@testitem "BaseParams / calc_temperature_correction / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    using Measurements
    # Basic Functionality
    @testset "Basic Functionality" begin
        # Example from docstring: alpha = 0.00393, Top = 75.0, T0 = 20.0
        k = calc_temperature_correction(0.00393, 75.0, 20.0)
        @test isapprox(k, 1.2161, atol = 1e-4)

        # Default T0 (should use T₀ constant)
        k2 = calc_temperature_correction(0.00393, 75.0)
        k2_ref = calc_temperature_correction(0.00393, 75.0, T₀)
        @test isapprox(k2, k2_ref, atol = TestNumerics.absolute_floor(Float64))
    end

    # Edge Cases
    @testset "Edge Cases" begin
        # Zero temperature difference
        @test isapprox(calc_temperature_correction(0.00393, 20.0, 20.0),
            1.0, atol = TestNumerics.absolute_floor(Float64))
        # Negative alpha (unusual, but mathematically valid)
        @test isapprox(calc_temperature_correction(-0.001, 30.0, 20.0),
            0.99, atol = TestNumerics.absolute_floor(Float64))
        # Large temperature difference within ΔTmax
        @test isapprox(calc_temperature_correction(0.00393, 20.0 + (ΔTmax - 1), 20.0),
            1 + 0.00393 * (ΔTmax - 1), atol = TestNumerics.absolute_floor(Float64))
    end

    # Numerical Consistency
    @testset "Numerical Consistency" begin
        # Float32
        kf = calc_temperature_correction(Float32(0.00393), Float32(75.0), Float32(20.0))
        @test isapprox(kf, 1.2161f0, atol = Float32(1e-4))
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
        @test isapprox(value(km), μ, atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(uncertainty(km), sqrt(σ2), atol = TestNumerics.absolute_floor(Float64))
    end
end
@testitem "BaseParams / calc_tubular_resistance / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    # Basic Functionality
    @testset "Basic Functionality" begin
        # Example from docstring
        r_in = 0.01
        r_ex = 0.02
        rho = 1.7241e-8
        alpha = 0.00393
        T0 = 20.0
        Top = 25.0
        expected = calc_temperature_correction(alpha, Top, T0) * rho /
                   (π * (r_ex^2 - r_in^2))
        R = calc_tubular_resistance(r_in, r_ex, rho, alpha, T0, Top)
        @test isapprox(R, expected, atol = TestNumerics.absolute_floor(Float64))
    end

    # Edge Cases
    @testset "Edge Cases" begin
        # Zero thickness (r_in == r_ex): cross-section = 0, expect Inf or error
        r_in = 0.01
        r_ext = 0.01
        rho = 1.7241e-8
        alpha = 0.00393
        T0 = 20.0
        Top = 25.0
        # Should return Inf (division by zero)
        R = calc_tubular_resistance(r_in, r_ext, rho, alpha, T0, Top)
        @test isinf(R)
        # Very thin tube (r_ex - r_in ≈ eps)
        r_in2 = 0.01
        r_ext2 = 0.01 + eps()
        R2 = calc_tubular_resistance(r_in2, r_ext2, rho, alpha, T0, Top)
        @test R2 > 0
        # Large radii
        R3 = calc_tubular_resistance(1.0, 2.0, rho, alpha, T0, Top)
        @test R3 < 1e-8
        # Negative temperature coefficient (mathematically valid)
        R4 = calc_tubular_resistance(0.01, 0.02, rho, -0.001, T0, Top)
        expected4 = calc_temperature_correction(-0.001, Top, T0) * rho /
                    (π * (0.02^2 - 0.01^2))
        @test isapprox(R4, expected4, atol = TestNumerics.absolute_floor(Float64))
    end

    # Numerical Consistency
    @testset "Numerical Consistency" begin
        # Float32
        Rf = calc_tubular_resistance(
            Float32(0.01),
            Float32(0.02),
            Float32(1.7241e-8),
            Float32(0.00393),
            Float32(20.0),
            Float32(25.0)
        )
        expectedf = calc_temperature_correction(Float32(0.00393), Float32(25.0), Float32(20.0)) *
                    Float32(1.7241e-8) / (π * (Float32(0.02)^2 - Float32(0.01)^2))
        @test isapprox(Rf, expectedf, atol = TestNumerics.absolute_floor(Float32))
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
        @test isapprox(value(Rm), μ, atol = TestNumerics.absolute_floor(Float64))
        # Uncertainty should be nonzero and scale with input uncertainties
        @test uncertainty(Rm) > 0
    end
end
@testitem "BaseParams / calc_strip_resistance / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    using Measurements
    # --- Basic Functionality ---
    @testset "Basic Functionality" begin
        thickness = 0.002
        width = 0.05
        rho = 1.7241e-8
        alpha = 0.00393
        T0 = 20.0
        Top = 25.0
        R = calc_strip_resistance(thickness, width, rho, alpha, T0, Top)
        @test isapprox(R, 0.00017579785649999996, atol = TestNumerics.absolute_floor(Float64))
    end

    # --- Edge Cases ---
    @testset "Edge Cases" begin
        # Zero thickness (should return Inf or error in physical context, but function will return Inf)
        R = calc_strip_resistance(0.0, 0.05, 1.7241e-8, 0.00393, 20.0, 25.0)
        @test isinf(R)
        # Zero width
        R = calc_strip_resistance(0.002, 0.0, 1.7241e-8, 0.00393, 20.0, 25.0)
        @test isinf(R)
        # Large temperature difference (within asserted range)
        R = calc_strip_resistance(0.002, 0.05, 1.7241e-8, 0.00393, 20.0, 100.0)
        @test R > 0.00017579785649999996
    end

    # --- Numerical Consistency ---
    @testset "Numerical Consistency" begin
        # Float32
        R = calc_strip_resistance(Float32(0.002), Float32(0.05), Float32(1.7241e-8),
            Float32(0.00393), Float32(20.0), Float32(25.0))
        @test isapprox(R, 0.00017579785649999996, atol = TestNumerics.absolute_floor(Float64))
    end

    # --- Physical Behavior ---
    @testset "Physical Behavior" begin
        # Resistance increases with temperature
        R1 = calc_strip_resistance(0.002, 0.05, 1.7241e-8, 0.00393, 20.0, 25.0)
        R2 = calc_strip_resistance(0.002, 0.05, 1.7241e-8, 0.00393, 20.0, 75.0)
        @test R2 > R1
        # Resistance decreases with increasing cross-section
        R3 = calc_strip_resistance(0.002, 0.05, 1.7241e-8, 0.00393, 20.0, 25.0)
        R4 = calc_strip_resistance(0.004, 0.05, 1.7241e-8, 0.00393, 20.0, 25.0)
        @test R4 < R3
    end

    # --- Type Stability & Promotion ---
    @testset "Type Stability & Promotion" begin
        # All Float64
        R = calc_strip_resistance(0.002, 0.05, 1.7241e-8, 0.00393, 20.0, 25.0)
        @test typeof(R) == Float64
        # All Measurement
        Rm = calc_strip_resistance(measurement(0.002, 1e-6), measurement(0.05, 1e-5),
            measurement(1.7241e-8, 1e-10), measurement(0.00393, 1e-6),
            measurement(20.0, 0.1), measurement(25.0, 0.1))
        @test Rm isa Measurement{Float64}
        # Mixed: thickness as Measurement
        R1 = calc_strip_resistance(
            measurement(0.002, 1e-6), 0.05, 1.7241e-8, 0.00393, 20.0, 25.0)
        @test R1 isa Measurement{Float64}
        # Mixed: alpha as Measurement
        R2 = calc_strip_resistance(
            0.002, 0.05, 1.7241e-8, measurement(0.00393, 1e-6), 20.0, 25.0)
        @test R2 isa Measurement{Float64}
    end

    # --- Uncertainty Quantification ---
    @testset "Uncertainty Quantification" begin
        t = measurement(0.002, 1e-6)
        w = measurement(0.05, 1e-5)
        r = measurement(1.7241e-8, 1e-10)
        a = measurement(0.00393, 1e-6)
        t0 = measurement(20.0, 0.1)
        top = measurement(25.0, 0.1)
        R = calc_strip_resistance(t, w, r, a, t0, top)
        @test uncertainty(R) > 0
    end
end
@testitem "BaseParams / calc_parallel_equivalent / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        # Test with real numbers (Float64)
        Z1_real = 5.0
        Z2_real = 10.0
        expected_real = 1 / (1 / Z1_real + 1 / Z2_real)
        result_real = calc_parallel_equivalent(Z1_real, Z2_real)
        @test isapprox(result_real, expected_real; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(result_real, 3.3333333333333335; atol = TestNumerics.absolute_floor(Float64))

        # Test with complex numbers (Complex{Float64})
        Z1_complex = 3.0 + 4.0im
        Z2_complex = 8.0 - 6.0im
        expected_complex = 1 / (1 / Z1_complex + 1 / Z2_complex)
        @test isapprox(calc_parallel_equivalent(Z1_complex, Z2_complex),
            expected_complex; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "Edge Cases" begin
        # Zero impedance (short circuit)
        @test isapprox(calc_parallel_equivalent(0.0, 10.0), 0.0; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(calc_parallel_equivalent(10.0, 0.0), 0.0; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(calc_parallel_equivalent(0.0, 0.0), 0.0; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(calc_parallel_equivalent(0.0 + 0.0im, 5.0 + 5.0im),
            0.0 + 0.0im; atol = TestNumerics.absolute_floor(Float64))

        # Infinite impedance (open circuit)
        @test isapprox(calc_parallel_equivalent(Inf, 10.0), 10.0; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(calc_parallel_equivalent(10.0, Inf), 10.0; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(calc_parallel_equivalent(Inf, Inf), Inf; atol = TestNumerics.absolute_floor(Float64))

        # NaN propagation
        @test isnan(calc_parallel_equivalent(NaN, 10.0))
        @test isnan(calc_parallel_equivalent(10.0, NaN))

        # Equal and opposite impedances (Z1 = -Z2), leading to singularity
        result_inf = calc_parallel_equivalent(10.0, -10.0)
        @test isinf(real(result_inf))
        result_nan = calc_parallel_equivalent(3.0 + 4.0im, -3.0 - 4.0im)
        @test isnan(real(result_nan)) && isnan(imag(result_nan))
    end

    @testset "Numerical Consistency" begin
        Z1f = 5.0
        Z2f = 10.0
        resultf = calc_parallel_equivalent(Z1f, Z2f)
        @test resultf isa Float64
        @test isapprox(resultf, 3.33333333; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "Physical Behavior" begin
        # Parallel resistance is always less than the smallest individual resistance
        @test calc_parallel_equivalent(10.0, 20.0) < 10.0

        # Symmetry: calc_parallel_equivalent(Z1, Z2) == calc_parallel_equivalent(Z2, Z1)
        @test isapprox(calc_parallel_equivalent(7.0, 13.0),
            calc_parallel_equivalent(13.0, 7.0); atol = TestNumerics.absolute_floor(Float64))

        # If Z1 == Z2, the result is Z1 / 2
        @test isapprox(calc_parallel_equivalent(8.0, 8.0), 4.0; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "Type Stability & Promotion" begin
        # Both Float64 -> Float64
        @test calc_parallel_equivalent(5.0, 10.0) isa Float64

        # Int and Float64 -> Float64
        result_mixed_real = calc_parallel_equivalent(5, 10.0)
        @test result_mixed_real isa Float64
        @test isapprox(result_mixed_real, 1 / (1 / 5.0 + 1 / 10.0); atol = TestNumerics.absolute_floor(Float64))

        # Float64 and Complex{Float64} -> Complex{Float64}
        result_mixed_complex = calc_parallel_equivalent(10.0, 3.0 + 4.0im)
        @test result_mixed_complex isa Complex{Float64}
        expected_mixed_complex = 1 / (1 / (10.0 + 0.0im) + 1 / (3.0 + 4.0im))
        @test isapprox(result_mixed_complex, expected_mixed_complex;
            atol = TestNumerics.absolute_floor(Float64))

        # Both Measurement -> Measurement
        Z1m = measurement(5.0, 0.1)
        Z2m = measurement(10.0, 0.2)
        @test calc_parallel_equivalent(Z1m, Z2m) isa Measurement

        # Mixed: Measurement and Float64 -> Measurement
        @test calc_parallel_equivalent(Z1m, 10.0) isa Measurement
        @test calc_parallel_equivalent(5.0, Z2m) isa Measurement
    end

    @testset "Uncertainty Quantification with Measurements.jl" begin
        # Mixed Case 1: First argument is a Measurement
        Z1_meas = measurement(5.0, 0.1)
        Z2_float = 10.0
        result_mixed1 = calc_parallel_equivalent(Z1_meas, Z2_float)
        expected_mixed1 = 1 / (1 / Z1_meas + 1 / Z2_float)
        @test result_mixed1 isa Measurement{Float64}
        @test isapprox(value(result_mixed1), value(expected_mixed1);
            atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(uncertainty(result_mixed1), uncertainty(expected_mixed1);
            atol = TestNumerics.absolute_floor(Float64))

        # Mixed Case 2: Second argument is a Measurement
        Z1_float = 5.0
        Z2_meas = measurement(10.0, 0.2)
        result_mixed2 = calc_parallel_equivalent(Z1_float, Z2_meas)
        expected_mixed2 = 1 / (1 / Z1_float + 1 / Z2_meas)
        @test result_mixed2 isa Measurement{Float64}
        @test isapprox(value(result_mixed2), value(expected_mixed2);
            atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(uncertainty(result_mixed2), uncertainty(expected_mixed2);
            atol = TestNumerics.absolute_floor(Float64))

        # Fully Promoted Case: Both inputs are Measurements
        result_full_meas = calc_parallel_equivalent(Z1_meas, Z2_meas)
        expected_full_meas = 1 / (1 / Z1_meas + 1 / Z2_meas)
        @test result_full_meas isa Measurement{Float64}
        @test isapprox(value(result_full_meas), value(expected_full_meas);
            atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(uncertainty(result_full_meas), uncertainty(expected_full_meas);
            atol = TestNumerics.absolute_floor(Float64))

        # Fully Promoted Complex Case
        Z1_cplx_meas = measurement(3.0, 0.1) + measurement(4.0, 0.2)im
        Z2_cplx_meas = measurement(8.0, 0.3) - measurement(6.0, 0.4)im
        result_cplx_meas = calc_parallel_equivalent(Z1_cplx_meas, Z2_cplx_meas)
        expected_cplx_meas = 1 / (1 / Z1_cplx_meas + 1 / Z2_cplx_meas)
        @test result_cplx_meas isa Complex{Measurement{Float64}}
        @test isapprox(value(real(result_cplx_meas)), value(real(expected_cplx_meas));
            atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(value(imag(result_cplx_meas)), value(imag(expected_cplx_meas));
            atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(uncertainty(real(result_cplx_meas)),
            uncertainty(real(expected_cplx_meas)); atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(uncertainty(imag(result_cplx_meas)),
            uncertainty(imag(expected_cplx_meas)); atol = TestNumerics.absolute_floor(Float64))
    end
end
@testitem "BaseParams / calc_equivalent_alpha / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "calc_equivalent_alpha: Basic Functionality (Copper & Aluminum)" begin
        alpha1 = 0.00393  # Copper
        R1 = 0.5
        alpha2 = 0.00403  # Aluminum
        R2 = 1.0
        expected = (alpha1 * R2 + alpha2 * R1) / (R1 + R2)
        result = calc_equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "calc_equivalent_alpha: Edge Case - Zero Resistance" begin
        alpha1 = 0.00393
        R1 = 0.0
        alpha2 = 0.00403
        R2 = 1.0
        expected = alpha1  # Only R2 matters
        result = calc_equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))

        alpha1 = 0.00393
        R1 = 0.5
        alpha2 = 0.00403
        R2 = 0.0
        expected = alpha2  # Only R1 matters
        result = calc_equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "calc_equivalent_alpha: Edge Case - Very Large Resistance" begin
        alpha1 = 0.00393
        R1 = 1e12
        alpha2 = 0.00403
        R2 = 1.0
        expected = (alpha1 * R2 + alpha2 * R1) / (R1 + R2)
        result = calc_equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "calc_equivalent_alpha: Type Stability & Promotion" begin
        alpha1 = 0.00393
        R1 = 0.5
        alpha2 = 0.00403
        R2 = 1
        result = calc_equivalent_alpha(alpha1, R1, alpha2, R2)
        @test typeof(result) == typeof(1.0)

        # Mixed Int/Float
        result2 = calc_equivalent_alpha(0, 1, 1, 1.0)
        @test typeof(result2) == typeof(1.0)
    end

    @testset "calc_equivalent_alpha: Uncertainty Quantification (Measurements.jl)" begin
        alpha1 = measurement(0.00393, 1e-5)
        R1 = measurement(0.5, 1e-3)
        alpha2 = measurement(0.00403, 1e-5)
        R2 = measurement(1.0, 1e-3)
        result = calc_equivalent_alpha(alpha1, R1, alpha2, R2)
        # Check value
        expected_val = (value(alpha1) * value(R2) + value(alpha2) * value(R1)) /
                       (value(R1) + value(R2))
        @test isapprox(value(result), expected_val; atol = TestNumerics.absolute_floor(Float64))
        # Check uncertainty propagation (should be nonzero)
        @test uncertainty(result) > 0
    end

    @testset "calc_equivalent_alpha: Equivalent temperature coefficient for parallel resistors" begin

        # Example values (Copper and Aluminum)
        alpha1 = 0.00393  # Copper
        R1 = 0.5
        alpha2 = 0.00403 # Aluminum
        R2 = 1.0

        # Analytical result
        expected = (alpha1 * R2 + alpha2 * R1) / (R1 + R2)
        result = calc_equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))

        # Edge case: Identical conductors
        alpha = 0.00393
        R = 1.0
        @test isapprox(calc_equivalent_alpha(alpha, R, alpha, R), alpha;
            atol = TestNumerics.absolute_floor(Float64))

        # Edge case: One resistance much larger than the other
        @test isapprox(calc_equivalent_alpha(0.003, 1e6, 0.005, 1.0),
            0.005; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(calc_equivalent_alpha(0.003, 1.0, 0.005, 1e6),
            0.003; atol = TestNumerics.absolute_floor(Float64))

        # Type promotion and Measurements.jl propagation
        using Measurements: ±, value, uncertainty
        m1 = 0.00393 ± 0.00001
        m2 = 0.00403 ± 0.00001
        r1 = 0.5 ± 0.01
        r2 = 1.0 ± 0.01

        @testset "Type Promotion with Measurements.jl" begin
            # Base case: All Float64
            @test calc_equivalent_alpha(alpha1, R1, alpha2, R2) isa Float64

            # Fully promoted: All Measurement
            res = calc_equivalent_alpha(m1, r1, m2, r2)
            @test res isa Measurement{Float64}
            @test isapprox(value(res), expected; atol = TestNumerics.absolute_floor(Float64))
            # Uncertainty should be nonzero
            @test uncertainty(res) > 0

            # Mixed case 1: First argument is Measurement
            res = calc_equivalent_alpha(m1, R1, alpha2, R2)
            @test res isa Measurement{Float64}
            @test isapprox(value(res), expected; atol = TestNumerics.absolute_floor(Float64))

            # Mixed case 2: Middle argument is Measurement
            res = calc_equivalent_alpha(alpha1, r1, alpha2, R2)
            @test res isa Measurement{Float64}
            @test isapprox(value(res), expected; atol = TestNumerics.absolute_floor(Float64))

            # Mixed case 3: Last argument is Measurement
            res = calc_equivalent_alpha(alpha1, R1, alpha2, r2)
            @test res isa Measurement{Float64}
            @test isapprox(value(res), expected; atol = TestNumerics.absolute_floor(Float64))
        end

        # Physically unusual but valid: zero resistance (should return NaN)
        @test isnan(calc_equivalent_alpha(0.003, 0.0, 0.005, 0.0))

        # Large values
        @test isapprox(calc_equivalent_alpha(1e-3, 1e6, 2e-3, 2e6),
            (1e-3 * 2e6 + 2e-3 * 1e6) / (1e6 + 2e6); atol = TestNumerics.absolute_floor(Float64))
    end
end # End of test file
@testitem "BaseParams / calc_equivalent_rho / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        # Example from docstring: R=0.01 Ω, r_ext=0.02 m, r_in=0.01 m
        result = calc_equivalent_rho(0.01, 0.02, 0.01)
        expected = 0.01 * π * (0.02^2 - 0.01^2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
        @test result > 0
    end

    @testset "Edge Cases" begin
        # Zero resistance
        result = calc_equivalent_rho(0.0, 0.02, 0.01)
        @test isapprox(result, 0.0; atol = TestNumerics.absolute_floor(Float64))
        # Zero thickness (r_ext == r_in)
        result = calc_equivalent_rho(0.01, 0.01, 0.01)
        @test isapprox(result, 0.0; atol = TestNumerics.absolute_floor(Float64))
        # Very large radii
        result = calc_equivalent_rho(0.01, 1e6, 1e3)
        expected = 0.01 * π * (1e6^2 - 1e3^2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
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
        @test isapprox(r, d; atol = TestNumerics.absolute_floor(Float64))
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
