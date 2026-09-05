@testitem "BaseParams / tubular_resistance / reference-state contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    using Measurements

    r_in=0.01
    r_ex=0.02
    rho=1.7241e-8
    expected=rho/(π*(r_ex^2-r_in^2))
    @test isapprox(
        tubular_resistance(r_in, r_ex, rho),
        expected;
        atol = TestNumerics.absolute_floor(Float64)
    )
    @test isinf(tubular_resistance(r_in, r_in, rho))
    @test tubular_resistance(r_in, r_in + eps(), rho) > 0
    @test tubular_resistance(r_in, 0.03, rho) <
          tubular_resistance(r_in, 0.015, rho)
    @test tubular_resistance(r_in, r_ex, 1e-7) >
          tubular_resistance(r_in, r_ex, 1e-8)

    result32=tubular_resistance(
        Float32(r_in),
        Float32(r_ex),
        Float32(rho)
    )
    @test result32 isa Float32

    uncertain=tubular_resistance(
        measurement(r_in, 1e-6),
        measurement(r_ex, 1e-6),
        measurement(rho, 1e-10)
    )
    @test uncertain isa Measurement{Float64}
    @test uncertainty(uncertain) > 0

    @test_throws MethodError tubular_resistance(
        r_in,
        r_ex,
        rho,
        0.00393,
        20.0,
        75.0
    )
end

@testitem "BaseParams / strip_resistance / reference-state contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    using Measurements

    thickness=0.002
    width=0.05
    rho=1.7241e-8
    expected=rho/(thickness*width)
    @test isapprox(
        strip_resistance(thickness, width, rho),
        expected;
        atol = TestNumerics.absolute_floor(Float64)
    )
    @test isinf(strip_resistance(0.0, width, rho))
    @test isinf(strip_resistance(thickness, 0.0, rho))
    @test strip_resistance(2thickness, width, rho) <
          strip_resistance(thickness, width, rho)
    @test strip_resistance(thickness, width, 1e-7) >
          strip_resistance(thickness, width, 1e-8)

    result32=strip_resistance(
        Float32(thickness),
        Float32(width),
        Float32(rho)
    )
    @test result32 isa Float32

    uncertain=strip_resistance(
        measurement(thickness, 1e-6),
        measurement(width, 1e-5),
        measurement(rho, 1e-10)
    )
    @test uncertain isa Measurement{Float64}
    @test uncertainty(uncertain) > 0

    @test_throws MethodError strip_resistance(
        thickness,
        width,
        rho,
        0.00393,
        20.0,
        75.0
    )
end

@testitem "BaseParams / parallel / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        # Test with real numbers (Float64)
        Z1_real = 5.0
        Z2_real = 10.0
        expected_real = 1 / (1 / Z1_real + 1 / Z2_real)
        result_real = parallel(Z1_real, Z2_real)
        @test isapprox(result_real, expected_real; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(result_real, 3.3333333333333335; atol = TestNumerics.absolute_floor(Float64))

        # Test with complex numbers (Complex{Float64})
        Z1_complex = 3.0 + 4.0im
        Z2_complex = 8.0 - 6.0im
        expected_complex = 1 / (1 / Z1_complex + 1 / Z2_complex)
        @test isapprox(parallel(Z1_complex, Z2_complex),
            expected_complex; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "Edge Cases" begin
        # Zero impedance (short circuit)
        @test isapprox(parallel(0.0, 10.0), 0.0; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(parallel(10.0, 0.0), 0.0; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(parallel(0.0, 0.0), 0.0; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(parallel(0.0 + 0.0im, 5.0 + 5.0im),
            0.0 + 0.0im; atol = TestNumerics.absolute_floor(Float64))

        # Infinite impedance (open circuit)
        @test isapprox(parallel(Inf, 10.0), 10.0; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(parallel(10.0, Inf), 10.0; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(parallel(Inf, Inf), Inf; atol = TestNumerics.absolute_floor(Float64))

        # NaN propagation
        @test isnan(parallel(NaN, 10.0))
        @test isnan(parallel(10.0, NaN))

        # Equal and opposite impedances (Z1 = -Z2), leading to singularity
        result_inf = parallel(10.0, -10.0)
        @test isinf(real(result_inf))
        result_nan = parallel(3.0 + 4.0im, -3.0 - 4.0im)
        @test isnan(real(result_nan)) && isnan(imag(result_nan))
    end

    @testset "Numerical Consistency" begin
        Z1f = 5.0
        Z2f = 10.0
        resultf = parallel(Z1f, Z2f)
        @test resultf isa Float64
        @test isapprox(resultf, 3.33333333; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "Physical Behavior" begin
        # Parallel resistance is always less than the smallest individual resistance
        @test parallel(10.0, 20.0) < 10.0

        # Symmetry: parallel(Z1, Z2) == parallel(Z2, Z1)
        @test isapprox(parallel(7.0, 13.0),
            parallel(13.0, 7.0); atol = TestNumerics.absolute_floor(Float64))

        # If Z1 == Z2, the result is Z1 / 2
        @test isapprox(parallel(8.0, 8.0), 4.0; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "Type Stability & Promotion" begin
        # Both Float64 -> Float64
        @test parallel(5.0, 10.0) isa Float64

        # Int and Float64 -> Float64
        result_mixed_real = parallel(5, 10.0)
        @test result_mixed_real isa Float64
        @test isapprox(result_mixed_real, 1 / (1 / 5.0 + 1 / 10.0); atol = TestNumerics.absolute_floor(Float64))

        # Float64 and Complex{Float64} -> Complex{Float64}
        result_mixed_complex = parallel(10.0, 3.0 + 4.0im)
        @test result_mixed_complex isa Complex{Float64}
        expected_mixed_complex = 1 / (1 / (10.0 + 0.0im) + 1 / (3.0 + 4.0im))
        @test isapprox(result_mixed_complex, expected_mixed_complex;
            atol = TestNumerics.absolute_floor(Float64))

        # Both Measurement -> Measurement
        Z1m = measurement(5.0, 0.1)
        Z2m = measurement(10.0, 0.2)
        @test parallel(Z1m, Z2m) isa Measurement

        # Mixed: Measurement and Float64 -> Measurement
        @test parallel(Z1m, 10.0) isa Measurement
        @test parallel(5.0, Z2m) isa Measurement
    end

    @testset "Uncertainty Quantification with Measurements.jl" begin
        # Mixed Case 1: First argument is a Measurement
        Z1_meas = measurement(5.0, 0.1)
        Z2_float = 10.0
        result_mixed1 = parallel(Z1_meas, Z2_float)
        expected_mixed1 = 1 / (1 / Z1_meas + 1 / Z2_float)
        @test result_mixed1 isa Measurement{Float64}
        @test isapprox(value(result_mixed1), value(expected_mixed1);
            atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(uncertainty(result_mixed1), uncertainty(expected_mixed1);
            atol = TestNumerics.absolute_floor(Float64))

        # Mixed Case 2: Second argument is a Measurement
        Z1_float = 5.0
        Z2_meas = measurement(10.0, 0.2)
        result_mixed2 = parallel(Z1_float, Z2_meas)
        expected_mixed2 = 1 / (1 / Z1_float + 1 / Z2_meas)
        @test result_mixed2 isa Measurement{Float64}
        @test isapprox(value(result_mixed2), value(expected_mixed2);
            atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(uncertainty(result_mixed2), uncertainty(expected_mixed2);
            atol = TestNumerics.absolute_floor(Float64))

        # Fully Promoted Case: Both inputs are Measurements
        result_full_meas = parallel(Z1_meas, Z2_meas)
        expected_full_meas = 1 / (1 / Z1_meas + 1 / Z2_meas)
        @test result_full_meas isa Measurement{Float64}
        @test isapprox(value(result_full_meas), value(expected_full_meas);
            atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(uncertainty(result_full_meas), uncertainty(expected_full_meas);
            atol = TestNumerics.absolute_floor(Float64))

        # Fully Promoted Complex Case
        Z1_cplx_meas = measurement(3.0, 0.1) + measurement(4.0, 0.2)im
        Z2_cplx_meas = measurement(8.0, 0.3) - measurement(6.0, 0.4)im
        result_cplx_meas = parallel(Z1_cplx_meas, Z2_cplx_meas)
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
@testitem "BaseParams / equivalent_alpha / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "equivalent_alpha: Basic Functionality (Copper & Aluminum)" begin
        alpha1 = 0.00393  # Copper
        R1 = 0.5
        alpha2 = 0.00403  # Aluminum
        R2 = 1.0
        expected = (alpha1 * R2 + alpha2 * R1) / (R1 + R2)
        result = equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "equivalent_alpha: Edge Case - Zero Resistance" begin
        alpha1 = 0.00393
        R1 = 0.0
        alpha2 = 0.00403
        R2 = 1.0
        expected = alpha1  # Only R2 matters
        result = equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))

        alpha1 = 0.00393
        R1 = 0.5
        alpha2 = 0.00403
        R2 = 0.0
        expected = alpha2  # Only R1 matters
        result = equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "equivalent_alpha: Edge Case - Very Large Resistance" begin
        alpha1 = 0.00393
        R1 = 1e12
        alpha2 = 0.00403
        R2 = 1.0
        expected = (alpha1 * R2 + alpha2 * R1) / (R1 + R2)
        result = equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "equivalent_alpha: Type Stability & Promotion" begin
        alpha1 = 0.00393
        R1 = 0.5
        alpha2 = 0.00403
        R2 = 1
        result = equivalent_alpha(alpha1, R1, alpha2, R2)
        @test typeof(result) == typeof(1.0)

        # Mixed Int/Float
        result2 = equivalent_alpha(0, 1, 1, 1.0)
        @test typeof(result2) == typeof(1.0)
    end

    @testset "equivalent_alpha: Uncertainty Quantification (Measurements.jl)" begin
        alpha1 = measurement(0.00393, 1e-5)
        R1 = measurement(0.5, 1e-3)
        alpha2 = measurement(0.00403, 1e-5)
        R2 = measurement(1.0, 1e-3)
        result = equivalent_alpha(alpha1, R1, alpha2, R2)
        # Check value
        expected_val = (value(alpha1) * value(R2) + value(alpha2) * value(R1)) /
                       (value(R1) + value(R2))
        @test isapprox(value(result), expected_val; atol = TestNumerics.absolute_floor(Float64))
        # Check uncertainty propagation (should be nonzero)
        @test uncertainty(result) > 0
    end

    @testset "equivalent_alpha: Equivalent temperature coefficient for parallel resistors" begin

        # Example values (Copper and Aluminum)
        alpha1 = 0.00393  # Copper
        R1 = 0.5
        alpha2 = 0.00403 # Aluminum
        R2 = 1.0

        # Analytical result
        expected = (alpha1 * R2 + alpha2 * R1) / (R1 + R2)
        result = equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))

        # Edge case: Identical conductors
        alpha = 0.00393
        R = 1.0
        @test isapprox(equivalent_alpha(alpha, R, alpha, R), alpha;
            atol = TestNumerics.absolute_floor(Float64))

        # Edge case: One resistance much larger than the other
        @test isapprox(equivalent_alpha(0.003, 1e6, 0.005, 1.0),
            0.005; atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(equivalent_alpha(0.003, 1.0, 0.005, 1e6),
            0.003; atol = TestNumerics.absolute_floor(Float64))

        # Type promotion and Measurements.jl propagation
        using Measurements: ±, value, uncertainty
        m1 = 0.00393 ± 0.00001
        m2 = 0.00403 ± 0.00001
        r1 = 0.5 ± 0.01
        r2 = 1.0 ± 0.01

        @testset "Type Promotion with Measurements.jl" begin
            # Base case: All Float64
            @test equivalent_alpha(alpha1, R1, alpha2, R2) isa Float64

            # Fully promoted: All Measurement
            res = equivalent_alpha(m1, r1, m2, r2)
            @test res isa Measurement{Float64}
            @test isapprox(value(res), expected; atol = TestNumerics.absolute_floor(Float64))
            # Uncertainty should be nonzero
            @test uncertainty(res) > 0

            # Mixed case 1: First argument is Measurement
            res = equivalent_alpha(m1, R1, alpha2, R2)
            @test res isa Measurement{Float64}
            @test isapprox(value(res), expected; atol = TestNumerics.absolute_floor(Float64))

            # Mixed case 2: Middle argument is Measurement
            res = equivalent_alpha(alpha1, r1, alpha2, R2)
            @test res isa Measurement{Float64}
            @test isapprox(value(res), expected; atol = TestNumerics.absolute_floor(Float64))

            # Mixed case 3: Last argument is Measurement
            res = equivalent_alpha(alpha1, R1, alpha2, r2)
            @test res isa Measurement{Float64}
            @test isapprox(value(res), expected; atol = TestNumerics.absolute_floor(Float64))
        end

        # Physically unusual but valid: zero resistance (should return NaN)
        @test isnan(equivalent_alpha(0.003, 0.0, 0.005, 0.0))

        # Large values
        @test isapprox(equivalent_alpha(1e-3, 1e6, 2e-3, 2e6),
            (1e-3 * 2e6 + 2e-3 * 1e6) / (1e6 + 2e6); atol = TestNumerics.absolute_floor(Float64))
    end
end # End of test file
@testitem "BaseParams / equivalent_rho / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        # Example from docstring: R=0.01 Ω, r_ext=0.02 m, r_in=0.01 m
        result = equivalent_rho(0.01, 0.02, 0.01)
        expected = 0.01 * π * (0.02^2 - 0.01^2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
        @test result > 0
    end

    @testset "Edge Cases" begin
        # Zero resistance
        result = equivalent_rho(0.0, 0.02, 0.01)
        @test isapprox(result, 0.0; atol = TestNumerics.absolute_floor(Float64))
        # Zero thickness (r_ext == r_in)
        result = equivalent_rho(0.01, 0.01, 0.01)
        @test isapprox(result, 0.0; atol = TestNumerics.absolute_floor(Float64))
        # Very large radii
        result = equivalent_rho(0.01, 1e6, 1e3)
        expected = 0.01 * π * (1e6^2 - 1e3^2)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
        # Inf/NaN
        @test isnan(equivalent_rho(NaN, 0.02, 0.01))
        @test isnan(equivalent_rho(0.01, NaN, 0.01))
        @test isnan(equivalent_rho(0.01, 0.02, NaN))
        @test isinf(equivalent_rho(Inf, 0.02, 0.01))
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        r = equivalent_rho(Float32(0.01), Float32(0.02), Float32(0.01))
        d = equivalent_rho(0.01, 0.02, 0.01)
        @test isapprox(r, d; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "Physical Behavior" begin
        # Increases with resistance
        r1 = equivalent_rho(0.01, 0.02, 0.01)
        r2 = equivalent_rho(0.02, 0.02, 0.01)
        @test r2 > r1
        # Increases with area
        r3 = equivalent_rho(0.01, 0.03, 0.01)
        @test r3 > r1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        # All Float64
        r1 = equivalent_rho(0.01, 0.02, 0.01)
        @test typeof(r1) == Float64
        # All Measurement
        r2 = equivalent_rho(measurement(0.01, 1e-4), measurement(0.02, 1e-5), measurement(0.01, 1e-5))
        @test r2 isa Measurement{Float64}
        # Mixed: R as Measurement
        r3 = equivalent_rho(measurement(0.01, 1e-4), 0.02, 0.01)
        @test r3 isa Measurement{Float64}
        # Mixed: radius_ext_con as Measurement
        r4 = equivalent_rho(0.01, measurement(0.02, 1e-5), 0.01)
        @test r4 isa Measurement{Float64}
        # Mixed: radius_in_con as Measurement
        r5 = equivalent_rho(0.01, 0.02, measurement(0.01, 1e-5))
        @test r5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        R = measurement(0.01, 1e-4)
        r_ext = measurement(0.02, 1e-5)
        r_in = measurement(0.01, 1e-5)
        result = equivalent_rho(R, r_ext, r_in)
        @test result isa Measurement{Float64}
        @test uncertainty(result) > 0
    end
end
