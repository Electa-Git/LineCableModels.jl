@testitem "BaseParams: calc_parallel_equivalent unit tests" setup = [defaults] begin

    @testset "Basic Functionality" begin
        # Test with real numbers (Float64)
        Z1_real = 5.0
        Z2_real = 10.0
        expected_real = 1 / (1 / Z1_real + 1 / Z2_real)
        result_real = calc_parallel_equivalent(Z1_real, Z2_real)
        @test isapprox(result_real, expected_real; atol=TEST_TOL)
        @test isapprox(result_real, 3.3333333333333335; atol=TEST_TOL)

        # Test with complex numbers (Complex{Float64})
        Z1_complex = 3.0 + 4.0im
        Z2_complex = 8.0 - 6.0im
        expected_complex = 1 / (1 / Z1_complex + 1 / Z2_complex)
        @test isapprox(calc_parallel_equivalent(Z1_complex, Z2_complex), expected_complex; atol=TEST_TOL)
    end

    @testset "Edge Cases" begin
        # Zero impedance (short circuit)
        @test isapprox(calc_parallel_equivalent(0.0, 10.0), 0.0; atol=TEST_TOL)
        @test isapprox(calc_parallel_equivalent(10.0, 0.0), 0.0; atol=TEST_TOL)
        @test isapprox(calc_parallel_equivalent(0.0, 0.0), 0.0; atol=TEST_TOL)
        @test isapprox(calc_parallel_equivalent(0.0 + 0.0im, 5.0 + 5.0im), 0.0 + 0.0im; atol=TEST_TOL)

        # Infinite impedance (open circuit)
        @test isapprox(calc_parallel_equivalent(Inf, 10.0), 10.0; atol=TEST_TOL)
        @test isapprox(calc_parallel_equivalent(10.0, Inf), 10.0; atol=TEST_TOL)
        @test isapprox(calc_parallel_equivalent(Inf, Inf), Inf; atol=TEST_TOL)

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
        @test isapprox(resultf, 3.33333333; atol=TEST_TOL)
    end

    @testset "Physical Behavior" begin
        # Parallel resistance is always less than the smallest individual resistance
        @test calc_parallel_equivalent(10.0, 20.0) < 10.0

        # Symmetry: calc_parallel_equivalent(Z1, Z2) == calc_parallel_equivalent(Z2, Z1)
        @test isapprox(calc_parallel_equivalent(7.0, 13.0), calc_parallel_equivalent(13.0, 7.0); atol=TEST_TOL)

        # If Z1 == Z2, the result is Z1 / 2
        @test isapprox(calc_parallel_equivalent(8.0, 8.0), 4.0; atol=TEST_TOL)
    end

    @testset "Type Stability & Promotion" begin
        # Both Float64 -> Float64
        @test calc_parallel_equivalent(5.0, 10.0) isa Float64

        # Int and Float64 -> Float64
        result_mixed_real = calc_parallel_equivalent(5, 10.0)
        @test result_mixed_real isa Float64
        @test isapprox(result_mixed_real, 1 / (1 / 5.0 + 1 / 10.0); atol=TEST_TOL)

        # Float64 and Complex{Float64} -> Complex{Float64}
        result_mixed_complex = calc_parallel_equivalent(10.0, 3.0 + 4.0im)
        @test result_mixed_complex isa Complex{Float64}
        expected_mixed_complex = 1 / (1 / (10.0 + 0.0im) + 1 / (3.0 + 4.0im))
        @test isapprox(result_mixed_complex, expected_mixed_complex; atol=TEST_TOL)

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
        @test isapprox(value(result_mixed1), value(expected_mixed1); atol=TEST_TOL)
        @test isapprox(uncertainty(result_mixed1), uncertainty(expected_mixed1); atol=TEST_TOL)

        # Mixed Case 2: Second argument is a Measurement
        Z1_float = 5.0
        Z2_meas = measurement(10.0, 0.2)
        result_mixed2 = calc_parallel_equivalent(Z1_float, Z2_meas)
        expected_mixed2 = 1 / (1 / Z1_float + 1 / Z2_meas)
        @test result_mixed2 isa Measurement{Float64}
        @test isapprox(value(result_mixed2), value(expected_mixed2); atol=TEST_TOL)
        @test isapprox(uncertainty(result_mixed2), uncertainty(expected_mixed2); atol=TEST_TOL)

        # Fully Promoted Case: Both inputs are Measurements
        result_full_meas = calc_parallel_equivalent(Z1_meas, Z2_meas)
        expected_full_meas = 1 / (1 / Z1_meas + 1 / Z2_meas)
        @test result_full_meas isa Measurement{Float64}
        @test isapprox(value(result_full_meas), value(expected_full_meas); atol=TEST_TOL)
        @test isapprox(uncertainty(result_full_meas), uncertainty(expected_full_meas); atol=TEST_TOL)

        # Fully Promoted Complex Case
        Z1_cplx_meas = measurement(3.0, 0.1) + measurement(4.0, 0.2)im
        Z2_cplx_meas = measurement(8.0, 0.3) - measurement(6.0, 0.4)im
        result_cplx_meas = calc_parallel_equivalent(Z1_cplx_meas, Z2_cplx_meas)
        expected_cplx_meas = 1 / (1 / Z1_cplx_meas + 1 / Z2_cplx_meas)
        @test result_cplx_meas isa Complex{Measurement{Float64}}
        @test isapprox(value(real(result_cplx_meas)), value(real(expected_cplx_meas)); atol=TEST_TOL)
        @test isapprox(value(imag(result_cplx_meas)), value(imag(expected_cplx_meas)); atol=TEST_TOL)
        @test isapprox(uncertainty(real(result_cplx_meas)), uncertainty(real(expected_cplx_meas)); atol=TEST_TOL)
        @test isapprox(uncertainty(imag(result_cplx_meas)), uncertainty(imag(expected_cplx_meas)); atol=TEST_TOL)
    end
end