@testitem "BaseParams: calc_equivalent_alpha unit tests" setup = [defaults] begin

    @testset "calc_equivalent_alpha: Basic Functionality (Copper & Aluminum)" begin
        alpha1 = 0.00393  # Copper
        R1 = 0.5
        alpha2 = 0.00403  # Aluminum
        R2 = 1.0
        expected = (alpha1 * R2 + alpha2 * R1) / (R1 + R2)
        result = calc_equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol=TEST_TOL)
    end

    @testset "calc_equivalent_alpha: Edge Case - Zero Resistance" begin
        alpha1 = 0.00393
        R1 = 0.0
        alpha2 = 0.00403
        R2 = 1.0
        expected = alpha1  # Only R2 matters
        result = calc_equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol=TEST_TOL)

        alpha1 = 0.00393
        R1 = 0.5
        alpha2 = 0.00403
        R2 = 0.0
        expected = alpha2  # Only R1 matters
        result = calc_equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol=TEST_TOL)
    end

    @testset "calc_equivalent_alpha: Edge Case - Very Large Resistance" begin
        alpha1 = 0.00393
        R1 = 1e12
        alpha2 = 0.00403
        R2 = 1.0
        expected = (alpha1 * R2 + alpha2 * R1) / (R1 + R2)
        result = calc_equivalent_alpha(alpha1, R1, alpha2, R2)
        @test isapprox(result, expected; atol=TEST_TOL)
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
        expected_val = (value(alpha1) * value(R2) + value(alpha2) * value(R1)) / (value(R1) + value(R2))
        @test isapprox(value(result), expected_val; atol=TEST_TOL)
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
        @test isapprox(result, expected; atol=TEST_TOL)

        # Edge case: Identical conductors
        alpha = 0.00393
        R = 1.0
        @test isapprox(calc_equivalent_alpha(alpha, R, alpha, R), alpha; atol=TEST_TOL)

        # Edge case: One resistance much larger than the other
        @test isapprox(calc_equivalent_alpha(0.003, 1e6, 0.005, 1.0), 0.005; atol=TEST_TOL)
        @test isapprox(calc_equivalent_alpha(0.003, 1.0, 0.005, 1e6), 0.003; atol=TEST_TOL)

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
            @test isapprox(value(res), expected; atol=TEST_TOL)
            # Uncertainty should be nonzero
            @test uncertainty(res) > 0

            # Mixed case 1: First argument is Measurement
            res = calc_equivalent_alpha(m1, R1, alpha2, R2)
            @test res isa Measurement{Float64}
            @test isapprox(value(res), expected; atol=TEST_TOL)

            # Mixed case 2: Middle argument is Measurement
            res = calc_equivalent_alpha(alpha1, r1, alpha2, R2)
            @test res isa Measurement{Float64}
            @test isapprox(value(res), expected; atol=TEST_TOL)

            # Mixed case 3: Last argument is Measurement
            res = calc_equivalent_alpha(alpha1, R1, alpha2, r2)
            @test res isa Measurement{Float64}
            @test isapprox(value(res), expected; atol=TEST_TOL)
        end

        # Physically unusual but valid: zero resistance (should return NaN)
        @test isnan(calc_equivalent_alpha(0.003, 0.0, 0.005, 0.0))

        # Large values
        @test isapprox(calc_equivalent_alpha(1e-3, 1e6, 2e-3, 2e6), (1e-3 * 2e6 + 2e-3 * 1e6) / (1e6 + 2e6); atol=TEST_TOL)
    end

end # End of test file
