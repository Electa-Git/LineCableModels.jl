@testitem "BaseParams: calc_solenoid_correction unit tests" setup = [defaults] begin
    @testset "Basic Functionality" begin
        # Example from docstring: 10 turns/m, conductor radius 5 mm, insulator radius 10 mm
        result = calc_solenoid_correction(10.0, 0.005, 0.01)
        expected = 1.0 + 2 * 10.0^2 * pi^2 * (0.01^2 - 0.005^2) / log(0.01 / 0.005)
        @test isapprox(result, expected; atol=TEST_TOL)
        @test result > 1.0

        # Non-helical cable (NaN turns)
        result = calc_solenoid_correction(NaN, 0.005, 0.01)
        @test result == 1.0
    end

    @testset "Edge Cases" begin
        # Zero turns (should be 1.0)
        result = calc_solenoid_correction(0.0, 0.005, 0.01)
        @test isapprox(result, 1.0; atol=TEST_TOL)

        # Collapsing geometry: radii nearly equal
        result = calc_solenoid_correction(10.0, 0.01, 0.010001)
        @test result > 1.0

        # Very large number of turns
        result = calc_solenoid_correction(1e6, 0.005, 0.01)
        @test result > 1e9

        # Inf/NaN radii
        @test isnan(calc_solenoid_correction(10.0, NaN, 0.01))
        @test isnan(calc_solenoid_correction(10.0, 0.005, NaN))
        @test isinf(calc_solenoid_correction(10.0, 0.0, 0.01)) == false  # log(0.01/0) = Inf, but numerator is finite
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        r = calc_solenoid_correction(Float32(10.0), Float32(0.005), Float32(0.01))
        d = calc_solenoid_correction(10.0, 0.005, 0.01)
        @test isapprox(r, d; atol=TEST_TOL)
    end

    @testset "Physical Behavior" begin
        # Correction increases with more turns
        c1 = calc_solenoid_correction(5.0, 0.005, 0.01)
        c2 = calc_solenoid_correction(10.0, 0.005, 0.01)
        @test c2 > c1
        # Correction increases with larger insulator radius
        c3 = calc_solenoid_correction(10.0, 0.005, 0.02)
        @test c3 > c2
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        # All Float64
        r1 = calc_solenoid_correction(10.0, 0.005, 0.01)
        @test typeof(r1) == Float64
        # All Measurement
        r2 = calc_solenoid_correction(measurement(10.0, 0.1), measurement(0.005, 1e-5), measurement(0.01, 1e-5))
        @test r2 isa Measurement{Float64}
        # Mixed: num_turns as Measurement
        r3 = calc_solenoid_correction(measurement(10.0, 0.1), 0.005, 0.01)
        @test r3 isa Measurement{Float64}
        # Mixed: radius_ext_con as Measurement
        r4 = calc_solenoid_correction(10.0, measurement(0.005, 1e-5), 0.01)
        @test r4 isa Measurement{Float64}
        # Mixed: radius_ext_ins as Measurement
        r5 = calc_solenoid_correction(10.0, 0.005, measurement(0.01, 1e-5))
        @test r5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        num_turns = measurement(10.0, 0.1)
        radius_ext_con = measurement(0.005, 1e-5)
        radius_ext_ins = measurement(0.01, 1e-5)
        result = calc_solenoid_correction(num_turns, radius_ext_con, radius_ext_ins)
        # Should propagate uncertainty
        @test result isa Measurement{Float64}
        # Uncertainty should be nonzero
        @test uncertainty(result) > 0
    end
end
