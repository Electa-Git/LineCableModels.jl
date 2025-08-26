@testitem "BaseParams: calc_strip_resistance unit tests" setup = [defaults] begin
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
        @test isapprox(R, 0.00017579785649999996, atol=TEST_TOL)
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
        R = calc_strip_resistance(Float32(0.002), Float32(0.05), Float32(1.7241e-8), Float32(0.00393), Float32(20.0), Float32(25.0))
        @test isapprox(R, 0.00017579785649999996, atol=TEST_TOL)
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
        Rm = calc_strip_resistance(measurement(0.002, 1e-6), measurement(0.05, 1e-5), measurement(1.7241e-8, 1e-10), measurement(0.00393, 1e-6), measurement(20.0, 0.1), measurement(25.0, 0.1))
        @test Rm isa Measurement{Float64}
        # Mixed: thickness as Measurement
        R1 = calc_strip_resistance(measurement(0.002, 1e-6), 0.05, 1.7241e-8, 0.00393, 20.0, 25.0)
        @test R1 isa Measurement{Float64}
        # Mixed: alpha as Measurement
        R2 = calc_strip_resistance(0.002, 0.05, 1.7241e-8, measurement(0.00393, 1e-6), 20.0, 25.0)
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
