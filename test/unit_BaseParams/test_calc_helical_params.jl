@testitem "BaseParams: calc_helical_params unit tests" setup = [defaults] begin
    using Measurements
    # --- Basic Functionality ---
    @testset "Basic Functionality" begin
        radius_in = 0.01
        radius_ext = 0.015
        lay_ratio = 12.0
        mean_diam, pitch, overlength = calc_helical_params(radius_in, radius_ext, lay_ratio)
        @test isapprox(mean_diam, 0.025, atol=TEST_TOL)
        @test isapprox(pitch, 0.3, atol=TEST_TOL)
        @test overlength > 1.0
    end

    # --- Edge Cases ---
    @testset "Edge Cases" begin
        # Zero lay ratio (pitch_length = 0)
        m, p, o = calc_helical_params(0.01, 0.015, 0.0)
        @test isapprox(m, 0.025, atol=TEST_TOL)
        @test isapprox(p, 0.0, atol=TEST_TOL)
        @test isapprox(o, 1.0, atol=TEST_TOL)

        # Collapsing geometry (radius_in == radius_ext)
        m2, p2, o2 = calc_helical_params(0.02, 0.02, 10.0)
        @test isapprox(m2, 0.04, atol=TEST_TOL)
        @test isapprox(p2, 0.4, atol=TEST_TOL)
        @test o2 > 1.0

        # Very large lay ratio
        m3, p3, o3 = calc_helical_params(0.01, 0.015, 1e6)
        @test isapprox(m3, 0.025, atol=TEST_TOL)
        @test isapprox(p3, 25000.0, atol=TEST_TOL)
        @test isapprox(o3, 1.0, atol=TEST_TOL)
    end

    # --- Numerical Consistency ---
    @testset "Numerical Consistency" begin
        # Float32
        m, p, o = calc_helical_params(Float32(0.01), Float32(0.015), Float32(12.0))
        @test isapprox(m, 0.025, atol=TEST_TOL)
        @test isapprox(p, 0.3, atol=TEST_TOL)
        @test o > 1.0
    end

    # --- Physical Behavior ---
    @testset "Physical Behavior" begin
        # Increasing lay_ratio increases pitch_length
        _, p1, _ = calc_helical_params(0.01, 0.015, 10.0)
        _, p2, _ = calc_helical_params(0.01, 0.015, 20.0)
        @test p2 > p1
        # Overlength approaches 1 as lay_ratio increases
        _, _, o1 = calc_helical_params(0.01, 0.015, 1e3)
        @test isapprox(o1, 1.0, atol=1e-5)
    end

    # --- Type Stability & Promotion ---
    @testset "Type Stability & Promotion" begin
        # All Float64
        m, p, o = calc_helical_params(0.01, 0.015, 12.0)
        @test typeof(m) == Float64
        @test typeof(p) == Float64
        @test typeof(o) == Float64
        # All Measurement
        mM, pM, oM = calc_helical_params(measurement(0.01, 1e-5), measurement(0.015, 1e-5), measurement(12.0, 0.1))
        @test mM isa Measurement{Float64}
        @test pM isa Measurement{Float64}
        @test oM isa Measurement{Float64}
        # Mixed: radius_in as Measurement
        m1, p1, o1 = calc_helical_params(measurement(0.01, 1e-5), 0.015, 12.0)
        @test m1 isa Measurement{Float64}
        @test p1 isa Measurement{Float64}
        @test o1 isa Measurement{Float64}
        # Mixed: lay_ratio as Measurement
        m2, p2, o2 = calc_helical_params(0.01, 0.015, measurement(12.0, 0.1))
        @test m2 isa Measurement{Float64}
        @test p2 isa Measurement{Float64}
        @test o2 isa Measurement{Float64}
    end

    # --- Uncertainty Quantification ---
    @testset "Uncertainty Quantification" begin
        rin = measurement(0.01, 1e-5)
        rext = measurement(0.015, 1e-5)
        lrat = measurement(12.0, 0.1)
        m, p, o = calc_helical_params(rin, rext, lrat)
        # Check propagated uncertainties are nonzero
        @test uncertainty(m) > 0
        @test uncertainty(p) > 0
        @test uncertainty(o) > 0
    end
end
