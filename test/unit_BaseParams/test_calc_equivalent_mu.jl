@testitem "BaseParams: calc_equivalent_mu unit tests" setup = [defaults] begin
    @testset "Basic Functionality" begin
        # Example from docstring
        gmr = 0.015
        radius_ext = 0.02
        radius_in = 0.01
        mu_r = calc_equivalent_mu(gmr, radius_ext, radius_in)
        @test isapprox(mu_r, 1.79409188, atol=TEST_TOL)

        # Solid conductor (radius_in = 0)
        radius_ext_solid = 0.0135
        radius_in_solid = 0.0
        gmr_solid = calc_tubular_gmr(radius_ext_solid, radius_in_solid, 1.0)
        mu_r_solid = calc_equivalent_mu(gmr_solid, radius_ext_solid, radius_in_solid)
        @test isapprox(mu_r_solid, 1.0, atol=TEST_TOL)
        radius_ext = -0.01
        radius_in = 0.01
        @test_throws ArgumentError calc_equivalent_mu(gmr, radius_ext, radius_in)
    end

    @testset "Edge Cases" begin
        # Collapsing geometry: radius_in -> radius_ext, should be 0 if == gmr
        gmr = 0.02
        radius_ext = 0.02
        radius_in = 0.02
        mu_r = calc_equivalent_mu(gmr, radius_ext, radius_in)
        @test isapprox(mu_r, 0.0, atol=TEST_TOL)

        # Very large radii
        gmr = 1e3
        radius_ext = 1e3
        radius_in = 1e2
        mu_r = calc_equivalent_mu(gmr, radius_ext, radius_in)
        @test isfinite(mu_r)

        # Inf/NaN input
        @test isnan(calc_equivalent_mu(NaN, 0.02, 0.01))
        @test isnan(calc_equivalent_mu(0.015, NaN, 0.01))
        @test isnan(calc_equivalent_mu(0.015, 0.02, NaN))
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        gmr = Float32(0.015)
        radius_ext = Float32(0.02)
        radius_in = Float32(0.01)
        mu_r_f32 = calc_equivalent_mu(gmr, radius_ext, radius_in)
        mu_r_f64 = calc_equivalent_mu(Float64(gmr), Float64(radius_ext), Float64(radius_in))
        @test isapprox(mu_r_f32, mu_r_f64, atol=TEST_TOL)
    end

    @testset "Physical Behavior" begin
        # mu_r increases as gmr decreases (for fixed radii)
        mu1 = calc_equivalent_mu(0.015, 0.02, 0.01)
        mu2 = calc_equivalent_mu(0.012, 0.02, 0.01)
        @test mu2 > mu1
        # mu_r decreases as gmr increases
        mu3 = calc_equivalent_mu(0.018, 0.02, 0.01)
        @test mu3 < mu1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        gmr = 0.015
        radius_ext = 0.02
        radius_in = 0.01
        mgmr = measurement(gmr, 1e-4)
        mrex = measurement(radius_ext, 1e-4)
        mrin = measurement(radius_in, 1e-4)

        # All Float64
        res1 = calc_equivalent_mu(gmr, radius_ext, radius_in)
        @test typeof(res1) == Float64
        # All Measurement
        res2 = calc_equivalent_mu(mgmr, mrex, mrin)
        @test res2 isa Measurement{Float64}
        # Mixed: first argument Measurement
        res3 = calc_equivalent_mu(mgmr, radius_ext, radius_in)
        @test res3 isa Measurement{Float64}
        # Mixed: second argument Measurement
        res4 = calc_equivalent_mu(gmr, mrex, radius_in)
        @test res4 isa Measurement{Float64}
        # Mixed: third argument Measurement
        res5 = calc_equivalent_mu(gmr, radius_ext, mrin)
        @test res5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        gmr = measurement(0.015, 1e-4)
        radius_ext = measurement(0.02, 1e-4)
        radius_in = measurement(0.01, 1e-4)
        mu_r = calc_equivalent_mu(gmr, radius_ext, radius_in)
        # Should propagate uncertainty
        @test mu_r isa Measurement{Float64}
        @test uncertainty(mu_r) > 0
    end

    @testset "Error Handling" begin
        # Only error thrown is for radius_ext < radius_in
        @test_throws ArgumentError calc_equivalent_mu(0.015, 0.01, 0.02)
    end
end
