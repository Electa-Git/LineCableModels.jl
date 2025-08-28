@testitem "BaseParams: calc_wirearray_gmr unit tests" setup = [defaults] begin
    using Measurements: measurement, value, uncertainty

    @testset "Basic Functionality" begin
        # Example from docstring
        lay_rad = 0.05
        N = 7
        rad_wire = 0.002
        mu_r = 1.0
        gmr = calc_wirearray_gmr(lay_rad, N, rad_wire, mu_r)
        expected = exp((log(rad_wire * exp(-mu_r / 4) * N * lay_rad^(N - 1)) / N))
        @test isapprox(gmr, expected; atol=TEST_TOL)
        @test gmr > 0
    end

    @testset "Edge Cases" begin
        # N = 1 (single wire)
        lay_rad = 0.05
        N = 1
        rad_wire = 0.002
        mu_r = 1.0
        gmr = calc_wirearray_gmr(lay_rad, N, rad_wire, mu_r)
        expected = rad_wire * exp(-mu_r / 4)
        @test isapprox(gmr, expected; atol=TEST_TOL)

        # mu_r = 0 (non-magnetic)
        lay_rad = 0.05
        N = 7
        rad_wire = 0.002
        mu_r = 0.0
        gmr = calc_wirearray_gmr(lay_rad, N, rad_wire, mu_r)
        expected = exp((log(rad_wire * N * lay_rad^(N - 1)) / N))
        @test isapprox(gmr, expected; atol=TEST_TOL)

        # rad_wire = 0 (degenerate wire)
        lay_rad = 0.05
        N = 7
        rad_wire = 0.0
        mu_r = 1.0
        gmr = calc_wirearray_gmr(lay_rad, N, rad_wire, mu_r)
        @test gmr == 0.0

        # lay_rad = 0 (all wires at center)
        lay_rad = 0.0
        N = 7
        rad_wire = 0.002
        mu_r = 1.0
        gmr = calc_wirearray_gmr(lay_rad, N, rad_wire, mu_r)
        expected = exp((log(rad_wire * exp(-mu_r / 4) * N * 0.0^(N - 1)) / N))
        @test gmr == 0.0
    end

    @testset "Numerical Consistency" begin
        # Float64
        gmr1 = calc_wirearray_gmr(0.05, 7, 0.002, 1.0)
        # Measurement{Float64}
        gmr2 = calc_wirearray_gmr(measurement(0.05, 1e-4), 7, 0.002, 1.0)
        @test isapprox(value(gmr2), gmr1; atol=TEST_TOL)
        @test uncertainty(gmr2) > 0
    end

    @testset "Physical Behavior" begin
        # GMR increases with lay_rad
        gmr1 = calc_wirearray_gmr(0.01, 7, 0.002, 1.0)
        gmr2 = calc_wirearray_gmr(0.05, 7, 0.002, 1.0)
        @test gmr2 > gmr1
        # GMR decreases with mu_r
        gmr1 = calc_wirearray_gmr(0.05, 7, 0.002, 0.5)
        gmr2 = calc_wirearray_gmr(0.05, 7, 0.002, 2.0)
        @test gmr2 < gmr1
    end

    @testset "Type Stability & Promotion" begin
        # All Float64
        gmr = calc_wirearray_gmr(0.05, 7, 0.002, 1.0)
        @test typeof(gmr) == Float64
        # All Measurement
        gmr = calc_wirearray_gmr(measurement(0.05, 1e-4), 7, measurement(0.002, 1e-5), measurement(1.0, 0.01))
        @test gmr isa Measurement{Float64}
        # Mixed: lay_rad as Measurement
        gmr = calc_wirearray_gmr(measurement(0.05, 1e-4), 7, 0.002, 1.0)
        @test gmr isa Measurement{Float64}
        # Mixed: rad_wire as Measurement
        gmr = calc_wirearray_gmr(0.05, 7, measurement(0.002, 1e-5), 1.0)
        @test gmr isa Measurement{Float64}
        # Mixed: mu_r as Measurement
        gmr = calc_wirearray_gmr(0.05, 7, 0.002, measurement(1.0, 0.01))
        @test gmr isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        lay_rad = measurement(0.05, 1e-4)
        rad_wire = measurement(0.002, 1e-5)
        mu_r = measurement(1.0, 0.01)
        gmr = calc_wirearray_gmr(lay_rad, 7, rad_wire, mu_r)
        @test gmr isa Measurement{Float64}
        @test uncertainty(gmr) > 0
    end
end
