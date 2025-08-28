@testitem "BaseParams: calc_shunt_conductance unit tests" setup = [defaults] begin
    @testset "Basic Functionality" begin
        # Example from docstring
        radius_in = 0.01
        radius_ext = 0.02
        rho = 1e9
        g = calc_shunt_conductance(radius_in, radius_ext, rho)
        @test isapprox(g, 2.7169e-9, atol=TEST_TOL)
        # Lower resistivity increases conductance
        g2 = calc_shunt_conductance(0.01, 0.02, 1e8)
        @test g2 > g
    end

    @testset "Edge Cases" begin
        # Collapsing geometry: radius_in -> radius_ext
        g = calc_shunt_conductance(0.02, 0.02, 1e9)
        @test isinf(g) || isnan(g)
        # Very large radii
        g = calc_shunt_conductance(1e2, 1e3, 1e9)
        @test isfinite(g)
        # Inf/NaN input
        @test isnan(calc_shunt_conductance(NaN, 0.02, 1e9))
        @test isnan(calc_shunt_conductance(0.01, NaN, 1e9))
        @test isnan(calc_shunt_conductance(0.01, 0.02, NaN))
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        g_f32 = calc_shunt_conductance(Float32(0.01), Float32(0.02), Float32(1e9))
        g_f64 = calc_shunt_conductance(0.01, 0.02, 1e9)
        @test isapprox(g_f32, g_f64, atol=TEST_TOL)
    end

    @testset "Physical Behavior" begin
        # Conductance increases as rho decreases
        g1 = calc_shunt_conductance(0.01, 0.02, 1e9)
        g2 = calc_shunt_conductance(0.01, 0.02, 1e8)
        @test g2 > g1
        # Conductance increases as radii get closer
        g3 = calc_shunt_conductance(0.01, 0.011, 1e9)
        @test g3 > g1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        radius_in = 0.01
        radius_ext = 0.02
        rho = 1e9
        min = measurement(radius_in, 1e-4)
        mex = measurement(radius_ext, 1e-4)
        mrho = measurement(rho, 1e7)
        # All Float64
        res1 = calc_shunt_conductance(radius_in, radius_ext, rho)
        @test typeof(res1) == Float64
        # All Measurement
        res2 = calc_shunt_conductance(min, mex, mrho)
        @test res2 isa Measurement{Float64}
        # Mixed: first argument Measurement
        res3 = calc_shunt_conductance(min, radius_ext, rho)
        @test res3 isa Measurement{Float64}
        # Mixed: second argument Measurement
        res4 = calc_shunt_conductance(radius_in, mex, rho)
        @test res4 isa Measurement{Float64}
        # Mixed: third argument Measurement
        res5 = calc_shunt_conductance(radius_in, radius_ext, mrho)
        @test res5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        min = measurement(0.01, 1e-4)
        mex = measurement(0.02, 1e-4)
        mrho = measurement(1e9, 1e7)
        g = calc_shunt_conductance(min, mex, mrho)
        @test g isa Measurement{Float64}
        @test uncertainty(g) > 0
    end
end
