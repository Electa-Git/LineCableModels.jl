@testitem "BaseParams / calc_equivalent_eps / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        # Example from docstring: C_eq=1e-10 F/m, r_ext=0.01 m, r_in=0.005 m
        result = calc_equivalent_eps(1e-10, 0.01, 0.005)
        expected = (1e-10 * log(0.01 / 0.005)) / (2 * pi) / ε₀
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
        @test result > 0
    end

    @testset "Edge Cases" begin
        # Zero capacitance
        result = calc_equivalent_eps(0.0, 0.01, 0.005)
        @test isapprox(result, 0.0; atol = TestNumerics.absolute_floor(Float64))
        # Collapsing geometry: r_ext == r_in
        result = calc_equivalent_eps(1e-10, 0.01, 0.01)
        @test isapprox(result, 0.0; atol = TestNumerics.absolute_floor(Float64))
        # Very large radii
        result = calc_equivalent_eps(1e-10, 1e6, 1e3)
        expected = (1e-10 * log(1e6 / 1e3)) / (2 * pi) / ε₀
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
        # Inf/NaN
        @test isnan(calc_equivalent_eps(NaN, 0.01, 0.005))
        @test isnan(calc_equivalent_eps(1e-10, NaN, 0.005))
        @test isnan(calc_equivalent_eps(1e-10, 0.01, NaN))
        @test isinf(calc_equivalent_eps(Inf, 0.01, 0.005))
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        r = calc_equivalent_eps(Float32(1e-10), Float32(0.01), Float32(0.005))
        d = calc_equivalent_eps(1e-10, 0.01, 0.005)
        @test isapprox(r, d; atol = 1e-6)
    end

    @testset "Physical Behavior" begin
        # Increases with capacitance
        r1 = calc_equivalent_eps(1e-10, 0.01, 0.005)
        r2 = calc_equivalent_eps(2e-10, 0.01, 0.005)
        @test r2 > r1
        # Increases with log(r_ext/r_in)
        r3 = calc_equivalent_eps(1e-10, 0.02, 0.005)
        @test r3 > r1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        # All Float64
        r1 = calc_equivalent_eps(1e-10, 0.01, 0.005)
        @test typeof(r1) == Float64
        # All Measurement
        r2 = calc_equivalent_eps(
            measurement(1e-10, 1e-12),
            measurement(0.01, 1e-5),
            measurement(0.005, 1e-5)
        )
        @test r2 isa Measurement{Float64}
        # Mixed: C_eq as Measurement
        r3 = calc_equivalent_eps(measurement(1e-10, 1e-12), 0.01, 0.005)
        @test r3 isa Measurement{Float64}
        # Mixed: r_ex as Measurement
        r4 = calc_equivalent_eps(1e-10, measurement(0.01, 1e-5), 0.005)
        @test r4 isa Measurement{Float64}
        # Mixed: r_in as Measurement
        r5 = calc_equivalent_eps(1e-10, 0.01, measurement(0.005, 1e-5))
        @test r5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        C_eq = measurement(1e-10, 1e-12)
        r_ext = measurement(0.01, 1e-5)
        r_in = measurement(0.005, 1e-5)
        result = calc_equivalent_eps(C_eq, r_ext, r_in)
        @test result isa Measurement{Float64}
        @test uncertainty(result) > 0
    end
end
@testitem "BaseParams / calc_equivalent_lossfact / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        # Example from docstring: G_eq=1e-8 S·m, C_eq=1e-10 F/m, ω=2π*50
        G_eq = 1e-8
        C_eq = 1e-10
        ω = 2 * pi * 50
        result = calc_equivalent_lossfact(G_eq, C_eq, ω)
        expected = G_eq / (ω * C_eq)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
        @test result > 0
    end

    @testset "Edge Cases" begin
        # Zero conductance
        result = calc_equivalent_lossfact(0.0, 1e-10, 2 * pi * 50)
        @test isapprox(result, 0.0; atol = TestNumerics.absolute_floor(Float64))
        # Zero capacitance (should be Inf)
        result = calc_equivalent_lossfact(1e-8, 0.0, 2 * pi * 50)
        @test isinf(result)
        # Zero frequency (should be Inf)
        result = calc_equivalent_lossfact(1e-8, 1e-10, 0.0)
        @test isinf(result)
        # Inf/NaN
        @test isnan(calc_equivalent_lossfact(NaN, 1e-10, 2 * pi * 50))
        @test isnan(calc_equivalent_lossfact(1e-8, NaN, 2 * pi * 50))
        @test isnan(calc_equivalent_lossfact(1e-8, 1e-10, NaN))
        @test isinf(calc_equivalent_lossfact(Inf, 1e-10, 2 * pi * 50))
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        r = calc_equivalent_lossfact(Float32(1e-8), Float32(1e-10), Float32(2 * pi * 50))
        d = calc_equivalent_lossfact(1e-8, 1e-10, 2 * pi * 50)
        @test isapprox(r, d; atol = 1e-6)
    end

    @testset "Physical Behavior" begin
        # Increases with G_eq
        r1 = calc_equivalent_lossfact(1e-8, 1e-10, 2 * pi * 50)
        r2 = calc_equivalent_lossfact(2e-8, 1e-10, 2 * pi * 50)
        @test r2 > r1
        # Decreases with C_eq
        r3 = calc_equivalent_lossfact(1e-8, 2e-10, 2 * pi * 50)
        @test r3 < r1
        # Decreases with ω
        r4 = calc_equivalent_lossfact(1e-8, 1e-10, 2 * pi * 100)
        @test r4 < r1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        # All Float64
        r1 = calc_equivalent_lossfact(1e-8, 1e-10, 2 * pi * 50)
        @test typeof(r1) == Float64
        # All Measurement
        r2 = calc_equivalent_lossfact(measurement(1e-8, 1e-10), measurement(1e-10, 1e-12),
            measurement(2 * pi * 50, 0.1))
        @test r2 isa Measurement{Float64}
        # Mixed: G_eq as Measurement
        r3 = calc_equivalent_lossfact(measurement(1e-8, 1e-10), 1e-10, 2 * pi * 50)
        @test r3 isa Measurement{Float64}
        # Mixed: C_eq as Measurement
        r4 = calc_equivalent_lossfact(1e-8, measurement(1e-10, 1e-12), 2 * pi * 50)
        @test r4 isa Measurement{Float64}
        # Mixed: ω as Measurement
        r5 = calc_equivalent_lossfact(1e-8, 1e-10, measurement(2 * pi * 50, 0.1))
        @test r5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        G_eq = measurement(1e-8, 1e-10)
        C_eq = measurement(1e-10, 1e-12)
        ω = measurement(2 * pi * 50, 0.1)
        result = calc_equivalent_lossfact(G_eq, C_eq, ω)
        @test result isa Measurement{Float64}
        @test uncertainty(result) > 0
    end
end
@testitem "BaseParams / calc_shunt_capacitance / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        # Example from docstring
        r_in = 0.01
        r_ex = 0.02
        epsr = 2.3
        cap = calc_shunt_capacitance(r_in, r_ex, epsr)
        @test isapprox(cap, 1.241e-10, atol = TestNumerics.absolute_floor(Float64))
        # Vacuum (epsr = 1)
        cap_vac = calc_shunt_capacitance(0.01, 0.02, 1.0)
        @test cap_vac < cap
    end

    @testset "Edge Cases" begin
        # Collapsing geometry: r_in -> r_ex
        cap = calc_shunt_capacitance(0.02, 0.02, 2.3)
        @test isinf(cap) || isnan(cap)
        # Very large radii
        cap = calc_shunt_capacitance(1e2, 1e3, 2.3)
        @test isfinite(cap)
        # Inf/NaN input
        @test isnan(calc_shunt_capacitance(NaN, 0.02, 2.3))
        @test isnan(calc_shunt_capacitance(0.01, NaN, 2.3))
        @test isnan(calc_shunt_capacitance(0.01, 0.02, NaN))
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        cap_f32 = calc_shunt_capacitance(Float32(0.01), Float32(0.02), Float32(2.3))
        cap_f64 = calc_shunt_capacitance(0.01, 0.02, 2.3)
        @test isapprox(cap_f32, cap_f64, atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "Physical Behavior" begin
        # Capacitance increases with epsr
        c1 = calc_shunt_capacitance(0.01, 0.02, 2.3)
        c2 = calc_shunt_capacitance(0.01, 0.02, 3.0)
        @test c2 > c1
        # Capacitance decreases as radii get closer
        c3 = calc_shunt_capacitance(0.01, 0.011, 2.3)
        @test c3 > c1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        r_in = 0.01
        r_ex = 0.02
        epsr = 2.3
        min = measurement(r_in, 1e-4)
        mex = measurement(r_ex, 1e-4)
        mepsr = measurement(epsr, 1e-2)
        # All Float64
        res1 = calc_shunt_capacitance(r_in, r_ex, epsr)
        @test typeof(res1) == Float64
        # All Measurement
        res2 = calc_shunt_capacitance(min, mex, mepsr)
        @test res2 isa Measurement{Float64}
        # Mixed: first argument Measurement
        res3 = calc_shunt_capacitance(min, r_ex, epsr)
        @test res3 isa Measurement{Float64}
        # Mixed: second argument Measurement
        res4 = calc_shunt_capacitance(r_in, mex, epsr)
        @test res4 isa Measurement{Float64}
        # Mixed: third argument Measurement
        res5 = calc_shunt_capacitance(r_in, r_ex, mepsr)
        @test res5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        min = measurement(0.01, 1e-4)
        mex = measurement(0.02, 1e-4)
        mepsr = measurement(2.3, 1e-2)
        cap = calc_shunt_capacitance(min, mex, mepsr)
        @test cap isa Measurement{Float64}
        @test uncertainty(cap) > 0
    end
end
@testitem "BaseParams / calc_shunt_conductance / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        # Example from docstring
        r_in = 0.01
        r_ex = 0.02
        rho = 1e9
        g = calc_shunt_conductance(r_in, r_ex, rho)
        @test isapprox(g, 2.7169e-9, atol = TestNumerics.absolute_floor(Float64))
        # Lower resistivity increases conductance
        g2 = calc_shunt_conductance(0.01, 0.02, 1e8)
        @test g2 > g
    end

    @testset "Edge Cases" begin
        # Collapsing geometry: r_in -> r_ex
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
        @test isapprox(g_f32, g_f64, atol = TestNumerics.absolute_floor(Float64))
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
        r_in = 0.01
        r_ex = 0.02
        rho = 1e9
        min = measurement(r_in, 1e-4)
        mex = measurement(r_ex, 1e-4)
        mrho = measurement(rho, 1e7)
        # All Float64
        res1 = calc_shunt_conductance(r_in, r_ex, rho)
        @test typeof(res1) == Float64
        # All Measurement
        res2 = calc_shunt_conductance(min, mex, mrho)
        @test res2 isa Measurement{Float64}
        # Mixed: first argument Measurement
        res3 = calc_shunt_conductance(min, r_ex, rho)
        @test res3 isa Measurement{Float64}
        # Mixed: second argument Measurement
        res4 = calc_shunt_conductance(r_in, mex, rho)
        @test res4 isa Measurement{Float64}
        # Mixed: third argument Measurement
        res5 = calc_shunt_conductance(r_in, r_ex, mrho)
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
@testitem "BaseParams / calc_sigma_lossfact / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        # Example from docstring: G_eq=2.7169e-9 S·m, r_in=0.01 m, r_ext=0.02 m
        G_eq = 2.7169e-9
        r_in = 0.01
        r_ext = 0.02
        result = calc_sigma_lossfact(G_eq, r_in, r_ext)
        expected = G_eq * log(r_ext / r_in) / (2 * pi)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
        @test result > 0
    end

    @testset "Edge Cases" begin
        # Zero conductance
        result = calc_sigma_lossfact(0.0, 0.01, 0.02)
        @test isapprox(result, 0.0; atol = TestNumerics.absolute_floor(Float64))
        # Collapsing geometry: r_ext == r_in
        result = calc_sigma_lossfact(1e-9, 0.01, 0.01)
        @test isapprox(result, 0.0; atol = TestNumerics.absolute_floor(Float64))
        # Very large radii
        result = calc_sigma_lossfact(1e-9, 1e3, 1e6)
        expected = 1e-9 * log(1e6 / 1e3) / (2 * pi)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
        # Inf/NaN
        @test isnan(calc_sigma_lossfact(NaN, 0.01, 0.02))
        @test isnan(calc_sigma_lossfact(1e-9, NaN, 0.02))
        @test isnan(calc_sigma_lossfact(1e-9, 0.01, NaN))
        @test isinf(calc_sigma_lossfact(Inf, 0.01, 0.02))
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        r = calc_sigma_lossfact(Float32(1e-9), Float32(0.01), Float32(0.02))
        d = calc_sigma_lossfact(1e-9, 0.01, 0.02)
        @test isapprox(r, d; atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "Physical Behavior" begin
        # Increases with G_eq
        r1 = calc_sigma_lossfact(1e-9, 0.01, 0.02)
        r2 = calc_sigma_lossfact(2e-9, 0.01, 0.02)
        @test r2 > r1
        # Increases with log(r_ext/r_in)
        r3 = calc_sigma_lossfact(1e-9, 0.01, 0.04)
        @test r3 > r1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        # All Float64
        r1 = calc_sigma_lossfact(1e-9, 0.01, 0.02)
        @test typeof(r1) == Float64
        # All Measurement
        r2 = calc_sigma_lossfact(
            measurement(1e-9, 1e-11),
            measurement(0.01, 1e-5),
            measurement(0.02, 1e-5)
        )
        @test r2 isa Measurement{Float64}
        # Mixed: G_eq as Measurement
        r3 = calc_sigma_lossfact(measurement(1e-9, 1e-11), 0.01, 0.02)
        @test r3 isa Measurement{Float64}
        # Mixed: r_in as Measurement
        r4 = calc_sigma_lossfact(1e-9, measurement(0.01, 1e-5), 0.02)
        @test r4 isa Measurement{Float64}
        # Mixed: r_ex as Measurement
        r5 = calc_sigma_lossfact(1e-9, 0.01, measurement(0.02, 1e-5))
        @test r5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        G_eq = measurement(1e-9, 1e-11)
        r_in = measurement(0.01, 1e-5)
        r_ext = measurement(0.02, 1e-5)
        result = calc_sigma_lossfact(G_eq, r_in, r_ext)
        @test result isa Measurement{Float64}
        @test uncertainty(result) > 0
    end
end
