@testitem "BaseParams / calc_equivalent_mu / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        # Example from docstring
        gmr = 0.015
        r_ex = 0.02
        r_in = 0.01
        mu_r = calc_equivalent_mu(gmr, r_ex, r_in)
        @test isapprox(mu_r, 1.79409188, atol = TestNumerics.absolute_floor(Float64))

        # Solid conductor (r_in = 0)
        radius_ext_solid = 0.0135
        radius_in_solid = 0.0
        gmr_solid = calc_tubular_gmr(radius_ext_solid, radius_in_solid, 1.0)
        mu_r_solid = calc_equivalent_mu(gmr_solid, radius_ext_solid, radius_in_solid)
        @test isapprox(mu_r_solid, 1.0, atol = TestNumerics.absolute_floor(Float64))
        r_ex = -0.01
        r_in = 0.01
        @test_throws ArgumentError calc_equivalent_mu(gmr, r_ex, r_in)
    end

    @testset "Edge Cases" begin
        # Collapsing geometry: r_in -> r_ex, should be 0 if == gmr
        gmr = 0.02
        r_ex = 0.02
        r_in = 0.02
        mu_r = calc_equivalent_mu(gmr, r_ex, r_in)
        @test isapprox(mu_r, 0.0, atol = TestNumerics.absolute_floor(Float64))

        # Very large radii
        gmr = 1e3
        r_ex = 1e3
        r_in = 1e2
        mu_r = calc_equivalent_mu(gmr, r_ex, r_in)
        @test isfinite(mu_r)

        # Inf/NaN input
        @test isnan(calc_equivalent_mu(NaN, 0.02, 0.01))
        @test isnan(calc_equivalent_mu(0.015, NaN, 0.01))
        @test isnan(calc_equivalent_mu(0.015, 0.02, NaN))
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        gmr = Float32(0.015)
        r_ex = Float32(0.02)
        r_in = Float32(0.01)
        mu_r_f32 = calc_equivalent_mu(gmr, r_ex, r_in)
        mu_r_f64 = calc_equivalent_mu(Float64(gmr), Float64(r_ex), Float64(r_in))
        @test isapprox(mu_r_f32, mu_r_f64; rtol = sqrt(eps(Float32)))
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
        r_ex = 0.02
        r_in = 0.01
        mgmr = measurement(gmr, 1e-4)
        mrex = measurement(r_ex, 1e-4)
        mrin = measurement(r_in, 1e-4)

        # All Float64
        res1 = calc_equivalent_mu(gmr, r_ex, r_in)
        @test typeof(res1) == Float64
        # All Measurement
        res2 = calc_equivalent_mu(mgmr, mrex, mrin)
        @test res2 isa Measurement{Float64}
        # Mixed: first argument Measurement
        res3 = calc_equivalent_mu(mgmr, r_ex, r_in)
        @test res3 isa Measurement{Float64}
        # Mixed: second argument Measurement
        res4 = calc_equivalent_mu(gmr, mrex, r_in)
        @test res4 isa Measurement{Float64}
        # Mixed: third argument Measurement
        res5 = calc_equivalent_mu(gmr, r_ex, mrin)
        @test res5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        gmr = measurement(0.015, 1e-4)
        r_ex = measurement(0.02, 1e-4)
        r_in = measurement(0.01, 1e-4)
        mu_r = calc_equivalent_mu(gmr, r_ex, r_in)
        # Should propagate uncertainty
        @test mu_r isa Measurement{Float64}
        @test uncertainty(mu_r) > 0
    end

    @testset "Error Handling" begin
        # Only error thrown is for r_ex < r_in
        @test_throws ArgumentError calc_equivalent_mu(0.015, 0.01, 0.02)
    end
end
@testitem "BaseParams / calc_tubular_inductance / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    using Measurements
    using TOML

    reference=TOML.parsefile(joinpath(
        pkgdir(LineCableModels),
        "test",
        "fixtures",
        "reference",
        "coaxial_capacitance.toml"
    ))
    # Basic Functionality
    @testset "Basic Functionality" begin
        # Example from docstring: r_in = 0.01, r_ex = 0.02, mu_r = 1.0
        L = calc_tubular_inductance(0.01, 0.02, 1.0)
        expected = 1.0 * μ₀ / (2 * π) * log(0.02 / 0.01)
        @test isapprox(L, expected, atol = TestNumerics.absolute_floor(Float64))
    end

    # Edge Cases
    @testset "Edge Cases" begin
        # Very thin tube (r_ex ≈ r_in)
        r_in = 0.01
        r_ext = 0.010001
        L_thin = calc_tubular_inductance(r_in, r_ext, 1.0)
        expected_thin = 1.0 * μ₀ / (2 * π) * log(r_ext / r_in)
        @test isapprox(L_thin, expected_thin, atol = TestNumerics.absolute_floor(Float64))
        # Large radii
        L_large = calc_tubular_inductance(1e3, 2e3, 1.0)
        expected_large = 1.0 * μ₀ / (2 * π) * log(2e3 / 1e3)
        @test isapprox(L_large, expected_large, atol = TestNumerics.absolute_floor(Float64))
        # mu_r = 0 (non-magnetic)
        @test isapprox(calc_tubular_inductance(0.01, 0.02, 0.0), 0.0,
            atol = TestNumerics.absolute_floor(Float64))
    end

    # Numerical Consistency
    @testset "Numerical Consistency" begin
        # Float32
        Lf = calc_tubular_inductance(Float32(0.01), Float32(0.02), Float32(1.0))
        expectedf = Float32(μ₀) / (2.0f0 * Float32(π)) * log(Float32(0.02) / Float32(0.01))
        @test isapprox(Lf, expectedf, atol = TestNumerics.absolute_floor(Float32))
        # High precision remains available for this transcendental kernel.
        setprecision(BigFloat, 128) do
            Lbig = calc_tubular_inductance(big"0.01", big"0.02", big"1.0")
            @test Lbig isa BigFloat
            @test isapprox(
                Lbig,
                parse(BigFloat, reference["tubular_inductance"]["value"]),
                rtol = 64eps(BigFloat)
            )
        end
    end

    # Physical Behavior
    @testset "Physical Behavior" begin
        # L increases with mu_r
        L1 = calc_tubular_inductance(0.01, 0.02, 1.0)
        L2 = calc_tubular_inductance(0.01, 0.02, 2.0)
        @test L2 > L1
        # L increases with r_ex
        L3 = calc_tubular_inductance(0.01, 0.03, 1.0)
        @test L3 > L1
        # L decreases with r_in
        L4 = calc_tubular_inductance(0.02, 0.03, 1.0)
        @test L4 < L3
    end

    # Type Stability & Promotion
    @testset "Type Stability & Promotion" begin
        # All Float64
        Lf = calc_tubular_inductance(0.01, 0.02, 1.0)
        @test typeof(Lf) == Float64
        # All Measurement
        rinm = measurement(0.01, 1e-5)
        rextm = measurement(0.02, 1e-5)
        murm = measurement(1.0, 1e-3)
        Lm = calc_tubular_inductance(rinm, rextm, murm)
        @test Lm isa Measurement{Float64}
        # Mixed: r_in as Measurement
        Lmix1 = calc_tubular_inductance(rinm, 0.02, 1.0)
        @test Lmix1 isa Measurement{Float64}
        # Mixed: r_ex as Measurement
        Lmix2 = calc_tubular_inductance(0.01, rextm, 1.0)
        @test Lmix2 isa Measurement{Float64}
        # Mixed: mu_r as Measurement
        Lmix3 = calc_tubular_inductance(0.01, 0.02, murm)
        @test Lmix3 isa Measurement{Float64}
    end

    # Uncertainty Quantification
    @testset "Uncertainty Quantification" begin
        rinm = measurement(0.01, 1e-5)
        rextm = measurement(0.02, 1e-5)
        murm = measurement(1.0, 1e-3)
        Lm = calc_tubular_inductance(rinm, rextm, murm)
        # Analytical propagation: L = mu_r * μ₀ / (2π) * log(r_ext/r_in)
        μ = 1.0 * μ₀ / (2 * π) * log(0.02 / 0.01)
        # Partial derivatives
        dL_drin = -murm * μ₀ / (2 * π) * (1 / rinm) / (rextm / rinm)
        dL_drext = murm * μ₀ / (2 * π) * (1 / rextm) / (rextm / rinm)
        dL_dmurm = μ₀ / (2 * π) * log(0.02 / 0.01)
        σ2 = (value(dL_drin) * uncertainty(rinm))^2 +
             (value(dL_drext) * uncertainty(rextm))^2 +
             (value(dL_dmurm) * uncertainty(murm))^2
        @test isapprox(value(Lm), μ, atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(uncertainty(Lm), sqrt(σ2), atol = TestNumerics.absolute_floor(Float64))
    end
end
# test/unit_BaseParams/test_calc_inductance_trifoil.jl

@testitem "BaseParams / calc_inductance_trifoil / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    #=
    ## Test Case Setup
    Parameters are explicitly separated into positional and keyword arguments
    to match the function signature. This makes all test calls clean and robust.
    =#
    const REFERENCE_POS_ARGS=(
        r_in_co = 10e-3,
        r_ext_co = 15e-3,
        rho_co = 1.72e-8,
        mu_r_co = 1.0,
        r_in_scr = 20e-3,
        r_ext_scr = 25e-3,
        rho_scr = 2.82e-8,
        mu_r_scr = 1.0,
        S = 100e-3
    )

    const REFERENCE_KW_ARGS=(rho_e = 100.0, f = 50.0)

    @testset "Basic functionality: reference example" begin
        L = calc_inductance_trifoil(values(REFERENCE_POS_ARGS)...; REFERENCE_KW_ARGS...)
        expected_L = 1.573964832699787e-7 # H/m
        @test L ≈ expected_L atol = TestNumerics.absolute_floor(Float64)
    end

    @testset "Physical behavior" begin
        L_base = calc_inductance_trifoil(values(REFERENCE_POS_ARGS)...; REFERENCE_KW_ARGS...)

        pos_args_better_screen = merge(REFERENCE_POS_ARGS, (rho_scr = REFERENCE_POS_ARGS.rho_scr /
                                                                      10,))
        L_better_screen = calc_inductance_trifoil(values(pos_args_better_screen)...; REFERENCE_KW_ARGS...)
        @test L_better_screen < L_base

        pos_args_higher_mu = merge(REFERENCE_POS_ARGS, (mu_r_co = REFERENCE_POS_ARGS.mu_r_co *
                                                                  2,))
        L_higher_mu = calc_inductance_trifoil(values(pos_args_higher_mu)...; REFERENCE_KW_ARGS...)
        @test L_higher_mu > L_base

        # Override a keyword argument directly in the call
        L_60Hz = calc_inductance_trifoil(values(REFERENCE_POS_ARGS)...; REFERENCE_KW_ARGS..., f = 60.0)
        @test L_60Hz < L_base
    end

    @testset "Edge cases" begin
        pos_args_solid_core = merge(REFERENCE_POS_ARGS, (r_in_co = 0.0,))
        L_solid_core = calc_inductance_trifoil(values(pos_args_solid_core)...; REFERENCE_KW_ARGS...)
        @test isfinite(L_solid_core)
        @test L_solid_core > 0.0

        pos_args_perfect_screen = merge(REFERENCE_POS_ARGS, (rho_scr = 0.0,))
        L_perfect_screen = calc_inductance_trifoil(values(pos_args_perfect_screen)...; REFERENCE_KW_ARGS...)
        @test isfinite(L_perfect_screen)
        @test L_perfect_screen <
              calc_inductance_trifoil(values(REFERENCE_POS_ARGS)...; REFERENCE_KW_ARGS...)
    end

    @testset "Type stability and promotion with Measurements.jl" begin
        # Base values
        p_pos = REFERENCE_POS_ARGS
        p_kw = REFERENCE_KW_ARGS

        # Create Measurement versions
        p_pos_meas = map(x -> x ± (x * 0.01), p_pos)
        p_kw_meas = map(x -> x ± (x * 0.01), p_kw)

        # Float case
        L_float = calc_inductance_trifoil(values(p_pos)...; p_kw...)
        @test L_float isa Float64

        # Fully promoted case
        L_meas_all = calc_inductance_trifoil(values(p_pos_meas)...; p_kw_meas...)
        @test L_meas_all isa Measurement{Float64}
        @test L_meas_all.val ≈ L_float atol = TestNumerics.absolute_floor(Float64)
        @test L_meas_all.err > 0.0

        # Mixed case (manual call for clarity)
        L_meas_rho_e = calc_inductance_trifoil(
            p_pos.r_in_co ± 0.0, p_pos.r_ext_co, p_pos.rho_co, p_pos.mu_r_co,
            p_pos.r_in_scr, p_pos.r_ext_scr, p_pos.rho_scr, p_pos.mu_r_scr, p_pos.S;
            rho_e = p_kw.rho_e ± 10.0, f = p_kw.f
        )
        @test L_meas_rho_e isa Measurement{Float64}
        @test L_meas_rho_e.val ≈ L_float atol = TestNumerics.absolute_floor(Float64)
        @test L_meas_rho_e.err > 0.0
    end
end
