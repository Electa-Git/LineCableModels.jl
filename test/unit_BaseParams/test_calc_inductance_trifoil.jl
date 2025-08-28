# test/unit_BaseParams/test_calc_inductance_trifoil.jl

@testitem "BaseParams: calc_inductance_trifoil unit tests" setup = [defaults] begin

    #=
    ## Test Case Setup
    Parameters are explicitly separated into positional and keyword arguments
    to match the function signature. This makes all test calls clean and robust.
    =#
    const CANONICAL_POS_ARGS = (
        r_in_co=10e-3,
        r_ext_co=15e-3,
        rho_co=1.72e-8,
        mu_r_co=1.0,
        r_in_scr=20e-3,
        r_ext_scr=25e-3,
        rho_scr=2.82e-8,
        mu_r_scr=1.0,
        S=100e-3,
    )

    const CANONICAL_KW_ARGS = (rho_e=100.0, f=50.0)

    @testset "Basic functionality: canonical example" begin
        L = calc_inductance_trifoil(values(CANONICAL_POS_ARGS)...; CANONICAL_KW_ARGS...)
        expected_L = 1.573964832699787e-7 # H/m
        @test L ≈ expected_L atol = TEST_TOL
    end

    @testset "Physical behavior" begin
        L_base = calc_inductance_trifoil(values(CANONICAL_POS_ARGS)...; CANONICAL_KW_ARGS...)

        pos_args_better_screen = merge(CANONICAL_POS_ARGS, (rho_scr=CANONICAL_POS_ARGS.rho_scr / 10,))
        L_better_screen = calc_inductance_trifoil(values(pos_args_better_screen)...; CANONICAL_KW_ARGS...)
        @test L_better_screen < L_base

        pos_args_higher_mu = merge(CANONICAL_POS_ARGS, (mu_r_co=CANONICAL_POS_ARGS.mu_r_co * 2,))
        L_higher_mu = calc_inductance_trifoil(values(pos_args_higher_mu)...; CANONICAL_KW_ARGS...)
        @test L_higher_mu > L_base

        # Override a keyword argument directly in the call
        L_60Hz = calc_inductance_trifoil(values(CANONICAL_POS_ARGS)...; CANONICAL_KW_ARGS..., f=60.0)
        @test L_60Hz < L_base
    end

    @testset "Edge cases" begin

        pos_args_solid_core = merge(CANONICAL_POS_ARGS, (r_in_co=0.0,))
        L_solid_core = calc_inductance_trifoil(values(pos_args_solid_core)...; CANONICAL_KW_ARGS...)
        @test isfinite(L_solid_core)
        @test L_solid_core > 0.0

        pos_args_perfect_screen = merge(CANONICAL_POS_ARGS, (rho_scr=0.0,))
        L_perfect_screen = calc_inductance_trifoil(values(pos_args_perfect_screen)...; CANONICAL_KW_ARGS...)
        @test isfinite(L_perfect_screen)
        @test L_perfect_screen < calc_inductance_trifoil(values(CANONICAL_POS_ARGS)...; CANONICAL_KW_ARGS...)
    end

    @testset "Type stability and promotion with Measurements.jl" begin
        # Base values
        p_pos = CANONICAL_POS_ARGS
        p_kw = CANONICAL_KW_ARGS

        # Create Measurement versions
        p_pos_meas = map(x -> x ± (x * 0.01), p_pos)
        p_kw_meas = map(x -> x ± (x * 0.01), p_kw)

        # Float case
        L_float = calc_inductance_trifoil(values(p_pos)...; p_kw...)
        @test L_float isa Float64

        # Fully promoted case
        L_meas_all = calc_inductance_trifoil(values(p_pos_meas)...; p_kw_meas...)
        @test L_meas_all isa Measurement{Float64}
        @test L_meas_all.val ≈ L_float atol = TEST_TOL
        @test L_meas_all.err > 0.0

        # Mixed case (manual call for clarity)
        L_meas_rho_e = calc_inductance_trifoil(
            p_pos.r_in_co ± 0.0, p_pos.r_ext_co, p_pos.rho_co, p_pos.mu_r_co,
            p_pos.r_in_scr, p_pos.r_ext_scr, p_pos.rho_scr, p_pos.mu_r_scr, p_pos.S;
            rho_e=p_kw.rho_e ± 10.0, f=p_kw.f
        )
        @test L_meas_rho_e isa Measurement{Float64}
        @test L_meas_rho_e.val ≈ L_float atol = TEST_TOL
        @test L_meas_rho_e.err > 0.0
    end
end