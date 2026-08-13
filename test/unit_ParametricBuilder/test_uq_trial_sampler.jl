@testitem "UQ.mc: custom trial sampler preserves shared primitive draws" setup = [defaults] begin
    using Random
    using LineCableModels.ParametricBuilder:
                                             CableBuilder, Conductor, Insulator, Material,
                                             build
    using LineCableModels.UQ

    conductor = Material(rho = 1.7241e-8, eps_r = 1.0, mu_r = 1.0,
        T0 = 20.0, alpha = 0.00393)
    insulation = Material(rho = 1.0e14, eps_r = 2.3, mu_r = 1.0,
        T0 = 20.0, alpha = 0.0)
    base = CableBuilder(
        "correlated-trial-sampler",
        Conductor.Solid(:core; d = (0.02, 10.0), m = conductor),
        Insulator.Tubular(:core; layers = 1, t = (0.005, 10.0), m = insulation),
        Conductor.Tubular(:sheath; layers = 1, t = 0.001, m = conductor),
        Insulator.Tubular(:sheath; layers = 1, t = 0.001, m = insulation),
        nominal = NominalData(designation_code = "correlated-trial-sampler")
    )

    shared_draws = Tuple{Float64, Float64}[]
    function correlated_sampler(_, trial, distribution)
        distribution == :normal || error("unexpected distribution")
        u = randn()
        d = 0.02 * (1 + 0.01u)
        t = 0.005 * (1 + 0.01u)
        push!(shared_draws, (d, t))
        trial_builder = CableBuilder(
            "correlated-trial-$trial",
            Conductor.Solid(:core; d = d, m = conductor),
            Insulator.Tubular(:core; layers = 1, t = t, m = insulation),
            Conductor.Tubular(:sheath; layers = 1, t = 0.001, m = conductor),
            Insulator.Tubular(:sheath; layers = 1, t = 0.001, m = insulation),
            nominal = NominalData(designation_code = "correlated-trial-$trial")
        )
        return only(build(trial_builder))
    end

    result = UQ.mc(
        base;
        trials = 8,
        seed = 2048,
        trial_sampler = correlated_sampler,
        return_samples = true,
        print_step = 100
    )

    @test length(shared_draws) == 8
    @test result.samples !== nothing
    @test result.stats.R.n == 8
    @test all(isapprox(d / 0.02, t / 0.005; atol = 32eps()) for (d, t) in shared_draws)
    @test length(unique(first.(shared_draws))) > 1
end
