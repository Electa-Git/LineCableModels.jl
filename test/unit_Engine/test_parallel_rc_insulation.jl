@testitem "Engine(ParallelRC): analytical layer and lossless limit" setup = [defaults] begin
    using LinearAlgebra

    formulation = InsulationAdmittance.ParallelRC()
    lossless = InsulationAdmittance.Lossless()

    r_in = 0.010
    r_ex = 0.018
    rho = 2.0e11
    eps_r = 2.4
    s = Complex(0.0, 2π * 50.0)

    log_ratio = log(r_ex / r_in)
    capacitance = 2π * ε₀ * eps_r / log_ratio
    conductance = 2π / (rho * log_ratio)
    p = formulation(r_in, r_ex, rho, eps_r, s)

    @test s / p ≈ conductance + s * capacitance
    @test real(s / p) > 0
    @test imag(s / p) > 0

    # Radially stacked layers are series admittances, hence their potential
    # coefficients add before the solver constructs the Maxwell matrix.
    r_mid = 0.013
    rho_2 = 8.0e10
    eps_r_2 = 3.6
    p_1 = formulation(r_in, r_mid, rho, eps_r, s)
    p_2 = formulation(r_mid, r_ex, rho_2, eps_r_2, s)
    y_1 = s / p_1
    y_2 = s / p_2
    y_series = inv(inv(y_1) + inv(y_2))
    @test s / (p_1 + p_2) ≈ y_series

    @test formulation(r_in, r_ex, Inf, eps_r, s) ≈
          lossless(r_in, r_ex, eps_r, s, 0.0)
end

@testitem "Engine(ParallelRC): physical layers reach the full EMT solve" setup = [defaults] begin
    using LinearAlgebra

    function two_terminal_problem(; uncertain = false)
        copper = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
        rho_1 = uncertain ? measurement(2.0e11, 2.0e10) : 2.0e11
        eps_1 = uncertain ? measurement(2.4, 0.12) : 2.4
        t_1 = uncertain ? measurement(3.0e-3, 1.5e-4) : 3.0e-3
        dielectric_1 = Material(rho_1, eps_1, 1.0, 20.0, 0.0)
        dielectric_2 = Material(8.0e10, 3.6, 1.0, 20.0, 0.0)
        outer_dielectric = Material(1.0e14, 2.3, 1.0, 20.0, 0.0)

        core_conductor = ConductorGroup(Tubular(0.0, Diameter(0.020), copper))
        core_insulation = InsulatorGroup(
            Insulator(core_conductor, Thickness(t_1), dielectric_1),
        )
        core_insulation = add!(
            core_insulation,
            Insulator,
            Thickness(4.0e-3),
            dielectric_2
        )
        core = CableComponent("core", core_conductor, core_insulation)

        sheath_conductor = ConductorGroup(
            Tubular(core_insulation, Thickness(1.0e-3), copper),
        )
        sheath_insulation = InsulatorGroup(
            Insulator(sheath_conductor, Thickness(2.0e-3), outer_dielectric),
        )
        sheath = CableComponent("sheath", sheath_conductor, sheath_insulation)

        design = CableDesign("parallel-rc-test", core)
        design = add!(design, sheath)
        position = CablePosition(
            design,
            0.0,
            -1.0,
            Dict("core" => 1, "sheath" => 2)
        )
        system = LineCableSystem("parallel-rc-test", 1000.0, position)
        frequencies = [1.0e-3, 50.0, 1.0e6]
        earth = EarthModel(frequencies, 100.0, 10.0, 1.0)
        return LineParametersProblem(
            system;
            temperature = 20.0,
            earth_props = earth,
            frequencies = frequencies
        )
    end

    formulation = FormulationSet(
        Val(:EMT);
        insulation_admittance = InsulationAdmittance.ParallelRC(),
        modal_transform = nothing,
        options = (
            reduce_bundle = false,
            kron_reduction = false,
            ideal_transposition = true,
            store_primitive_matrices = true
        )
    )

    problem = two_terminal_problem()
    ws, line_parameters = compute!(problem, formulation)
    pre_refactor = include(joinpath(
        pkgdir(LineCableModels),
        "test",
        "reference",
        "result_fixtures.jl"
    ))

    @test ws.insulator_layer_ranges == [1:2, 3:3]
    @test ws.r_ins_layer_in[2] ≈ ws.r_ins_layer_ext[1]
    @test size(line_parameters.Y) == (2, 2, 3)
    @test line_parameters.Z.values ≈ pre_refactor.deterministic.Z rtol = 8eps()
    @test line_parameters.Y.values ≈ pre_refactor.deterministic.Y rtol = 8eps()

    for k in eachindex(ws.freq)
        s = ws.jω[k]
        expected_p = sum(
            formulation.insulation_admittance(
                ws.r_ins_layer_in[layer_idx],
                ws.r_ins_layer_ext[layer_idx],
                ws.rho_ins_layer[layer_idx],
                ws.eps_ins_layer[layer_idx],
                s
            )
        for layer_idx in ws.insulator_layer_ranges[1]
        )
        @test ws.Pin[1, 1, k] ≈ expected_p

        Y = line_parameters.Y.values[:, :, k]
        @test Y ≈ transpose(Y)
        @test all(isfinite, real.(Y))
        @test all(isfinite, imag.(Y))
        G = real.(Y)
        tolerance = 1.0e-9 * max(opnorm(G), eps())
        @test minimum(eigvals(Symmetric(G))) >= -tolerance
    end

    # Conductive leakage dominates the inner radial branch near DC, whereas
    # capacitive current dominates it at the top of the frequency sweep.
    y_near_dc = ws.jω[1] / ws.Pin[1, 1, 1]
    y_high = ws.jω[end] / ws.Pin[1, 1, end]
    @test real(y_near_dc) > imag(y_near_dc)
    @test imag(y_high) > real(y_high)

    # LEP must survive primitives -> components -> design -> system -> Z/Y.
    ws_uq, line_parameters_uq = compute!(two_terminal_problem(uncertain = true), formulation)
    y_gap_uq = ws_uq.jω[2] / ws_uq.Pin[1, 1, 2]
    @test uncertainty(real(y_gap_uq)) > 0
    @test uncertainty(imag(y_gap_uq)) > 0
    @test uncertainty(real(line_parameters_uq.Y[1, 1, 2])) > 0
    @test uncertainty(imag(line_parameters_uq.Y[1, 1, 2])) > 0
end

@testitem "UQ.mc: ParallelRC is sampled before line-parameter assembly" setup = [defaults] begin
    using Random
    using Statistics
    using LineCableModels.ParametricBuilder:
                                             CableBuilder, Conductor, Insulator, Material,
                                             Earth, SystemBuilder, at
    import LineCableModels.UQ

    copper = Material(
        rho = 1.7241e-8,
        eps_r = 1.0,
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.00393
    )
    dielectric = Material(
        rho = (2.0e11, 10.0),
        eps_r = (2.4, 5.0),
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.0
    )
    outer_dielectric = Material(
        rho = 1.0e14,
        eps_r = 2.3,
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.0
    )

    parts = [
        Conductor.Solid(:core; d = 0.020, m = copper),
        Insulator.Tubular(:core; layers = 1, t = (7.0e-3, 5.0), m = dielectric),
        Conductor.Tubular(:sheath; layers = 1, t = 1.0e-3, m = copper),
        Insulator.Tubular(:sheath; layers = 1, t = 2.0e-3, m = outer_dielectric)
    ]
    builder = CableBuilder("parallel-rc-mc", parts; nominal = NominalData())
    spec = SystemBuilder(
        "parallel-rc-mc",
        builder,
        at(
            x = 0.0,
            y = -1.0,
            phases = (:core => 1, :sheath => 2)
        );
        length = 1000.0,
        temperature = 20.0,
        earth = Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0),
        f = [50.0, 500.0]
    )
    formulation = FormulationSet(
        Val(:EMT);
        insulation_admittance = InsulationAdmittance.ParallelRC(),
        modal_transform = nothing,
        options = (
            reduce_bundle = false,
            kron_reduction = false,
            ideal_transposition = true
        )
    )

    result = UQ.mc(
        spec,
        formulation;
        trials = 12,
        distribution = :normal,
        seed = 20260812,
        print_step = 1000,
        return_samples = true
    )
    pre_refactor = include(joinpath(
        pkgdir(LineCableModels),
        "test",
        "reference",
        "result_fixtures.jl"
    ))

    @test size(surrogate(result).Y) == (2, 2, 2)
    @test size(result.samples.G) == (2, 2, 2, 12)
    @test std(result.samples.G[1, 1, 1, :]) > 0
    @test std(result.samples.C[1, 1, 1, :]) > 0
    @test uncertainty(real(surrogate(result).Y[1, 1, 1])) > 0
    @test uncertainty(imag(surrogate(result).Y[1, 1, 1])) > 0
    for component in (:R, :L, :G, :C)
        retained = getproperty(result.samples, component)
        @test retained[:, :, :, 1] ≈
              getproperty(pre_refactor.monte_carlo.first, component) rtol = 64eps()
        @test retained[:, :, :, end] ≈
              getproperty(pre_refactor.monte_carlo.last, component) rtol = 64eps()
    end

    ω = reshape(2π .* frequencies(result), 1, 1, :)
    Rmeas = real.(surrogate(result).Z.values)
    Lmeas = imag.(surrogate(result).Z.values) ./ ω
    Gmeas = real.(surrogate(result).Y.values)
    Cmeas = imag.(surrogate(result).Y.values) ./ ω
    X = vcat(
        reshape(result.samples.R, :, 12),
        reshape(result.samples.L, :, 12),
        reshape(result.samples.G, :, 12),
        reshape(result.samples.C, :, 12)
    )
    joint = vcat(vec(Rmeas), vec(Lmeas), vec(Gmeas), vec(Cmeas))
    μ = vec(mean(X; dims = 2))
    centered = X .- μ
    empirical_covariance = centered * transpose(centered) / 11

    @test value.(joint) ≈ μ
    @test uncertainty.(joint) ≈ vec(std(X; dims = 2))
    @test Measurements.cov(joint) ≈ empirical_covariance

    for t in (1, 6, 12)
        lp_trial = UQ.trial(result, t)
        expected_Z = result.samples.R[:, :, :, t] .+
                     im .* ω .* result.samples.L[:, :, :, t]
        expected_Y = result.samples.G[:, :, :, t] .+
                     im .* ω .* result.samples.C[:, :, :, t]
        @test domain(lp_trial) === domain(result)
        @test frequencies(lp_trial) == frequencies(result)
        @test basis(lp_trial) == basis(result)
        @test lp_trial.Z.values == expected_Z
        @test lp_trial.Y.values == expected_Y
    end
    @test_throws BoundsError UQ.trial(result, 0)
    @test_throws BoundsError UQ.trial(result, 13)

    index_rng = MersenneTwister(9182)
    adapter_rng = MersenneTwister(9182)
    for _ in 1:5
        t = rand(index_rng, axes(result.samples.R, 4))
        expected = UQ.trial(result, t)
        actual = rand(adapter_rng, result)
        @test actual.Z.values == expected.Z.values
        @test actual.Y.values == expected.Y.values
    end
    @test rand(result) isa LineParameters

    for component in (:R, :L, :G, :C)
        samples = getproperty(result.samples, component)
        stats = getproperty(result.statistics, component)
        for index in CartesianIndices(stats)
            trials = @view samples[index, :]
            @test stats[index].mean ≈ mean(trials)
            @test stats[index].std ≈ std(trials)
        end
    end

    repeated = UQ.mc(
        spec,
        formulation;
        trials = 12,
        distribution = :normal,
        seed = 20260812,
        print_step = 1000,
        return_samples = true
    )
    @test repeated.samples == result.samples
    @test value.(surrogate(repeated).Z.values) == value.(surrogate(result).Z.values)
    @test value.(surrogate(repeated).Y.values) == value.(surrogate(result).Y.values)
    @test Measurements.cov(vcat(
        vec(real.(surrogate(repeated).Z.values)),
        vec(imag.(surrogate(repeated).Z.values) ./ ω),
        vec(real.(surrogate(repeated).Y.values)),
        vec(imag.(surrogate(repeated).Y.values) ./ ω)
    )) ≈ empirical_covariance

    single = UQ.mc(
        spec,
        formulation;
        trials = 1,
        seed = 20260812,
        print_step = 1000
    )
    for component in (:R, :L, :G, :C)
        @test all(iszero(stat.std) for stat in getproperty(single.statistics, component))
    end
    @test all(iszero ∘ uncertainty ∘ real, surrogate(single).Z.values)
    @test all(iszero ∘ uncertainty ∘ imag, surrogate(single).Z.values)
    @test all(iszero ∘ uncertainty ∘ real, surrogate(single).Y.values)
    @test all(iszero ∘ uncertainty ∘ imag, surrogate(single).Y.values)
    @test_throws ArgumentError UQ.trial(single, 1)
    @test_throws ArgumentError rand(MersenneTwister(1), single)

    total = UQ.mc(
        spec,
        formulation;
        trials = 1,
        seed = 20260812,
        print_step = 1000,
        return_samples = true,
        per_length = false
    )
    @test basis(single) === :per_length
    @test basis(total) === :total
    for component in (:R, :L, :G, :C)
        @test mean(total, component, 1, 1, 1) ≈
              1000.0 * mean(single, component, 1, 1, 1)
    end
    @test basis(UQ.trial(total, 1)) === :total
end
