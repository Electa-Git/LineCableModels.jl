@testitem "Engine / ParallelRC / analytical layer and lossless limit" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using TOML

    formulation=InsulationAdmittance.ParallelRC()
    lossless=InsulationAdmittance.Lossless()

    reference=TOML.parsefile(joinpath(
        pkgdir(LineCableModels),
        "test",
        "fixtures",
        "reference",
        "coaxial_capacitance.toml"
    ))
    for T in (Float32, Float64, BigFloat)
        setprecision(BigFloat, 128) do
            typed(value) = parse(T, value)
            r_inner=typed("0.01")
            r_outer=typed("0.02")
            resistivity=typed("2.0e11")
            relative_permittivity=typed("2.3")
            angular_frequency=T(2)*T(π)*T(50)
            frequency_point=Complex{T}(zero(T), angular_frequency)
            coefficient=formulation(
                r_inner,
                r_outer,
                resistivity,
                relative_permittivity,
                frequency_point
            )
            recovered_admittance=frequency_point/coefficient
            expected_conductance=T(parse(
                BigFloat,
                reference["parallel_rc"]["conductance"]
            ))
            expected_capacitance=T(parse(BigFloat, reference["value"]))
            @test TestNumerics.isapprox_scaled(
                real(recovered_admittance),
                expected_conductance
            )
            @test TestNumerics.isapprox_scaled(
                imag(recovered_admittance)/angular_frequency,
                expected_capacitance
            )
        end
    end

    r_in=0.010
    r_ex=0.018
    rho=2.0e11
    eps_r=2.4
    s=Complex(0.0, 2π*50.0)
    epsilon0=8.8541878128e-12

    log_ratio=log(r_ex/r_in)
    capacitance=2π*epsilon0*eps_r/log_ratio
    conductance=2π/(rho*log_ratio)
    coefficient=formulation(r_in, r_ex, rho, eps_r, s)

    @test s / coefficient ≈ conductance + s * capacitance
    @test real(s / coefficient) > 0
    @test imag(s / coefficient) > 0

    r_mid=0.013
    second_rho=8.0e10
    second_eps=3.6
    first_layer=formulation(r_in, r_mid, rho, eps_r, s)
    second_layer=formulation(r_mid, r_ex, second_rho, second_eps, s)
    first_admittance=s/first_layer
    second_admittance=s/second_layer
    series_admittance=inv(inv(first_admittance)+inv(second_admittance))
    @test s / (first_layer + second_layer) ≈ series_admittance

    @test formulation(r_in, r_ex, Inf, eps_r, s) ≈
          lossless(r_in, r_ex, eps_r, s, 0.0)
end

@testitem "Engine / ParallelRC / strict layers reach direct compute!" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using LinearAlgebra

    function two_terminal_problem(; uncertain = false)
        copper=Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
        rho_1=uncertain ? measurement(2.0e11, 2.0e10) : 2.0e11
        eps_1=uncertain ? measurement(2.4, 0.12) : 2.4
        thickness_1=uncertain ? measurement(3.0e-3, 1.5e-4) : 3.0e-3
        dielectric_1=Material(rho_1, eps_1, 1.0, 20.0, 0.0)
        dielectric_2=Material(8.0e10, 3.6, 1.0, 20.0, 0.0)
        outer_dielectric=Material(1.0e14, 2.3, 1.0, 20.0, 0.0)

        core_conductor=ConductorGroup(Tubular(0.0, 0.010, copper))
        core_insulation=InsulatorGroup(Insulator(
            core_conductor.r_ex,
            core_conductor.r_ex+thickness_1,
            dielectric_1
        ))
        core_insulation=add!(
            core_insulation,
            Insulator(core_insulation.r_ex, core_insulation.r_ex+4.0e-3, dielectric_2)
        )
        core=CableComponent("core", core_conductor, core_insulation)

        sheath_conductor=ConductorGroup(Tubular(
            core_insulation.r_ex,
            core_insulation.r_ex+1.0e-3,
            copper
        ))
        sheath_insulation=InsulatorGroup(Insulator(
            sheath_conductor.r_ex,
            sheath_conductor.r_ex+2.0e-3,
            outer_dielectric
        ))
        sheath=CableComponent("sheath", sheath_conductor, sheath_insulation)

        design=add!(CableDesign("parallel-rc-test", core), sheath)
        position=CablePosition(
            design,
            0.0,
            -1.0,
            Dict("core"=>1, "sheath"=>2)
        )
        system=LineCableSystem("parallel-rc-test", 1000.0, position)
        frequency=[1.0e-3, 50.0, 1.0e6]
        earth=EarthModel(100.0, 10.0, 1.0)
        return LineParametersProblem(
            system;
            temperature = 20.0,
            earth_props = earth,
            frequencies = frequency
        )
    end

    formulation=Formulation(
        Val(:EMT);
        insulation_admittance = InsulationAdmittance.ParallelRC(),
        modal_transform = nothing,
        options = (
            reduce_bundle = false,
            kron_reduction = false,
            ideal_transposition = true
        )
    )
    problem=two_terminal_problem()
    input=Engine.EMTInput(problem)
    trace=compute!(problem, formulation; inspect = true)
    parameters=trace.result
    public_parameters=compute!(problem, formulation)
    @test public_parameters.Z.values == trace.result.Z.values
    @test public_parameters.Y.values == trace.result.Y.values
    @test input.insulator_layer_ranges == [1:2, 3:3]
    @test input.r_ins_layer_in[2] ≈ input.r_ins_layer_ext[1]
    @test size(parameters.Y) == (2, 2, 3)

    for frequency_index in eachindex(input.freq)
        s=input.jω[frequency_index]
        inner=sum(
            formulation.insulation_admittance(
                input.r_ins_layer_in[layer_index],
                input.r_ins_layer_ext[layer_index],
                input.rho_ins_layer[layer_index],
                input.eps_ins_layer[layer_index],
                s
            )
        for layer_index in input.insulator_layer_ranges[1]
        )
        outer=sum(
            formulation.insulation_admittance(
                input.r_ins_layer_in[layer_index],
                input.r_ins_layer_ext[layer_index],
                input.rho_ins_layer[layer_index],
                input.eps_ins_layer[layer_index],
                s
            )
        for layer_index in input.insulator_layer_ranges[2]
        )
        @test trace.Pin[:, :, frequency_index] ≈ [inner+outer outer
                                                  outer outer]

        admittance=parameters.Y.values[:, :, frequency_index]
        @test admittance ≈ transpose(admittance)
        @test all(isfinite, real.(admittance))
        @test all(isfinite, imag.(admittance))
        conductance=real.(admittance)
        tolerance=1.0e-9*max(opnorm(conductance), eps())
        @test minimum(eigvals(Symmetric(conductance))) >= -tolerance
    end

    function component_coefficient(component_index, frequency_index)
        s=input.jω[frequency_index]
        return sum(
            formulation.insulation_admittance(
                input.r_ins_layer_in[layer_index],
                input.r_ins_layer_ext[layer_index],
                input.rho_ins_layer[layer_index],
                input.eps_ins_layer[layer_index],
                s
            )
        for layer_index in input.insulator_layer_ranges[component_index]
        )
    end
    near_dc=input.jω[1]/component_coefficient(1, 1)
    high_frequency=input.jω[end]/component_coefficient(1, length(input.jω))
    @test real(near_dc) > imag(near_dc)
    @test imag(high_frequency) > real(high_frequency)

    uncertain_parameters=compute!(two_terminal_problem(uncertain = true), formulation)
    @test uncertainty(real(uncertain_parameters.Y[1, 1, 2])) > 0
    @test uncertainty(imag(uncertain_parameters.Y[1, 1, 2])) > 0

    total=compute!(
        problem,
        formulation;
        options = (output_basis = :total,)
    )
    @test basis(parameters) === :per_length
    @test basis(total) === :total
    @test total.Z.values ≈ problem.system.line_length .* parameters.Z.values
    @test total.Y.values ≈ problem.system.line_length .* parameters.Y.values
end

@testitem "Engine / ParallelRC / Gridspace Monte Carlo samples before assembly" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using Statistics
    import LineCableModels.ParametricBuilder as PB

    copper=PB.Material(
        rho = 1.7241e-8,
        eps_r = 1.0,
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.00393
    )
    dielectric=PB.Material(
        rho = Grid(2.0e11, 10.0),
        eps_r = Grid(2.4, 5.0)
    )
    outer_dielectric=PB.Material(rho = 1.0e14, eps_r = 2.3)
    design=PB.CableBuilder(
        "parallel-rc-mc",
        PB.Conductor.Solid(:core; radius = 0.010, material = copper),
        PB.Insulator.Tubular(
            :core;
            thickness = Grid(7.0e-3, 5.0),
            material = dielectric
        ),
        PB.Conductor.Tubular(:sheath; thickness = 1.0e-3, material = copper),
        PB.Insulator.Tubular(
            :sheath;
            thickness = 2.0e-3,
            material = outer_dielectric
        )
    )
    problem=PB.SystemBuilder(
        "parallel-rc-mc",
        design,
        PB.at(
            x = 0.0,
            y = -1.0,
            phases = (:core=>1, :sheath=>2)
        );
        length = 1000.0,
        temperature = 20.0,
        earth = PB.Earth(rho = 100.0, eps_r = 10.0),
        frequencies = [50.0, 500.0]
    )
    formulation=Formulation(
        Val(:EMT);
        insulation_admittance = InsulationAdmittance.ParallelRC(),
        modal_transform = nothing,
        options = (
            reduce_bundle = false,
            kron_reduction = false,
            ideal_transposition = true
        )
    )
    policy=MonteCarlo(
        trials = 12,
        distribution = :normal,
        seed = 20260812,
        return_samples = true
    )

    sampled=compute!(problem, formulation; run = policy)
    @test sampled isa MonteCarloResult{<:LineParameters}
    @test size(samples(sampled).G) == (2, 2, 2, 12)
    @test std(samples(sampled).G[1, 1, 1, :]) > 0
    @test std(samples(sampled).C[1, 1, 1, :]) > 0

    reconstructed=Measurements.measurement(sampled)
    @test uncertainty(real(reconstructed.Y[1, 1, 1])) > 0
    @test uncertainty(imag(reconstructed.Y[1, 1, 1])) > 0

    angular=reshape(2π .* frequencies(reconstructed), 1, 1, :)
    joint=vcat(
        vec(real.(reconstructed.Z.values)),
        vec(imag.(reconstructed.Z.values) ./ angular),
        vec(real.(reconstructed.Y.values)),
        vec(imag.(reconstructed.Y.values) ./ angular)
    )
    retained=samples(sampled)
    sample_matrix=vcat(
        reshape(retained.R, :, 12),
        reshape(retained.L, :, 12),
        reshape(retained.G, :, 12),
        reshape(retained.C, :, 12)
    )
    means=vec(mean(sample_matrix; dims = 2))
    centered=sample_matrix .- means
    empirical_covariance=centered*transpose(centered)/11
    @test value.(joint) ≈ means
    @test uncertainty.(joint) ≈ vec(std(sample_matrix; dims = 2))
    @test Measurements.cov(joint) ≈ empirical_covariance

    repeated=compute!(problem, formulation; run = policy)
    for component in (:R, :L, :G, :C)
        @test getproperty(samples(repeated), component) ==
              getproperty(samples(sampled), component)
    end

    total=compute!(
        problem,
        formulation;
        run = policy,
        options = (output_basis = :total,)
    )
    @test basis(LineCableModels.result(sampled)) === :per_length
    @test basis(LineCableModels.result(total)) === :total
    for component in (:R, :L, :G, :C)
        @test getproperty(samples(total), component) ≈
              1000.0 .* getproperty(samples(sampled), component)
    end
end
