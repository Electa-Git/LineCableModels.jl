@testitem "Engine / dielectric admittance / analytical layer and lossless limit" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using TOML

    formulation=InsulationAdmittance.Formula(:Ametani2004)
    lossless=InsulationAdmittance.Formula(:Gustavsen2013)
    semicon=SemiconAdmittance.Formula(:Ametani2004)

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
            material=Material(
                :insulator, resistivity, relative_permittivity, one(T),
                T(20), zero(T)
            )
            evaluated=formulation(material, T(50), T(20))
            coefficient=LineCableModels.Engine.potential_coefficient(
                r_inner, r_outer, evaluated, frequency_point
            )
            ε0=one(T)*T(88541878128)*(one(T)*T(10))^(-22)
            circumference=T(2)*(one(T)*T(π))
            log_ratio=log(r_outer/r_inner)
            capacitance=circumference*ε0*relative_permittivity/log_ratio
            conductance=circumference*inv(resistivity)/log_ratio
            previous=frequency_point/(conductance+frequency_point*capacitance)
            T===Float64&&@test coefficient == previous
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
    material=Material(:insulator, rho, eps_r, 1.0, 20.0, 0.0)
    coefficient=LineCableModels.Engine.potential_coefficient(
        r_in, r_ex, formulation(material, 50.0, 20.0), s)

    @test s / coefficient ≈ conductance + s * capacitance
    @test real(s / coefficient) > 0
    @test imag(s / coefficient) > 0

    r_mid=0.013
    second_rho=8.0e10
    second_eps=3.6
    first_material=Material(:semicon, rho, eps_r, 1.0, 20.0, 0.0)
    second_material=Material(:insulator, second_rho, second_eps, 1.0, 20.0, 0.0)
    first_layer=LineCableModels.Engine.potential_coefficient(
        r_in, r_mid, semicon(first_material, 50.0, 20.0), s)
    second_layer=LineCableModels.Engine.potential_coefficient(
        r_mid, r_ex, formulation(second_material, 50.0, 20.0), s)
    first_admittance=s/first_layer
    second_admittance=s/second_layer
    series_admittance=inv(inv(first_admittance)+inv(second_admittance))
    @test s / (first_layer + second_layer) ≈ series_admittance

    lossless_κ=lossless(material, 50.0, 20.0)
    @test iszero(real(lossless_κ))
    lossless_coefficient=LineCableModels.Engine.potential_coefficient(
        r_in, r_ex, lossless_κ, s
    )
    formula_epsilon0=one(r_in)*88541878128*(one(r_in)*10)^(-22)
    @test lossless_coefficient ≈ Complex(
        log(r_ex/r_in)/(2*(one(r_in)*π)*formula_epsilon0*eps_r)
    )
end

@testitem "Engine / dielectric admittance / strict layers reach direct computation" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using LinearAlgebra

    function two_terminal_problem(;
            uncertain = false,
            frequencies = [1.0e-3, 50.0, 1.0e6]
    )
        copper=Material(:conductor, 1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
        rho_1=uncertain ? measurement(2.0e11, 2.0e10) : 2.0e11
        eps_1=uncertain ? measurement(2.4, 0.12) : 2.4
        thickness_1=uncertain ? measurement(3.0e-3, 1.5e-4) : 3.0e-3
        dielectric_1=Material(:semicon, rho_1, eps_1, 1.0, 20.0, 0.0)
        dielectric_2=Material(
            :insulator, 8.0e10, 3.6, 1.0, 20.0, 0.0;
            tan_delta = 0.012
        )
        outer_dielectric=Material(:insulator, 1.0e14, 2.3, 1.0, 20.0, 0.0)

        core_outer=0.010
        first_outer=core_outer+thickness_1
        insulation_outer=first_outer+4.0e-3
        sheath_outer=insulation_outer+1.0e-3
        jacket_outer=sheath_outer+2.0e-3
        design=build(
            CableDesign,
            "parallel-rc-test",
            Stack(AbstractCablePart[
                Group(:core, Region(:core_conductor, Annulus(0.0, core_outer), copper)),
                Region(:insulation_1, Annulus(core_outer, first_outer), dielectric_1),
                Region(:insulation_2, Annulus(first_outer, insulation_outer), dielectric_2),
                Group(
                    :sheath,
                    Region(
                        :sheath_conductor,
                        Annulus(insulation_outer, sheath_outer),
                        copper
                    )),
                Region(
                    :outer_dielectric,
                    Annulus(sheath_outer, jacket_outer),
                    outer_dielectric
                )
            ])
        )
        system=build(
            LineCableSystem,
            design,
            Pose2(0.0, -1.0, 0.0);
            connections = Dict(:core=>1, :sheath=>2),
            system_id = "parallel-rc-test",
            line_length = 1000.0
        )
        earth=EarthModel(100.0, 10.0, 1.0)
        return LineParametersProblem(
            system;
            temperature = 20.0,
            earth_props = earth,
            frequencies
        )
    end

    formulation=Formulation(;
        insulation_admittance = formula(:Gustavsen2013),
        semicon_admittance = formula(:Ametani2004),
        earth_admittance = :IdealGround,
        options = (
            reduce_bundle = false,
            kron_reduction = false,
            ideal_transposition = true
        )
    )
    problem=two_terminal_problem()
    execution=computation_options(Val(LineCableModelsCoaxial), (;))
    workspace=LineParametersWorkspace(
        LineCableModelsCoaxial(), problem, formulation, execution)
    input=workspace.input
    @test !hasproperty(input, :rho_ins)
    @test !hasproperty(input, :eps_ins)
    @test !hasproperty(input, :tan_ins)
    @test getproperty.(input.dielectric_materials, :kind) ==
          [:semicon, :insulator, :insulator]
    @test input.semicon_layer_indices == [1]
    @test input.insulation_layer_indices == [2, 3]
    @test input.dielectric_materials[2].tan_delta == 0.012
    @test @inferred(LineCableModels.Engine.dielectric!(
        workspace.buffers.layer_coefficients,
        workspace,
        1,
        formulation
    )) === workspace.buffers.layer_coefficients
    parameters=compute(problem, formulation; options = (trace = true,))
    trace=details(parameters).trace
    public_parameters=compute(problem, formulation)
    @test public_parameters.Z.values == parameters.Z.values
    @test public_parameters.Y.values == parameters.Y.values
    @test input.insulator_layer_ranges == [1:2, 3:3]
    @test input.r_ins_layer_in[2] ≈ input.r_ins_layer_ext[1]
    @test size(parameters.Y) == (2, 2, 3)

    for frequency_index in eachindex(input.freq)
        s=input.jω[frequency_index]
        function layer_coefficient(layer_index)
            material=input.dielectric_materials[layer_index]
            relation=material.kind===:semicon ?
                     formulation.methods.semicon_admittance :
                     formulation.methods.insulation_admittance
            evaluated=constitutive(
                relation, material, input.freq[frequency_index], input.temperature)
            return LineCableModels.Engine.potential_coefficient(
                input.r_ins_layer_in[layer_index],
                input.r_ins_layer_ext[layer_index],
                evaluated,
                s
            )
        end
        inner=sum(layer_coefficient, input.insulator_layer_ranges[1])
        outer=sum(layer_coefficient, input.insulator_layer_ranges[2])
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
        coefficient=zero(s)
        for layer_index in input.insulator_layer_ranges[component_index]
            material=input.dielectric_materials[layer_index]
            relation=material.kind===:semicon ?
                     formulation.methods.semicon_admittance :
                     formulation.methods.insulation_admittance
            evaluated=constitutive(
                relation, material, input.freq[frequency_index], input.temperature)
            coefficient+=LineCableModels.Engine.potential_coefficient(
                input.r_ins_layer_in[layer_index],
                input.r_ins_layer_ext[layer_index],
                evaluated,
                s
            )
        end
        return coefficient
    end
    near_dc=input.jω[1]/component_coefficient(1, 1)
    high_frequency=input.jω[end]/component_coefficient(1, length(input.jω))
    @test real(near_dc) > 0
    @test imag(near_dc) > 0
    @test imag(high_frequency) > real(high_frequency)

    uncertain_parameters=compute(two_terminal_problem(uncertain = true), formulation)
    @test uncertainty(real(uncertain_parameters.Y[1, 1, 2])) > 0
    @test uncertainty(imag(uncertain_parameters.Y[1, 1, 2])) > 0

    shifted_reference=compute(
        two_terminal_problem(frequencies = [1.0e-3, 50.0]),
        formulation
    )
    direct_reference=compute(
        two_terminal_problem(frequencies = [50.0]),
        formulation
    )
    @test shifted_reference.Y.values[:, :, 2] ==
          direct_reference.Y.values[:, :, 1]

    total=compute(
        problem,
        formulation;
        options = (output_basis = :total,)
    )
    @test basis(parameters) === :pul
    @test basis(total) === :total
    @test total.Z.values ≈ problem.system.line_length .* parameters.Z.values
    @test total.Y.values ≈ problem.system.line_length .* parameters.Y.values
end

@testitem "Engine / Ametani2004 / Gridspace Monte Carlo samples before assembly" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using Statistics
    import LineCableModels.ParametricBuilder as PB

    copper=PB.Material(
        kind = :conductor,
        rho = 1.7241e-8,
        eps_r = 1.0,
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.00393
    )
    dielectric=PB.Material(
        kind = :insulator,
        rho = Grid(2.0e11, 10.0),
        eps_r = Grid(2.4, 5.0)
    )
    outer_dielectric=PB.Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)
    thickness=Grid(7.0e-3, 5.0)
    build_design=function (resolved_dielectric, resolved_thickness)
        core_outer=0.010
        insulation_outer=core_outer+resolved_thickness
        sheath_outer=insulation_outer+1.0e-3
        jacket_outer=sheath_outer+2.0e-3
        return build(
            CableDesign,
            "parallel-rc-mc",
            Stack(AbstractCablePart[
                Group(:core, Region(:core_conductor, Annulus(0.0, core_outer), copper)),
                Region(
                    :core_insulation,
                    Annulus(core_outer, insulation_outer),
                    resolved_dielectric
                ),
                Group(
                    :sheath,
                    Region(
                        :sheath_conductor,
                        Annulus(insulation_outer, sheath_outer),
                        copper
                    )),
                Region(
                    :outer_insulation,
                    Annulus(sheath_outer, jacket_outer),
                    outer_dielectric
                )
            ])
        )
    end
    design=Gridspace{CableDesign}(build_design, (dielectric, thickness))
    system=build(
        LineCableSystem,
        design,
        Pose2(0.0, -1.0, 0.0);
        connections = Dict(:core=>1, :sheath=>2),
        system_id = "parallel-rc-mc",
        line_length = 1000.0
    )
    problem=LineParametersProblem(
        system;
        temperature = 20.0,
        earth_props = PB.Earth(rho = 100.0, eps_r = 10.0),
        frequencies = [50.0, 500.0]
    )
    formulation=Formulation(;
        insulation_admittance = formula(:Ametani2004),
        options = (
            reduce_bundle = false,
            kron_reduction = false,
            ideal_transposition = true
        )
    )
    policy=MonteCarlo(
        formulation;
        trials = 12,
        distribution = :normal,
        seed = 20260812,
        return_samples = true
    )

    parametric_problem=ParametricProblem(problem)
    sampled=compute(parametric_problem, policy)
    @test sampled isa MonteCarloResult{<:LineParameters}
    @test size(only(samples(sampled)).G) == (2, 2, 2, 12)
    @test std(only(samples(sampled)).G[1, 1, 1, :]) > 0
    @test std(only(samples(sampled)).C[1, 1, 1, :]) > 0

    retained=samples(sampled)
    @test !applicable(Measurements.measurement, sampled)
    sample_matrix=vcat(
        reshape(only(retained).R, :, 12),
        reshape(only(retained).L, :, 12),
        reshape(only(retained).G, :, 12),
        reshape(only(retained).C, :, 12)
    )
    means=vec(mean(sample_matrix; dims = 2))
    centered=sample_matrix .- means
    empirical_covariance=centered*transpose(centered)/11
    @test all(isfinite, empirical_covariance)
    @test any(!iszero, empirical_covariance)

    repeated=compute(parametric_problem, policy)
    for component in (:R, :L, :G, :C)
        @test getproperty(only(samples(repeated)), component) ==
              getproperty(only(samples(sampled)), component)
    end

    total=compute(
        ParametricProblem(problem, (output_basis = :total,)),
        policy
    )
    @test basis(only(LineCableModels.result(sampled))) === :pul
    @test basis(only(LineCableModels.result(total))) === :total
    for component in (:R, :L, :G, :C)
        @test getproperty(only(samples(total)), component) ≈
              1000.0 .* getproperty(only(samples(sampled)), component)
    end
end
