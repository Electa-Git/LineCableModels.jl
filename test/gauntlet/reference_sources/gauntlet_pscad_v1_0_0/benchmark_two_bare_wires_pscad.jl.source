@testitem "PSCAD benchmark / 2 bare wires" tags=[:gauntlet] setup=[GauntletSupport] begin
    using Test
    using DataFrames
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport
    using .GauntletSupport.PSCADBenchmarks

    artificial_conductor=Material(
        rho = eps(Float64), eps_r = 1.0, mu_r = 1.0, T0 = 20.0, alpha = 0.0
    )
    artificial_insulation=Material(
        rho = 1.97e14, eps_r = 2.3, mu_r = 1.0, T0 = 20.0, alpha = 0.0
    )
    core_radius=0.0425
    insulation_thickness=1.0e-3

    core_parts=(
        Conductor.Solid(
            :core;
            radius = core_radius,
            material = artificial_conductor
        ),
        Insulator.Tubular(
            :core; thickness = insulation_thickness, material = artificial_insulation
        )
    )
    nominal_data=LineCableModels.DataModel.NominalData()
    design=CableBuilder(
        "two_bare_wires",
        core_parts;
        nominal = nominal_data
    )

    earth=Earth(rho = 0.1, eps_r = 1.0, mu_r = 1.0)
    frequencies_value=collect(10.0 .^ range(0, stop = 6, length = 101))
    positions=(
        at(
            x = 0.0,
            y = -1.0,
            phases = (:core=>1,)
        ),
        at(
            x = 1.0,
            y = -1.0,
            phases = (:core=>2,)
        )
    )
    problem=SystemBuilder(
        "benchmark_two_bare_wires_pscad",
        design,
        positions;
        length = 1.0,
        temperature = 20.0,
        earth,
        frequencies = frequencies_value
    )
    reference_problem=LineParametersProblem(
        problem.system;
        temperature = problem.temperature,
        earth_props = problem.earth_props,
        frequencies = problem.frequencies
    )
    reference_formulation=Formulation(
        :pscad;
        earth_impedance = EarthImpedance.Wedepohl()
    )
    formulation=Formulation(
        earth_impedance = EarthImpedance.Pollaczek(),
        earth_admittance = EarthAdmittance.IdealGround(),
        insulation_admittance = InsulationAdmittance.Lossless(),
        options = (
            kron_reduction = false,
            reduce_bundle = false,
            ideal_transposition = false
        )
    )
    tolerances=(
        reference = (
            Z = (absolute = 1.0e-6, relative = 5.0e-2),
            Y = (absolute = 1.0e-9, relative = 5.0e-2)
        ),
        regression = (
            Z = (absolute = 1.0e-9, relative = 1.0e-6),
            Y = (absolute = 1.0e-12, relative = 1.0e-6)
        ),
        performance = (
            median_time_ratio = 1.20,
            bytes_ratio = 1.05,
            allocations_ratio = 1.05
        )
    )
    case=GauntletCase(
        :benchmark_two_bare_wires_pscad,
        :pscad,
        @__FILE__,
        reference_problem,
        reference_formulation,
        problem,
        formulation,
        ["cable:1:core", "cable:2:core"],
        (2, 2, length(frequencies_value)),
        tolerances
    )

    if !ISHEADLESS
        display(DataFrame(design, :components))
        display(DataFrame(problem.system))
        display(DataFrame(problem.earth_props))
    end

    inferred=@inferred compute(problem, formulation)
    @test size(Z(inferred)) == (2, 2, 101)
    execution_options=(reference = (
        output_stem = "two_bare_wires",
        verbosity = (default = 0, PSCAD = 2)
    ),)
    outcome=run_case(case; options = execution_options)
    @test outcome.reference isa LineParameters
    @test outcome.candidate isa LineParameters
    @test size(Z(outcome.reference)) == (2, 2, 101)
    @test size(Y(outcome.reference)) == (2, 2, 101)
    @test frequencies(outcome.reference) == frequencies_value

    if !ISHEADLESS
        errors=DataFrame(
            outcome.comparison;
            zero_atol = (
                Z = tolerances.reference.Z.absolute,
                Y = tolerances.reference.Y.absolute
            )
        )
        sort!(errors, [:quantity, :rms_absolute], rev = [false, true])
        display(errors)
    end

    if outcome.mode===:snapshot
        @test comparison_passes(outcome.regression.Z, tolerances.regression.Z)
        @test comparison_passes(outcome.regression.Y, tolerances.regression.Y)
        if outcome.performance.comparable
            @test outcome.performance.passes
        end
    else
        @test outcome.reference_execution.backend === :pscad
        @test outcome.reference_execution.exit_code == 0
        @test isfile(outcome.reference_execution.console_path)
        @test outcome.timings.reference >= 0
        @test outcome.timings.julia.samples == 10
    end
end
