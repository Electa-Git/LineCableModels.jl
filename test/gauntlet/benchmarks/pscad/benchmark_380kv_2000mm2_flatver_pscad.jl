@testitem "PSCAD benchmark / 380 kV 2000 mm² flat vertical" tags=[:gauntlet, :pscad] setup=[GauntletSupport] begin
    using Test
    using DataFrames
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport
    using .GauntletSupport.PSCADBenchmarks

    model=load_case(:cable_380kv_2000mm2_flatver)
    reference_formulation=Formulation(
        :pscad; earth_impedance = :Wedepohl
    )
    candidate_formulation=Formulation(
        earth_impedance = :Pollaczek,
        earth_admittance = :IdealGround,
        insulation_admittance = formula(:Lossless),
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
    benchmark=GauntletCase(
        :benchmark_380kv_2000mm2_flatver_pscad,
        :pscad,
        @__FILE__,
        model,
        reference_formulation,
        candidate_formulation,
        tolerances
    )
    candidate=@inferred compute(model.nominal_problem, candidate_formulation)
    @test size(Z(candidate)) == model.expected_size
    outcome=run_case(
        benchmark;
        options = (reference = (
            output_stem = "380kV_2000_flatv",
            verbosity = (default = 0, PSCAD = 2)
        ),)
    )
    @test outcome.reference isa LineParameters
    @test outcome.candidate isa LineParameters
    @test size(Z(outcome.reference)) == model.expected_size
    @test size(Y(outcome.reference)) == model.expected_size
    @test frequencies(outcome.reference) == model.nominal_problem.frequencies
    if !ISHEADLESS
        display(DataFrame(outcome.comparison))
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
