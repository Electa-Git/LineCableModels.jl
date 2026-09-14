@testitem "Gauntlet / product and zip result spaces preserve every matrix and axis" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    model=load_case(:two_insulated_wires; variation = ExactOverrides(frequencies = [
        1.0, 37.0]))
    physical=(reduce_bundle = false, kron_reduction = false, ideal_transposition = false)
    reference=BenchmarkCalculation(:reference, model.problem, Formulation(;
        options = physical))
    for combine in (:product, :zip)
        space=Formulation(earth_impedance = Grid((:default, :Pollaczek1926)),
            insulation_admittance = Grid((:default, :Ametani2004)); options = physical, combine)
        observed=Int[]
        callback=(problem, index, result)->push!(observed, index)
        candidate=BenchmarkCalculation(:grid, model.problem, space; options = (on_result = callback,))
        declaration=benchmark_definition(
            :grid_comparison, model.id, :fixture, @__FILE__, model,
            reference, candidate, (; quantities = (:Z, :Y, :G)), (;))
        mktempdir() do directory
            value=run_benchmark(declaration; directory)
            expected=combine===:product ? 4 : 2
            @test length(value.candidate_result)==expected
            @test length(value.comparison)==15expected
            @test observed==collect(1:expected)
            @test length(value.candidate_result.axes.formulations)==expected
            for (index, formulation) in enumerate(space)
                scalar=compute(model.problem, formulation)
                @test value.candidate_result[index].Z == scalar.Z
                @test value.candidate_result[index].Y == scalar.Y
            end
            saved=read_calculation(joinpath(directory, "candidate", "calculation.jld2"))
            @test length(saved.result)==expected
            @test saved.result.axes.problems[1].terminal_order==model.problem.system.terminal_order
            @test all(saved.result[index].Z == value.candidate_result[index].Z
            for index in 1:expected)
            @test length(only(read_benchmark(directory).analyses)["reference_comparison"])==15expected
            retained=read_benchmark(directory)
            tables=report(BenchmarkTableDefinition(false),retained).table
            @test length(tables.comparisons.quantity)==15expected
            @test Set(tables.comparisons.candidate_point)==Set(1:expected)
            @test tables.calculations.axes[2]==saved.result.axes
            @test report(BenchmarkTableDefinition(),retained).illustration === nothing
            @test_throws r"Plotting is optional" LineCableModels.plot(retained, (R,))
            previous=length(observed)
            @test run_benchmark(declaration; directory).timings.execution.candidate.reused
            @test length(observed)==previous
        end
    end
end
