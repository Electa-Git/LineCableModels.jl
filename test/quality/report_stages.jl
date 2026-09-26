@testitem "Quality / report stage order / observed protocol calls" tags=[:quality] setup=[TestFixtures] begin
    using Measurements, Statistics
    const RB=LineCableModels.ReportBuilder
    line=TestFixtures.two_conductor_results()
    cable=CableConstants(0.2,3e-6,4e-10,5e-9)
    mc=TestFixtures.cable_monte_carlo_result()
    cases=((RB.CableConstantsTableDefinition(),ObservedResult(cable)),
        (RB.LineParametersTableDefinition((R,L,G,C),:base,:base,:base,false),ObservedResult(line,(R,L,G,C))),
        (RB.MonteCarloTableDefinition(:base,nothing),observables(mc,((LineCableModels.UQ.statistics,R,mean),(LineCableModels.UQ.statistics,R,std)))),
        (XLSXReportDefinition(),ObservedResult(line)),
        (TableReportDefinition((R,)),ObservedResult(cable)))
    for (definition,observed) in cases
        @test applicable(RB.select,definition,observed)
        @test applicable(RB.tabulate,definition,observed)
        @test RB.select(definition,observed) !== nothing
        @test RB.tabulate(definition,observed) !== nothing
        @test which(report,(typeof(definition),typeof(observed))).module === RB
    end
    benchmark=report(RB.BenchmarkTableDefinition(bands=(:all,)),(reference=line,candidate=line))
    @test RB.tabulate(RB.BenchmarkTableDefinition(bands=(:all,)),benchmark.observed;reference=benchmark.reference) !== nothing
    # Selection and tabulation are required; the remaining stages have explicit defaults.
    struct UnimplementedReport <: RB.AbstractReportDefinition end
    @test_throws MethodError report(UnimplementedReport(),ObservedResult(line))
    @test_throws MethodError report(UnimplementedReport(),line)
end
