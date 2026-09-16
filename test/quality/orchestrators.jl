@testitem "Quality / report orchestration / concrete protocol calls" tags=[:quality] setup=[TestFixtures] begin
    using RequiredInterfaces, Measurements
    const RB=LineCableModels.ReportBuilder
    @test RequiredInterfaces.isInterface(RB.AbstractReportDefinition)
    @test Set(RequiredInterfaces.functions(RequiredInterfaces.getInterface(RB.AbstractReportDefinition))) ==
        Set((RB.select,RB.tabulate))
    line=TestFixtures.two_conductor_results()
    cable=CableConstants(0.2,3e-6,4e-10,5e-9)
    mc=TestFixtures.cable_monte_carlo_result()
    cases=((RB.CableConstantsTableDefinition(),cable),
        (RB.LineParametersTableDefinition((@observe(R[:,:,:]),@observe(L[:,:,:]),@observe(G[:,:,:]),@observe(C[:,:,:])),:base,:base,:base,false),line),
        (RB.BenchmarkTableDefinition(),(reference=line,candidate=line)),
        (RB.MonteCarloTableDefinition(:base,nothing),mc),
        (XLSXReportDefinition(),line), (TableReportDefinition((R,)),cable))
    for (definition,source) in cases
        @test applicable(RB.select,definition,source)
        published=RB.select(definition,source)
        @test applicable(RB.tabulate,definition,source,published)
        table=RB.tabulate(definition,source,published)
        @test table !== nothing
        @test which(report,(typeof(definition),typeof(source))).module === RB
    end
    # Mentioning an implementor in the wrong argument is no implementation of
    # this protocol. Execute the required position to expose the fallback.
    struct WrongPositionReport <: RB.AbstractReportDefinition end
    RB.select(source::Integer,::WrongPositionReport)=source
    RB.tabulate(source::Integer,published,::WrongPositionReport)=published
    @test_throws RequiredInterfaces.NotImplementedError RB.select(WrongPositionReport(),23)
    @test_throws RequiredInterfaces.NotImplementedError RB.tabulate(WrongPositionReport(),23,29)
end
