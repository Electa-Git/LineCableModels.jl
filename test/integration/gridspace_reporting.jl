@testitem "ReportBuilder / formulation axes, five bands and unavailable relative RMS" tags=[:integration] setup=[TestFixtures] begin
    using LineCableModels.Engine: compare
    using LineCableModels.ReportBuilder
    using DataFrames
    system=TestFixtures.three_phase_system()
    problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),frequencies=[0.1,1.,10.,100.,1e3,1e4,1e5,1e6,1e7])
    space=Formulation(earth_impedance=Grid((:default,:Pollaczek1926)))
    reference=compute(problem,Formulation())
    candidates=compute(problem,space)
    errors=compare(reference,candidates,Z)
    @test errors.axes === candidates.axes
    @test length(errors)==2
    for i in 1:2
        scalar=compare(reference,candidates[i],Z)
        @test isequal(errors[i].absolute,scalar.absolute)
        @test isequal(errors[i].relative,scalar.relative)
    end
    definition=BenchmarkTableDefinition()
    @test definition.settings.bands == (:all,:dc,:harmonic,:narrow,:wide)
    artifact=report(definition,(;reference,candidate=candidates))
    @test artifact.illustration === nothing
    @test nrow(artifact.table.maxima)==60
    @test nrow(artifact.table.summary)==10
    @test Set(artifact.table.terms.formulation_index)==Set((1,2))
    @test Set(zip(artifact.table.terms.row,artifact.table.terms.column)) ==
        Set((i,j) for i in 1:size(reference.Z,1),j in 1:size(reference.Z,2))
    for q in (Z,Y,R,L,G,C),band in definition.settings.bands,i in 1:2
        exact=compare(reference,candidates[i],q;band)
        row=only(filter(row -> row.quantity==nameof(q) && row.band==band && row.formulation_index==i,eachrow(artifact.table.comparisons)))
        @test isequal(row.absolute_rms,exact.absolute)
        @test isequal(row.relative_rms_percent,100 .* exact.relative)
    end
    @test nrow(report(definition,artifact.published).table.summary)==10
    damaged=deepcopy(candidates)
    reverse!(damaged[1].details.coordinates)
    @test_throws ArgumentError report(definition,(;reference,candidate=damaged))
    @test_throws ArgumentError compare(candidates,candidates,Z)
    paired=compare(candidates,candidates,Z;pairing=[(1,1),(2,2)])
    @test all(error -> all(iszero,error.absolute),paired)
    @test paired.axes === candidates.axes
    labels=artifact.table.formulations.label[artifact.table.formulations.role .== :candidate]
    @test labels[1] != labels[2]
    @test occursin("default",labels[1])
    @test occursin("Pollaczek1926",labels[2])
end
