@testitem "ReportBuilder / formulation axes, five bands and unavailable relative RMS" tags=[:integration] setup=[TestFixtures, FormulaContractModels] begin
    using LineCableModels.Engine: compare
    using LineCableModels.ReportBuilder
    using DataFrames
    system=TestFixtures.three_phase_system()
    problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),frequencies=[0.1,1.,10.,100.,1e3,1e4,1e5,1e6,1e7])
    space=Formulation(earth_impedance=Grid((formula(:default),
        FormulaContractModels.selection(LineCableModels.Engine.EarthImpedance; layers=2:2))))
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
    definition=BenchmarkTableDefinition(;bands=(:all,:dc,:harmonic,:narrow,:wide))
    @test definition.settings.bands == (:all,:dc,:harmonic,:narrow,:wide)
    artifact=report(definition,(;reference,candidate=candidates))
    @test artifact.illustration === nothing
    @test nrow(artifact.tables.maxima)==60
    # The three shunt quantities share one candidate; numerical records for
    # both completed calculations remain available to the saved-result writer.
    @test length(metadata(artifact.tables.maxima,"comparison_records"))==60
    @test nrow(artifact.tables.comparisons)==60
    @test nrow(artifact.tables.summary)==60
    @test Set(artifact.tables.terms.formulation_index)==Set((1,2))
    @test Set(zip(artifact.tables.terms.row,artifact.tables.terms.column)) ==
        Set((i,j) for i in 1:size(reference.Z,1),j in 1:size(reference.Z,2))
    for q in (Z,Y,R,L,G,C),band in definition.settings.bands,i in 1:2
        exact=compare(reference,candidates[i],q;band)
        row=only(filter(row -> row.quantity==nameof(q) && row.band==band && row.formulation_index==i,eachrow(artifact.tables.comparisons)))
        @test isequal(row.absolute_rms,exact.absolute)
        @test isequal(row.relative_rms_percent,100 .* exact.relative)
    end
    @test nrow(report(definition,artifact.observed).tables.summary)==60
    damaged=deepcopy(candidates)
    reverse!(damaged[1].details.data.coordinates)
    @test_throws ArgumentError report(definition,(;reference,candidate=damaged))
    @test_throws ArgumentError compare(candidates,candidates,Z)
    paired=compare(candidates,candidates,Z;pairing=[(1,1),(2,2)])
    @test all(error -> all(iszero,skipmissing(error.absolute)),paired)
    @test paired.axes === candidates.axes
    labels=artifact.tables.formulations.label[artifact.tables.formulations.role .== :candidate]
    @test labels[1] != labels[2]
    @test occursin("unified",lowercase(labels[1]))
    @test occursin("layer",lowercase(labels[2]))
end
