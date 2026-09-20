@testitem "ReportBuilder / explicit retained selections and concrete empty-band errors" tags=[:integration] begin
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    using LineCableModels.Engine: compare
    using DataFrames
    records=[NamedTuple(Formulation()),NamedTuple(Formulation(earth_admittance=:default))]
    samples=Float32[1,10]
    z=fill(ComplexF32(1,2),2,2,2)
    y=fill(ComplexF32(0,1e-6),2,2,2)
    reference=LineParameters(PhaseDomain,z,y,samples;details=ComputationDetails(;coordinates=["a","b"],))
    candidate=LineParameters(PhaseDomain,2z,2y,samples;details=ComputationDetails(;coordinates=["a","b"],))
    artifact=report(BenchmarkTableDefinition(quantities=(Z,G),bands=(:all,)),(;reference,candidate))
    @test nrow(report(BenchmarkTableDefinition(quantities=(G,)),artifact.observed).tables.maxima)==1
    @test_throws ArgumentError report(BenchmarkTableDefinition(bands=(:wide,)),artifact.observed)
    @test_throws ArgumentError report(BenchmarkTableDefinition(atol=(G=1f-8,)),artifact.observed)
    @test_throws r"unknown" BenchmarkTableDefinition(no_such_control=true)
    @test_throws r"pairing" BenchmarkTableDefinition(pairing=[(1,1),(2,1)])
    @test_throws r"scalar reference" report(BenchmarkTableDefinition(pairing=[(2,1)]),(;reference,candidate))
    @test all(ismissing,filter(row -> row.quantity === :G,artifact.tables.terms).relative_rms_percent)
    @test all(ismissing,filter(row -> row.quantity === :G,artifact.tables.terms).absolute_rms)
    far=LineParameters(PhaseDomain,z,y,Float32[1e7,2e7];details=ComputationDetails(;coordinates=["a","b"],))
    space=ParametricResult(nothing,[reference,far],(problems=[:near,:far],formulations=records[1:1]), ComputationDetails((;)))
    errors=compare(space,space,Z;pairing=[(1,1),(2,2)],band=:wide)
    @test isconcretetype(eltype(errors))
    @test typeof(errors[1]) === typeof(errors[2])
    @test Base.nonmissingtype(eltype(errors[2].absolute)) === Float32
    @test all(ismissing,errors[1].absolute)
    @test all(iszero,errors[2].absolute)
    @test errors.axes === space.axes
    @test errors[2].details.data.actual_bounds == (1f7,2f7)
end

@testitem "ReportBuilder / two-sided RMS eligibility survives tables and explicit reanalysis" tags=[:integration] begin
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    using LineCableModels.Engine: RMSError
    using DataFrames
    f = [50.0, 500.0]
    z = ones(ComplexF64, 1, 1, 2)
    reference = LineParameters(PhaseDomain,z,fill(1e-6+1e-4im,1,1,2),f;
        details=ComputationDetails(;coordinates=["a"],))
    candidate = LineParameters(PhaseDomain,z,fill(1e-14+1e-4im,1,1,2),f;
        details=ComputationDetails(;coordinates=["a"],))
    definition = BenchmarkTableDefinition(quantities=(G,Y),bands=(:all,))
    artifact = report(definition,(;reference,candidate))
    g = only(eachrow(filter(row -> row.quantity === :G,artifact.tables.terms)))
    @test ismissing(g.relative_rms_percent)
    @test g.status === :candidate_below_tolerance
    @test ismissing(g.absolute_rms)
    @test only(filter(row -> row.quantity === :G,artifact.tables.maxima).unavailable) == 1
    @test !ismissing(only(filter(row -> row.quantity === :Y,artifact.tables.terms).relative_rms_percent))
    # Reading recorded errors never reruns comparison. Explicit fresh comparison
    # operates on the raw operands and leaves the earlier observation unchanged.
    original=first(artifact.observed.errors)
    recorded=merge(original,(relative=fill(1.,1,1),))
    retained=ObservedResult(artifact.observed.gridpoint,artifact.observed.quantities,[recorded],artifact.observed.timings)
    @test only(report(BenchmarkTableDefinition(),retained).tables.terms.relative_rms_percent)==100
    refreshed=report(definition,(;reference,candidate))
    @test ismissing(only(filter(row -> row.quantity===:G,refreshed.tables.terms).relative_rms_percent))
    @test only(first(retained.errors).relative)==1
end
