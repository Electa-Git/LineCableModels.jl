@testitem "ReportBuilder / explicit retained selections and concrete empty-band errors" tags=[:integration] begin
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    using LineCableModels.Engine: compare
    using DataFrames
    records=[NamedTuple(Formulation()),NamedTuple(Formulation(earth_admittance=:default))]
    samples=Float32[1,10]
    z=fill(ComplexF32(1,2),2,2,2)
    y=fill(ComplexF32(0,1e-6),2,2,2)
    reference=LineParameters(PhaseDomain,z,y,samples;details=(coordinates=["a","b"],))
    candidate=LineParameters(PhaseDomain,2z,2y,samples;details=(coordinates=["a","b"],))
    artifact=report(BenchmarkTableDefinition(quantities=(Z,G),bands=(:all,)),(;reference,candidate))
    @test nrow(report(BenchmarkTableDefinition(quantities=(G,)),artifact.published).table.maxima)==1
    @test_throws r"reanalysis" report(BenchmarkTableDefinition(bands=(:wide,)),artifact.published)
    @test_throws r"reanalysis" report(BenchmarkTableDefinition(atol=(G=1f-8,)),artifact.published)
    @test_throws r"unknown" BenchmarkTableDefinition(no_such_control=true)
    @test_throws r"pairing" BenchmarkTableDefinition(pairing=[(1,1),(2,1)])
    @test_throws r"scalar reference" report(BenchmarkTableDefinition(pairing=[(2,1)]),(;reference,candidate))
    @test all(ismissing,filter(row -> row.quantity === :G,artifact.table.terms).relative_rms_percent)
    @test all(iszero,filter(row -> row.quantity === :G,artifact.table.terms).absolute_rms)
    far=LineParameters(PhaseDomain,z,y,Float32[1e7,2e7];details=(coordinates=["a","b"],))
    space=ParametricResult(nothing,[reference,far],(problems=[:near,:far],formulations=records[1:1]),(;))
    errors=compare(space,space,Z;pairing=[(1,1),(2,2)],band=:wide)
    @test isconcretetype(eltype(errors))
    @test typeof(errors[1]) === typeof(errors[2])
    @test Base.nonmissingtype(eltype(errors[2].absolute)) === Float32
    @test all(ismissing,errors[1].absolute)
    @test all(iszero,errors[2].absolute)
    @test errors.axes === space.axes
    @test errors[2].details.actual_bounds == (1f7,2f7)
end
