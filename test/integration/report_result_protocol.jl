@testitem "ReportBuilder / an independent result uses observation and coordinate grammar" tags=[:integration] begin
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    const observations=Ref(0)
    struct IndependentResult <: LineCableModels.AbstractCoreResult
        values::Array{ComplexF64,3}
        f::Vector{Float64}
        coordinates::Vector{String}
    end
    LineCableModels.frequencies(result::IndependentResult)=result.f
    LineCableModels.basis(::IndependentResult)=:pul
    LineCableModels.domain(::IndependentResult)=PhaseDomain
    LineCableModels.details(result::IndependentResult)=(coordinates=result.coordinates,)
    function LineCableModels.observe(result::IndependentResult,::typeof(Z))
        observations[]+=1
        return result.values
    end
    z=reshape(complex.(1.:12.,21.:32.),2,2,3)
    reference=IndependentResult(z,[1.,10.,100.],["a","b"])
    candidate=IndependentResult(2z,[1.,10.,100.],["a","b"])
    definition=BenchmarkTableDefinition(quantities=(Z,),bands=(:all,))
    artifact=report(definition,(;reference,candidate))
    @test all(==(1),only(artifact.published.comparisons).error.relative)
    @test length(artifact.table.terms.row)==4
    before=observations[]
    repeated=report(definition,artifact.published)
    # Retained selection observes each operand once to check dimensions. Another
    # RMS calculation would observe both operands a second time.
    @test observations[]-before == 2
    @test repeated.published.comparisons === artifact.published.comparisons ||
        only(repeated.published.comparisons).error === only(artifact.published.comparisons).error
    @test all(==(100),repeated.table.terms.relative_rms_percent)
end
