@testitem "ObservedResult / independent primary owner reuses tables and plotting" tags=[:integration] begin
    using CairoMakie, DataFrames
    using LineCableModels.Grammar: observation_quantity
    const reads=Ref(0)
    struct IndependentResult
        values::Array{Float64,3}
        frequencies::Vector{Float64}
    end
    LineCableModels.basis(::IndependentResult)=:pul
    LineCableModels.observables(::Type{IndependentResult})=(R,)
    function LineCableModels.Grammar.observation_quantity(result::IndependentResult,request;
            unit=nothing,clip=true,atol=nothing,frequencies=nothing)
        reads[]+=1
        n,m,k=size(result.values)
        target=something(unit,LineCableModels.Units.native_unit(R,:pul))
        native=LineCableModels.Units.native_unit(R,:pul)
        factor=LineCableModels.Units.scale_factor(native,target)
        return (request,quantity=LineCableModels.Units.quantity(R),family=:independent,statistic=:value,
            values=copy(result.values).*factor,unit=target,basis=:pul,
            coordinates=(kind=:matrix,indices=(:,:, :),rows=collect(1:n),columns=collect(1:m),
                samples=collect(1:k),frequencies=copy(result.frequencies),frequency_unit=LineCableModels.Units.units(:base,:hertz),
                extent=(n,m,k),labels=string.(1:n),domain=:PhaseDomain),
            thresholds=nothing,available=trues(n,m,k),engineering_zero=falses(n,m,k),clipped=false,missing_reason=nothing)
    end
    source=IndependentResult(reshape(collect(1.:12.),2,2,3),[1.,10.,100.])
    observed=ObservedResult(source)
    @test reads[]==1
    source.values .= NaN
    artifact=report(TableReportDefinition(),observed)
    @test size(artifact.tables.independent.R)==(3,5)
    @test artifact.tables.independent.R[1,Symbol("[1,2]")]==3000
    rendered=LineCableModels.plot(observed;backend=:cairo,display_plot=false,controls=false)
    @test length(rendered.axes)==4
    @test reads[]==1
    @test all(isfinite,observe(observed,R))
end
