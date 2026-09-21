@testitem "ObservedResult / primary representations and quantity tables" tags=[:unit] begin
    using DataFrames, LinearAlgebra, Measurements
    using LineCableModels.ReportBuilder: tabulate
    using LineCableModels.Grammar: observation_requests
    line=LineParameters(reshape(ComplexF64.(1:12),2,2,3),fill(3.0+4im,2,2,3),[0.,1.,2.])
    observed=ObservedResult(line)
    @test fieldnames(ObservedResult)==(:gridpoint,:quantities,:errors,:timings)
    @test_throws ArgumentError ObservedResult(line,(R,))
    @test_throws ArgumentError ObservedResult(line,(X,L))
    @test_throws DimensionMismatch ObservedResult(line,((R,1,1,:),(X,2,2,:)))
    @test observation_requests(line,(R,);complete_pairs=true).displayed==((R,:,:,:),)
    @test length(ObservedResult(line,(R,);complete_pairs=true).quantities)==4
    table=tabulate(observed,R)
    @test size(table)==(3,5)
    @test names(table)==["frequency","[1,1]","[1,2]","[2,1]","[2,2]"]
    @test table[1,Symbol("[1,2]")]==3000
    @test table[1,Symbol("[2,1]")]==2000
    line.Z.values .= -1
    @test observe(observed,R)[1]==1000
    @test tabulate(observed,(R,2,1,:))[!,2]==[2000,6000,10000]
    @test observe(observed,R,2,1,:)==[2000,6000,10000]
    for requests in ((R,X),(R,L),((Z,abs),(Z,angle)),(G,B),(G,C),((Y,abs),(Y,angle)))
        @test length(ObservedResult(line,requests).quantities)==4
    end
    dc=LineParameters(fill(1.0+0im,1,1,1),fill(1.0+0im,1,1,1),[0.])
    proxies=ObservedResult(dc,(R,L,G,C))
    @test ismissing(only(observe(proxies,L)))
    @test only(observe(proxies,R))==1000
    scalar=ObservedResult(dc,((R,1,1,1),(L,1,1,1),(G,1,1,1),(C,1,1,1)))
    @test ismissing(observe(scalar,L))
    shared=measurement(0.,.2)
    uncertain=LineParameters(fill(complex(shared,2shared),1,1,2),fill(1.0+2im,1,1,2),[1.,2.])
    polar=ObservedResult(uncertain,((Z,abs),(Z,angle)))
    @test all(ismissing,observe(polar,Z,abs))
    @test all(ismissing,observe(polar,Z,angle))
    @test polar.quantities[1].missing_reason[1]===:undefined_first_order_magnitude
    @test uncertainty(polar.quantities[1].unavailable_components.imaginary[1]-2shared)==0
    @test_throws ArgumentError observe(observed,L)
end

@testitem "ObservedResult / benchmark association survives filtering" tags=[:unit] begin
    using LineCableModels.Engine: compare, retain_gridpoint
    using LineCableModels.Grammar: gridpoint_id, observation_groups
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    base=LineParameters(fill(1.0+2im,2,2,3),fill(3.0+4im,2,2,3),[1.,2.,3.])
    reference=retain_gridpoint(base,gridpoint_id())
    source=gridpoint_id().source_id
    candidates=[retain_gridpoint(LineParameters(base.Z.*factor,base.Y,[1.,2.,3.]),
        gridpoint_id(source_id=source,problem_index=index)) for (index,factor) in enumerate((1.1,1.2))]
    errors=compare(reference,candidates,[R,L];bands=(:all,))
    timings=[(candidate_id=LineCableModels.Grammar.observation_gridpoint(candidate).id,seconds=index)
        for (index,candidate) in enumerate(candidates)]
    observed=observables(reverse(candidates);comparisons=errors,timings)
    @test length(observed)==2
    @test observed[1].gridpoint.id.problem_index==2
    @test observed[1].timings.seconds==2
    @test all(error -> error.candidate_id==observed[1].gridpoint.id,observed[1].errors)
    @test all(error -> error.maxima.absolute.value>0,observed[1].errors)
    ref=ObservedResult(reference)
    artifact=report(BenchmarkTableDefinition(),observed[1:1];reference=ref)
    @test artifact.reference===ref
    @test artifact.observed[1].gridpoint.id.problem_index==2
    @test size(artifact.tables.terms,1)==8
    @test all(==(2),artifact.tables.terms.problem_index)
    @test only(artifact.tables.execution.seconds)==2
    @test only(artifact.tables.execution.candidate_source)==string(observed[1].gridpoint.id.source_id)
    @test only(artifact.tables.execution.candidate_point)==2
    @test only(artifact.tables.execution.reference_source)==string(ref.gridpoint.id.source_id)
    equal_timings=[(candidate_id=point.gridpoint.id,seconds=1.,scope=:compute_call_wall) for point in observed]
    equal_observed=observables(reverse(candidates);comparisons=errors,timings=equal_timings)
    timing_report=report(BenchmarkTableDefinition(),equal_observed;reference=ref)
    @test timing_report.tables.execution.seconds==[1.,1.]
    @test timing_report.tables.execution.candidate_point==[2,1]
    @test timing_report.tables.execution.candidate_formulation==[1,1]
    @test_throws ArgumentError ObservedResult(first(candidates);timings=first(timings[2:2]))
    candidates[2].Z.values .= NaN
    @test all(isfinite,skipmissing(artifact.tables.terms.absolute_rms))
    @test_throws ArgumentError report(BenchmarkTableDefinition((C,)),observed;reference=ref)
    @test length(observation_groups(observed;request=R))==2
end

@testitem "ObservedResult / UQ products and archive-wide dependencies" tags=[:unit] begin
    using Statistics, Measurements
    using LineCableModels.ReportBuilder: tabulate
    using LineCableModels.Grammar: gridpoint_id
    f=[1.,2.,3.]
    base=LineParameters(fill(1.0+2im,2,2,3),fill(3.0+4im,2,2,3),f)
    trials=(R=reshape(collect(1.:60.),2,2,3,5),L=fill(0.001,2,2,3,5),C=fill(1e-9,2,2,3,5),G=fill(1e-6,2,2,3,5))
    summaries=map(array -> [SampleSummary(vec(array[i,j,k,:])) for i in 1:2,j in 1:2,k in 1:3],trials)
    core=LineCableModels.materialize(base,summaries)
    core=LineCableModels.Engine.retain_gridpoint(core,gridpoint_id())
    mc=MonteCarloResult(MonteCarlo(Formulation();trials=5,seed=12,return_samples=true),
        [core],[summaries],[trials],nothing,UInt64(12),UInt64[13],[5])
    product=ObservedResult(mc,1,((statistics,R,mean),(statistics,R,std)))
    @test length(product.quantities)==2
    @test observe(product,statistics,R,mean)≈mean.(summaries.R).*1000
    @test observe(product,statistics,R,std)≈std.(summaries.R).*1000
    @test product.gridpoint.sampling.mean_standard_error.R≈std.(summaries.R)./sqrt(5)
    sampled=ObservedResult(mc,1,((samples,R),);length_unit=:base)
    selected_samples=ObservedResult(sampled,((samples,R,2,1,[3,1],[5,2]),))
    sample_product=only(selected_samples.quantities)
    @test sample_product.values==trials.R[2,1,[3,1],[5,2]]
    @test sample_product.coordinates.samples==[3,1]
    @test sample_product.coordinates.trials==[5,2]
    @test sample_product.coordinates.frequencies==[3.,1.]
    all_trials=ObservedResult(mc,1,((samples,R,2,1,[3,1]),);length_unit=:base)
    @test only(all_trials.quantities).values==trials.R[2,1,[3,1],:]
    retained_trials=ObservedResult(sampled,((samples,R,2,1,[3,1]),))
    @test only(retained_trials.quantities).values==only(all_trials.quantities).values
    @test size(tabulate(selected_samples,only(selected_samples.quantities).request),1)==4
    histogram=ObservedResult(mc,1,((histograms,R,1,2,3,3),))
    @test length(histogram.quantities)==1
    distribution=only(histogram.quantities)
    @test sum(distribution.values.count)==5
    @test sum(distribution.values.probability)≈1
    @test distribution.distribution.qq!==nothing
    converted=only(ObservedResult(histogram;length_unit=:base).quantities)
    @test converted.values.lower≈distribution.values.lower./1000
    @test converted.values.density≈distribution.values.density.*1000
    @test converted.values.probability==distribution.values.probability
    @test converted.values.count==distribution.values.count
    @test converted.distribution.qq.sample≈distribution.distribution.qq.sample./1000
    @test converted.distribution.model_cdf.y==distribution.distribution.model_cdf.y
    @test size(tabulate(histogram,distribution.request),1)==3
    trials.R .= -1
    @test all(>(0),distribution.distribution.qq.sample)
    shared=measurement(1.,.2)
    uncertain=LineParameters(fill(complex(shared,2shared),2,2,3),fill(complex(3shared,4shared),2,2,3),f)
    observed=only(observables([uncertain]))
    @test length(observed.quantities)==4
    @test all(q -> q.statistic===:value,observed.quantities)
    @test !haskey(observed.gridpoint,:sampling)
    artifact=report(TableReportDefinition(),[observed];reference=observed)
    mktempdir() do directory
        for suffix in (".json",".jls")
            path=LineCableModels.save(artifact,joinpath(directory,"shared"*suffix))
            restored=LineCableModels.import_data(:observed,path)
            a=observe(restored.observed[1],R)[1]
            b=observe(restored.reference,R)[1]
            x=observe(restored.reference,X)[1]
            @test uncertainty(a-b)==0
            @test uncertainty(x-2a)==0
            @test uncertainty(a)>0
        end
    end
    setprecision(BigFloat,384) do
        precise=LineParameters(fill(complex(big"1.234567890123456789012345678901",big"2"),1,1,2),
            fill(complex(big"3",big"4"),1,1,2),BigFloat[1,2])
        original=ObservedResult(precise;length_unit=:base)
        encoded=LineCableModels.ImportExport.serialize_value(original)
        restored=setprecision(BigFloat,128) do
            LineCableModels.ImportExport.deserialize_value(encoded)
        end
        @test precision(observe(restored,R)[1])==384
        @test observe(restored,R)[1]==observe(original,R)[1]
    end
end
