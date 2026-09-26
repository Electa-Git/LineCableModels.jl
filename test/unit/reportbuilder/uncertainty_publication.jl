@testitem "UQ / publication / selected matrix coordinates retain axes and units" tags=[:unit] begin
    using DataFrames
    using Statistics
    using Measurements
    const UQ = LineCableModels.UQ
    const U = LineCableModels.Units
    frequencies = [50.0, 1000.0]
    positions = reshape(collect(1.0:8.0), 2, 2, 2)
    quantities = (R=positions .* 1e-4, L=positions .* 1e-7,
        C=positions .* 1e-10, G=positions .* 1e-9)
    omega = reshape(2pi .* frequencies, 1, 1, :)
    parameters = LineParameters(quantities.R .+ im .* omega .* quantities.L,
        quantities.G .+ im .* omega .* quantities.C, frequencies)
    summaries = map(quantities) do values
        map(value -> UQ.SampleSummary([0.5value, 1.5value]), values)
    end
    sampled_parameters = LineCableModels.materialize(parameters,summaries)
    source = MonteCarloResult(MonteCarlo(Formulation(); trials=2, seed=7),
        [sampled_parameters], [summaries], nothing, nothing, UInt64(7), UInt64[11], [2])
    RB=LineCableModels.ReportBuilder
    for indices in ((),(:,:,:),(2,:,:),(:,1,:),(:,:,2),(2,1,2),([2,1],[2],[2,1]))
        selected=isempty(indices) ? ([1,2],[1,2],[1,2]) :
            map(index -> index isa Colon ? [1,2] : index isa Integer ? [index] : collect(index),indices)
        observed=ObservedResult(source,1,((statistics,R,indices...),(statistics,C,indices...));
            frequency_unit=:kilo,length_unit=:base,quantity_units=Dict((statistics,C)=>:nano))
        @test length(observed.quantities)==14
        for q in observed.quantities
            @test q.coordinates.rows==selected[1]
            @test q.coordinates.columns==selected[2]
            @test q.coordinates.samples==selected[3]
            @test q.coordinates.frequencies==frequencies[selected[3]]./1000
            @test U.label(q.coordinates.frequency_unit)=="kHz"
            selector=LineCableModels.Grammar.request_identity(q.request)[2]
            transform=LineCableModels.Grammar.request_identity(q.request)[3]
            expected=observe(source,statistics,selector,transform,1,indices...)
            factor=selector===C ? 1e9 : 1.
            @test q.values≈factor.*expected
            table=RB.tabulate(observed,q.request)
            @test size(table)==(length(selected[3]),1+length(selected[1])*length(selected[2]))
            table[1,2]=-1.
            @test first(source.stats[1].R).mean==1e-4
        end
        @test observed.gridpoint.sampling.trials==2
        @test observed.gridpoint.sampling.point_seed==UInt64(11)
    end
    request=@observe (statistics,R,mean)[:,:,:]
    capacitance_request=@observe (statistics,C,mean)[[2,1],[2],[2,1]]
    overridden=ObservedResult(source,1,(capacitance_request,);length_unit=:base,
        quantity_units=Dict(capacitance_request=>:nano,(statistics,C)=>:pico))
    @test only(overridden.quantities).values≈quantities.C[[2,1],[2],[2,1]].*1e9
    explicit=ObservedResult(source,1,(capacitance_request,);length_unit=:base,
        units=(:pico,),quantity_units=Dict(capacitance_request=>:nano))
    @test only(explicit.quantities).values≈quantities.C[[2,1],[2],[2,1]].*1e12
    @test only(explicit.quantities).coordinates==only(overridden.quantities).coordinates
    @test request==(statistics,R,mean,Colon(),Colon(),Colon())
    @test (@observe source (statistics,R,mean)[1,:,:,:])==observe(source,statistics,R,mean,1,:,:,:)
    @test (@observe parameters (Y,abs)[:,:,:])==observe(parameters,Y,abs,:,:,:)
    selected=ObservedResult(source,1,(request,(statistics,R,std));length_unit=:base)
    @test observe(selected,statistics,R,mean)==quantities.R
    @test observe(selected,statistics,R,std)==std.(summaries.R)
    @test_throws ArgumentError observe(selected,statistics,X,mean)
    for (base,derived) in ((L,X),(C,B)), transform in (mean,std)
        @test observe(source,statistics,derived,transform,1)≈omega.*observe(source,statistics,base,transform,1)
        acquired=ObservedResult(source,1,((statistics,derived,transform),);length_unit=:base,quantity_units=:base)
        @test observe(acquired,statistics,derived,transform)≈observe(source,statistics,derived,transform,1)
    end
    q05=Base.Fix2(quantile,.05)
    @test observe(source,statistics,R,q05,1)==getproperty.(summaries.R,:q05)
    @test_throws ArgumentError quantile(first(summaries.R),.1)
    @test only(source)===sampled_parameters

    uncertain_z = measurement.(real.(parameters.Z), 1e-6) .+
        im .* measurement.(imag.(parameters.Z), 2e-6)
    lep_core = LineParameters(uncertain_z, parameters.Y, frequencies)
    lep = LinearErrorResult(LinearError(Formulation()), [lep_core])
    @test only(lep) === lep_core
    @test observe(lep, statistics, R, mean, 1) == quantities.R
    @test observe(lep, statistics, R, std, 1) == fill(1e-6, 2, 2, 2)
    @test observe(lep, statistics, Z, std, 1) ≈ fill(hypot(1e-6, 2e-6), 2, 2, 2)
    lep_observed=ObservedResult(lep,1,(request,(statistics,R,std));length_unit=:base)
    @test observe(lep_observed,statistics,R,mean)==quantities.R
    @test observe(lep_observed,statistics,R,std)==fill(1e-6,2,2,2)
    @test_throws ArgumentError ObservedResult(lep,1,((statistics,R,median),))
    IE=LineCableModels.ImportExport
    restored_mc=IE.deserialize_value(IE.serialize_value(source))
    @test observe(restored_mc,statistics,R,mean,1) == observe(source,statistics,R,mean,1)
    @test NamedTuple(first(only(UQ.statistics(restored_mc)).R)) == NamedTuple(first(summaries.R))
    @test UQ.confidence(restored_mc,1).trials == 2
    restored_lep=IE.deserialize_value(IE.serialize_value(lep))
    @test observe(restored_lep,statistics,R,std,1) == observe(lep,statistics,R,std,1)
    shared=measurement(1.0,0.1)
    correlated=LineParameters(fill(shared+im*shared,2,2,2),fill(2shared+im*shared,2,2,2),frequencies)
    restored=IE.deserialize_value(IE.serialize_value(LinearErrorResult(LinearError(Formulation()),[correlated])))
    @test uncertainty(real(only(restored).Z[1])-real(only(restored).Z[2])) == 0
    @test uncertainty(real(only(restored).Y[1])-2real(only(restored).Z[1])) == 0

    definition=RB.BenchmarkTableDefinition(((statistics,R,mean),(statistics,R,std),(statistics,B,mean));bands=(:all,))
    artifact=RB.report(definition,(reference=source,candidate=source))
    @test length(artifact.tables.features)==3
    @test names(first(artifact.tables.features).relative)==["formula","all"]
    @test size(artifact.tables.sampling,1)==2
    @test all(iszero,skipmissing(artifact.tables.terms.relative_rms_percent))
    for row in eachrow(artifact.tables.mean_sampling_precision)
        selector=getproperty(LineCableModels,row.quantity)
        expected=observe(source,statistics,selector,std,1)
        @test row.standard_error≈expected[row.index...]/sqrt(2)
        @test row.row==row.index[1] && row.column==row.index[2]
        @test row.frequency_Hz==frequencies[row.index[3]]
        @test row.unit==LineCableModels.Units.label(LineCableModels.Units.native_unit(
            LineCableModels.Units.quantity(selector),:pul))
    end
    @test all(==(:empirical),[point.gridpoint.uncertainty.estimator for point in artifact.observed])
end
