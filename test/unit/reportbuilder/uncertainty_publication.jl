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
    source = MonteCarloResult(MonteCarlo(Formulation(); trials=2, seed=7),
        [parameters], [summaries], nothing, nothing, UInt64(7), UInt64[11], [2])
    statistics_order = (:mean, :std, :min, :q05, :median, :q95, :max)

    for indices in ((), (Colon(), Colon(), Colon()), (2, Colon(), Colon()),
            (Colon(), 1, Colon()), (Colon(), Colon(), 2),
            (2, 1, 2), ([2, 1], [2], [2, 1]))
        selected = isempty(indices) ? ([1, 2], [1, 2], [1, 2]) :
            map(index -> index isa Colon ? [1, 2] :
                index isa Integer ? [index] : collect(index), indices)
        publication = observables(source,
            ((statistics, R, 1, indices...), (statistics, C, 1, indices...));
            frequency_unit=:kilo, length_unit=:base,
            quantity_units=Dict((statistics, C)=>:nano), clip=false)
        frame = DataFrame(publication)
        expected_keys = [(1, frequencies[f] / 1000, row, column, statistic)
            for f in selected[3] for row in selected[1] for column in selected[2]
            for statistic in statistics_order]
        @test collect(zip(frame.point, frame.frequency, frame.row, frame.column,
            frame.statistic)) == expected_keys
        @test frame.R ≈ [getproperty(summaries.R[row, column, f], statistic)
            for f in selected[3] for row in selected[1] for column in selected[2]
            for statistic in statistics_order]
        @test frame.C ≈ [1e9 * getproperty(summaries.C[row, column, f], statistic)
            for f in selected[3] for row in selected[1] for column in selected[2]
            for statistic in statistics_order]
        @test all(==(2), frame.trials)
        @test all(==(UInt64(11)), frame.point_seed)
        @test publication.metadata.row_order == (:point, :frequency, :row, :column, :statistic)
        @test U.label(publication.metadata.observation_columns.frequency.unit) == "kHz"
        @test U.label(publication.metadata.observation_columns.R.unit) == "Ω/m"
        @test U.label(publication.metadata.observation_columns.C.unit) == "nF/m"
        frame.R[1] = -1.0
        @test all(>(0), quantities.R)
        @test source.stats[1].R[1].mean == 1e-4
    end
    @test_throws DimensionMismatch observables(source,
        ((statistics, R, 1, [1], [1], [1]), (statistics, C, 1, [2], [1], [1])))
    @test_throws ArgumentError observables(source,
        ((statistics, R, 1), (statistics, R, 1)))

    request = @observe (statistics, R, mean)[1, :, :, :]
    @test (@observe source (statistics,R,mean)[1,:,:,:])==observe(source,statistics,R,mean,1,:,:,:)
    @test (@observe parameters R[:,:,:])==observe(parameters,R,:,:,:)
    @test (@observe parameters (Y,abs)[:,:,:])==observe(parameters,Y,abs,:,:,:)
    @test request == (statistics, R, mean, 1, Colon(), Colon(), Colon())
    @test LineCableModels.Grammar.request_identity(request) == (statistics, R, mean)
    @test LineCableModels.Grammar.request_indices(request) == (1, Colon(), Colon(), Colon())
    selected = observables(source,
        (request, (statistics, R, std, 1)); length_unit=:base, clip=false)
    frame = DataFrame(selected)
    @test unique(frame.statistic) == [:mean, :std]
    @test frame.R[1:2] ≈ [mean(summaries.R[1]), std(summaries.R[1])]
    @test length(selected.metadata.observation_columns.R.resolution) == 2
    @test observe(source, statistics, B, mean, 1) ≈ quantities.C .* omega
    @test observe(source, statistics, Z, std, 1) ≈
        hypot.(std.(summaries.R), omega .* std.(summaries.L))
    q05 = Base.Fix2(quantile, 0.05)
    @test observe(source, statistics, R, q05, 1) == getproperty.(summaries.R, :q05)
    @test DataFrame(observables(source, ((statistics, R, q05, 1),);
        length_unit=:base)).statistic == fill(:q05, 8)
    @test_throws ArgumentError quantile(first(summaries.R), 0.1)
    @test only(source) === parameters

    uncertain_z = measurement.(real.(parameters.Z), 1e-6) .+
        im .* measurement.(imag.(parameters.Z), 2e-6)
    lep_core = LineParameters(uncertain_z, parameters.Y, frequencies)
    lep = LinearErrorResult(LinearError(Formulation()), [lep_core])
    @test only(lep) === lep_core
    @test observe(lep, statistics, R, mean, 1) == quantities.R
    @test observe(lep, statistics, R, std, 1) == fill(1e-6, 2, 2, 2)
    @test observe(lep, statistics, Z, std, 1) ≈ fill(hypot(1e-6, 2e-6), 2, 2, 2)
    @test DataFrame(observables(lep, (request, (statistics, R, std, 1));
        length_unit=:base, clip=false)).R[1:2] == [quantities.R[1], 1e-6]
    @test_throws ArgumentError observables(lep, ((statistics, R, median, 1),))
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

    RB=LineCableModels.ReportBuilder
    metadata=(port_order=["one","two"],formulation=(equation=:fixture,),axes=nothing)
    definition=RB.BenchmarkTableDefinition(((statistics,R,mean),(statistics,R,std),(statistics,B,mean)))
    artifact=RB.report(definition,(reference=(result=source,metadata=metadata),candidate=(result=source,metadata=metadata)))
    @test length(artifact.table.features) == 3
    @test names(first(artifact.table.features).relative) == ["formula","all","dc","harmonic","narrow","wide"]
    @test size(artifact.table.sampling,1) == 2
    @test all(iszero,skipmissing(artifact.table.terms.relative_rms_percent))
    # Coordinate-specific precision must select std, not mean or another entry.
    # Distinct values above make any transpose/statistic swap observable.
    quantities_by_name=Dict(:R=>R,:B=>B)
    for row in eachrow(artifact.table.mean_sampling_precision)
        k=only(findall(==(row.frequency_Hz),frequencies))
        expected=observe(source,statistics,quantities_by_name[row.quantity],std,row.point)
        @test row.mean_standard_error≈expected[row.row,row.column,k]/sqrt(2)
        @test row.response==metadata.port_order[row.row]
        @test row.excitation==metadata.port_order[row.column]
    end
    @test all(ismissing,artifact.table.sampling.std_sampling_precision)
    @test Set(artifact.table.sampling.method)==Set(artifact.table.formulations.label)
end
