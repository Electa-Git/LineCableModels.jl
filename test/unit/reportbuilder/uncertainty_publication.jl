@testitem "UQ / publication / selected matrix coordinates retain axes and units" tags=[:unit] begin
    using DataFrames
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
end
