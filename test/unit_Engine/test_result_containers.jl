@testitem "Result containers: numeric parity and selector grammar" setup = [defaults] begin
    using DataFrames

    library = CablesLibrary()
    load!(library; file_name = joinpath(pkgdir(LineCableModels), "test", "cable_test.json"))
    design = first(values(library.data))
    constants = CableConstants(design)

    # Plain numerical fixture captured before the result-container migration.
    @test constants.R ≈ 2.7567652874268654e-5 rtol = 2eps()
    @test constants.L ≈ 2.8718381083175005e-7 rtol = 2eps()
    @test constants.C ≈ 4.1335723330313053e-10 rtol = 2eps()
    @test basis(constants) === :per_length
    @test resistance(constants) === constants.R
    @test inductance(constants) === constants.L
    @test capacitance(constants) === constants.C

    baseparams = DataFrame(design, :baseparams)
    @test baseparams.computed ≈ [
        constants.R * 1.0e3,
        constants.L * 1.0e6,
        constants.C * 1.0e9
    ]

    frequency = [50.0, 100.0, 200.0]
    omega = reshape(2π .* frequency, 1, 1, :)
    resistance_values = reshape(collect(1.0:12.0), 2, 2, 3) .* 1.0e-4
    inductance_values = reshape(collect(13.0:24.0), 2, 2, 3) .* 1.0e-7
    conductance_values = reshape(collect(25.0:36.0), 2, 2, 3) .* 1.0e-8
    capacitance_values = reshape(collect(37.0:48.0), 2, 2, 3) .* 1.0e-10
    impedance = complex.(resistance_values, inductance_values .* omega)
    admittance = complex.(conductance_values, capacitance_values .* omega)
    parameters = LineParameters(impedance, admittance, frequency; basis = :total)

    @test @inferred(basis(parameters)) === :total
    @test @inferred(frequencies(parameters)) == frequency
    @test @inferred(nconductors(parameters)) == 2
    @test @inferred(nfrequencies(parameters)) == 3
    @test Z(parameters) == impedance
    @test Y(parameters) == admittance
    @test Z(parameters, 1, 2) == impedance[1, 2, :]
    @test Z(parameters, 1, 2, :) == impedance[1, 2, :]
    @test Z(parameters, 1, 2, 2:3) == impedance[1, 2, 2:3]
    @test Z(parameters, 1, 2, 2) == impedance[1, 2, 2]
    @test Y(parameters, 2, 1, 1:2) == admittance[2, 1, 1:2]
    @test R(parameters, 1, 2, :) == resistance_values[1, 2, :]
    @test X(parameters, 1, 2, :) == imag.(impedance[1, 2, :])
    @test L(parameters, 1, 2, :) ≈ inductance_values[1, 2, :]
    @test G(parameters, 2, 1) == conductance_values[2, 1, :]
    @test B(parameters, 2, 1) == imag.(admittance[2, 1, :])
    @test C(parameters, 2, 1) ≈ capacitance_values[2, 1, :]
    @test series_impedance(parameters) === parameters.Z
    @test shunt_admittance(parameters) === parameters.Y
    @test resistance(parameters, 1, 1) == R(parameters, 1, 1)
    @test reactance(parameters, 1, 1) == X(parameters, 1, 1)
    @test inductance(parameters, 1, 1) == L(parameters, 1, 1)
    @test conductance(parameters, 1, 1) == G(parameters, 1, 1)
    @test susceptance(parameters, 1, 1) == B(parameters, 1, 1)
    @test capacitance(parameters, 1, 1) == C(parameters, 1, 1)
    @test abs.(Z(parameters, 1, 1)) == abs.(impedance[1, 1, :])
    @test angle.(Y(parameters, 1, 1)) == angle.(admittance[1, 1, :])

    selected = @inferred parameters[2:3]
    @test basis(selected) === :total
    @test domain(selected) === domain(parameters)
    @test frequencies(selected) == frequency[2:3]
    @test selected.Z.values == impedance[:, :, 2:3]
    @test basis(parameters[2]) === :total
    @test frequencies(parameters[:]) == frequency

    series = SeriesImpedance(impedance; basis = :total)
    shunt = ShuntAdmittance(admittance; basis = :total)
    @test Z(series, 1, 1, :) == impedance[1, 1, :]
    @test R(series, 1, 1) == resistance_values[1, 1, :]
    @test X(series, 1, 1) == imag.(impedance[1, 1, :])
    @test Y(shunt, 1, 1, :) == admittance[1, 1, :]
    @test G(shunt, 1, 1) == conductance_values[1, 1, :]
    @test B(shunt, 1, 1) == imag.(admittance[1, 1, :])

    series_frames, shunt_frames = DataFrame(parameters)
    @test DataFrames.metadata(series_frames[1, 1], "units")[:R] == "Ω"
    @test DataFrames.metadata(series_frames[1, 1], "units")[:L] == "mH"
    @test DataFrames.metadata(shunt_frames[1, 1], "units")[:C] == "μF"

    zero_frequency = LineParameters(
        impedance[:, :, 1:1],
        admittance[:, :, 1:1],
        [0.0]
    )
    @test R(zero_frequency, 1, 1, 1) == real(impedance[1, 1, 1])
    @test X(zero_frequency, 1, 1, 1) == imag(impedance[1, 1, 1])
    @test G(zero_frequency, 1, 1, 1) == real(admittance[1, 1, 1])
    @test B(zero_frequency, 1, 1, 1) == imag(admittance[1, 1, 1])
    @test_throws DomainError L(zero_frequency)
    @test_throws DomainError C(zero_frequency, 1, 1, 1)

    @test_throws DimensionMismatch LineParameters(
        zeros(ComplexF64, 2, 3, 1),
        zeros(ComplexF64, 2, 3, 1),
        [50.0]
    )
    @test_throws DimensionMismatch LineParameters(
        zeros(ComplexF64, 2, 2, 1),
        zeros(ComplexF64, 3, 3, 1),
        [50.0]
    )
    @test_throws DimensionMismatch LineParameters(
        zeros(ComplexF64, 2, 2, 2),
        zeros(ComplexF64, 2, 2, 2),
        [50.0]
    )
    @test_throws ArgumentError LineParameters(
        zeros(ComplexF64, 1, 1, 1),
        zeros(ComplexF64, 1, 1, 1),
        [Inf]
    )
    @test_throws ArgumentError LineParameters(
        SeriesImpedance(zeros(ComplexF64, 1, 1, 1); basis = :per_length),
        ShuntAdmittance(zeros(ComplexF64, 1, 1, 1); basis = :total),
        [50.0]
    )
end

@testitem "Monte Carlo containers: optional storage and joint trials" setup = [defaults] begin
    using Random
    using Measurements: measurement

    values = [1.0, 2.0, 3.0, 4.0, 5.0]
    summary = SampleSummary(values)
    @test summary == SampleSummary(3.0, sqrt(2.5), 1.0, 1.2, 3.0, 4.8, 5.0)
    @test SampleSummary([1, 2, 3]).mean === 2.0

    model = HistogramPDF([1.0, 3.0, 5.0], [0.25, 0.25])
    @test HistogramPDF([0, 1], [2]) isa HistogramPDF{Float64}
    @test_throws DomainError quantile(model, -0.1)
    constants_statistics = CableConstants(summary, summary, summary)
    constants_samples = CableConstants(values, 2values, 3values)
    constants_distributions = CableConstants(model, model, model)
    constants_surrogate = CableConstants(
        measurement(3.0, 1.0),
        measurement(6.0, 2.0),
        measurement(9.0, 3.0)
    )
    constants_result = CableConstantsMC(
        constants_statistics,
        constants_samples,
        constants_distributions,
        constants_surrogate,
        5,
        0.95
    )

    @test @inferred(statistics(constants_result)) === constants_statistics
    @test @inferred(statistics(constants_result, :R)) === summary
    @test @inferred(mean(constants_result, :R)) == 3.0
    @test @inferred(std(constants_result, :R)) == sqrt(2.5)
    @test quantile(constants_result, :R, 0.05) == 1.2
    @test quantile(constants_result, :R, 0.50) == 3.0
    @test quantile(constants_result, :R, 0.95) == 4.8
    @test quantile(constants_result, :R, 0.25) == 2.0
    @test has_samples(constants_result)
    @test has_distributions(constants_result)
    @test samples(constants_result, :L) == 2values
    @test distribution(constants_result, :C) === model
    @test surrogate(constants_result) === constants_surrogate
    @test basis(constants_result) === :per_length
    @test ntrials(constants_result) == 5
    @test confidence(constants_result) == 0.95
    @test @inferred(trial(constants_result, 3)) == CableConstants(3.0, 6.0, 9.0)
    @test rand(MersenneTwister(42), constants_result) isa CableConstants
    @test_throws BoundsError trial(constants_result, 0)
    @test_throws ArgumentError CableConstantsMC(
        constants_statistics,
        (R = values, L = values, C = values),
        constants_distributions,
        constants_surrogate,
        5,
        0.95
    )

    no_storage = CableConstantsMC(
        constants_statistics,
        nothing,
        nothing,
        constants_surrogate,
        5,
        0.95
    )
    @test !has_samples(no_storage)
    @test !has_distributions(no_storage)
    @test_throws ArgumentError samples(no_storage, :R)
    @test_throws ArgumentError distribution(no_storage, :R)
    @test_throws ArgumentError quantile(no_storage, :R, 0.25)
    @test_throws ArgumentError trial(no_storage, 1)
    @test_throws ArgumentError rand(MersenneTwister(1), no_storage)

    frequency = [50.0, 100.0]
    sample_count = 5
    resistance_samples = reshape(collect(1.0:10.0), 1, 1, 2, sample_count) .* 1.0e-3
    inductance_samples = resistance_samples .* 1.0e-3
    capacitance_samples = resistance_samples .* 1.0e-6
    conductance_samples = resistance_samples .* 1.0e-4
    summarize(samples_array) = map(
        index -> SampleSummary(view(samples_array, index.I..., :)),
        CartesianIndices(size(samples_array)[1:3])
    )
    line_statistics = RLCG(
        summarize(resistance_samples),
        summarize(inductance_samples),
        summarize(capacitance_samples),
        summarize(conductance_samples)
    )
    line_samples = RLCG(
        resistance_samples,
        inductance_samples,
        capacitance_samples,
        conductance_samples
    )
    line_distributions = RLCG(
        fill(model, 1, 1, 2),
        fill(model, 1, 1, 2),
        fill(model, 1, 1, 2),
        fill(model, 1, 1, 2)
    )
    mean_R = mean(resistance_samples; dims = 4)[:, :, :, 1]
    mean_L = mean(inductance_samples; dims = 4)[:, :, :, 1]
    mean_C = mean(capacitance_samples; dims = 4)[:, :, :, 1]
    mean_G = mean(conductance_samples; dims = 4)[:, :, :, 1]
    omega = reshape(2π .* frequency, 1, 1, :)
    surrogate_parameters = LineParameters(
        complex.(measurement.(mean_R, 0.0), measurement.(mean_L .* omega, 0.0)),
        complex.(measurement.(mean_G, 0.0), measurement.(mean_C .* omega, 0.0)),
        frequency;
        basis = :total
    )
    line_result = LineParametersMC(
        line_statistics,
        line_samples,
        line_distributions,
        surrogate_parameters,
        sample_count,
        0.90
    )

    @test basis(line_result) === :total
    @test domain(line_result) === PhaseDomain
    @test frequencies(line_result) == frequency
    @test nconductors(line_result) == 1
    @test nfrequencies(line_result) == 2
    @test statistics(line_result, :R, 1, 1, 2) == line_statistics.R[1, 1, 2]
    @test mean(line_result, :L, 1, 1, :) == getproperty.(line_statistics.L[1, 1, :], :mean)
    @test samples(line_result, :C, 1, 1, 2) == capacitance_samples[1, 1, 2, :]
    @test distribution(line_result, :G, 1, 1, 1) === model
    @test quantile(line_result, :R, 0.95, 1, 1, 1) == line_statistics.R[1, 1, 1].q95
    @test quantile(line_result, :R, 0.25, 1, 1, 1) ==
          quantile(resistance_samples[1, 1, 1, :], 0.25)

    retained_trial = @inferred trial(line_result, 4)
    @test basis(retained_trial) === :total
    @test domain(retained_trial) === domain(line_result)
    @test frequencies(retained_trial) == frequency
    @test R(retained_trial) == resistance_samples[:, :, :, 4]
    @test L(retained_trial) ≈ inductance_samples[:, :, :, 4]
    @test C(retained_trial) ≈ capacitance_samples[:, :, :, 4]
    @test G(retained_trial) == conductance_samples[:, :, :, 4]
    @test rand(MersenneTwister(9), line_result) isa LineParameters

    distribution_only = LineParametersMC(
        line_statistics,
        nothing,
        line_distributions,
        surrogate_parameters,
        sample_count,
        0.90
    )
    @test quantile(distribution_only, :R, 0.25, 1, 1, 1) == 2.0
    @test_throws ArgumentError trial(distribution_only, 1)
    @test_throws ArgumentError statistics(line_result, :Z)
    @test_throws ArgumentError samples(line_result, :Z)
    @test_throws ArgumentError distribution(line_result, :Z)
    @test_throws ArgumentError LineParametersMC(
        RLCG(zeros(1, 1, 2), zeros(1, 1, 2), zeros(1, 1, 2), zeros(1, 1, 2)),
        nothing,
        nothing,
        surrogate_parameters,
        sample_count,
        0.90
    )
end

@testitem "PlotBuilder: result RenderSpec semantics" setup = [defaults] begin
    frequency = [50.0, 500.0, 900.0]
    resistance_values = reshape([1.0, 0.2, 0.2, 2.0], 2, 2, 1) .*
                        ones(1, 1, length(frequency)) .* 1.0e-4
    inductance_values = fill(2.0e-7, 2, 2, length(frequency))
    conductance_values = fill(3.0e-9, 2, 2, length(frequency))
    capacitance_values = fill(4.0e-10, 2, 2, length(frequency))
    omega = reshape(2π .* frequency, 1, 1, :)
    parameters = LineParameters(
        complex.(resistance_values, inductance_values .* omega),
        complex.(conductance_values, capacitance_values .* omega),
        frequency
    )

    render = LineCableModels.PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        parameters;
        mode = :RLCG,
        coord = :cart,
        con = (1:2, 1:2)
    )
    @test render.spec === LineCableModels.Engine.LineParameterPlotSpec
    @test length(render.figures) == 4
    @test getproperty.(only.(getproperty.(render.figures, :views)), :title) == [
        "Series resistance",
        "Series inductance",
        "Shunt conductance",
        "Shunt capacitance"
    ]
    @test all(page -> page.layout === :single, render.figures)
    @test all(page -> only(page.views).xaxis.label == "Frequency [Hz]", render.figures)
    @test only(render.figures[1].views).yaxis.label == "Series resistance [Ω/km]"
    @test only(render.figures[2].views).yaxis.label == "Series inductance [mH/km]"
    @test length(only(render.figures[1].views).series) == 4
    @test all(series -> series.kind === :line, only(render.figures[1].views).series)
    @test first(only(render.figures[1].views).series).label == "R[1,1]"
    @test render.figures[1].kwargs.controls ==
          LineCableModels.PlotBuilder.control_definitions()
    @test render.figures[1].kwargs.configuration.mode === :RLCG
    @test render.figures[1].kwargs.configuration.coord === :cart
    @test render.figures[1].kwargs.configuration.length_unit === :kilo
    @test render.figures[1].kwargs.configuration.conductors == (1:2, 1:2)

    polar = LineCableModels.PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        parameters;
        mode = :ZY,
        coord = :polar
    )
    @test length(polar.figures) == 4
    @test occursin("impedance magnitude", lowercase(only(polar.figures[1].views).yaxis.label))
    @test occursin("angle", lowercase(only(polar.figures[2].views).yaxis.label))

    summary = SampleSummary([1.0, 2.0, 3.0, 4.0])
    model = HistogramPDF([1.0, 3.0, 5.0], [0.25, 0.25])
    mc_result = CableConstantsMC(
        CableConstants(summary, summary, summary),
        CableConstants([1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]),
        CableConstants(model, model, model),
        CableConstants(2.5, 2.5, 2.5),
        4,
        0.95
    )
    for mode in (:hist, :pdf, :ecdf, :qq)
        mc_render = LineCableModels.PlotBuilder.make_render(
            LineCableModels.UQ.MCDistributionPlotSpec,
            mc_result;
            quantity = :R,
            mode,
            data = :both
        )
        @test length(mc_render.figures) == 1
        @test only(mc_render.figures).kwargs.mode === mode
        @test only(mc_render.figures).kwargs.controls.export_svg
        @test only(mc_render.figures).kwargs.configuration.mode === mode
        @test !isempty(only(only(mc_render.figures).views).series)
    end
end
