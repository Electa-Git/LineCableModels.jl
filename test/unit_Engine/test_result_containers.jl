@testitem "Result containers: numeric parity and selector grammar" setup = [defaults] begin
    using DataFrames

    library = CablesLibrary()
    load!(library; file_name = joinpath(pkgdir(LineCableModels), "test", "cable_test.json"))
    design = first(values(library.data))
    constants = @inferred CableConstants(design)
    @test @inferred(CableConstants(design; S = 0.25, rho_e = 150.0)) isa
          CableConstants{Float64}

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
    parameters = LineParameters(
        impedance,
        admittance,
        frequency;
        basis = :total
    )
    typed_series = @inferred SeriesImpedance{ComplexF64, :total}(copy(impedance))
    typed_shunt = @inferred ShuntAdmittance{ComplexF64, :total}(copy(admittance))
    typed_parameters = @inferred LineParameters(typed_series, typed_shunt, frequency)
    @test basis(typed_parameters) === :total

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
    @test R(parameters, 1, 2, 2:3) == resistance_values[1, 2, 2:3]
    @test X(parameters, 1, 2, 2) == imag(impedance[1, 2, 2])
    @test L(parameters, 1, 2, 2:3) ≈ inductance_values[1, 2, 2:3]
    @test G(parameters, 2, 1, :) == conductance_values[2, 1, :]
    @test B(parameters, 2, 1, 2) == imag(admittance[2, 1, 2])
    @test C(parameters, 2, 1, 1:2) ≈ capacitance_values[2, 1, 1:2]
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

    unit_handler = LineCableModels.UnitHandler
    @test unit_handler.line_components(:series, :RLCG, :cart) == (:R, :L)
    @test unit_handler.line_components(:shunt, :RLCG, :polar) == (:G, :C)
    @test unit_handler.line_components(:series, :ZY, :polar) == (:Z_abs, :Z_angle)
    @test unit_handler.line_components(:shunt, :ZY, :cart) == (:Y_re, :Y_im)
    resistance_quantity = unit_handler.line_component_quantity(:R)
    @test resistance_quantity.semantic === :resistance
    @test resistance_quantity.accessor === R
    @test resistance_quantity.selector === nothing
    @test resistance_quantity.unit_name === :ohm
    @test resistance_quantity.prefix === :base
    @test unit_handler.quantity(Z) isa unit_handler.QuantityTag{:impedance}
    @test unit_handler.quantity(Y) isa unit_handler.QuantityTag{:admittance}
    @test unit_handler.quantity(R) isa unit_handler.QuantityTag{:resistance}
    @test unit_handler.quantity(L) isa unit_handler.QuantityTag{:inductance}
    @test unit_handler.quantity(Z, :re) isa unit_handler.QuantityTag{:resistance}
    @test unit_handler.quantity(Z, :im) isa unit_handler.QuantityTag{:reactance}
    @test unit_handler.quantity(Z, :abs) isa
          unit_handler.QuantityTag{(:impedance, :abs)}
    @test unit_handler.quantity(Y, :angle) isa
          unit_handler.QuantityTag{(:admittance, :angle)}
    @test unit_handler.line_component_quantity(:Z_abs).accessor === Z
    @test unit_handler.line_component_quantity(:Z_abs).selector === :abs
    inductance_unit = unit_handler.line_component_unit(:L, :per_length)
    @test unit_handler.get_label(inductance_unit.units) == "mH/km"
    @test inductance_unit.scale == 1.0e6
    total_resistance_unit = unit_handler.line_component_unit(:R, :total)
    @test unit_handler.get_label(total_resistance_unit.units) == "Ω"
    @test total_resistance_unit.scale == 1.0
    @test !isdefined(unit_handler, :line_component_values)
    @test_throws ArgumentError unit_handler.line_components(:invalid, :ZY, :cart)
    @test_throws ArgumentError unit_handler.line_component_quantity(:invalid)
    @test_throws ArgumentError unit_handler.quantity(Z, :magnitude)
    @test_throws ArgumentError unit_handler.quantity(R, :abs)
    @test_throws ArgumentError unit_handler.quantity(sin)
    @test_throws ArgumentError unit_handler.line_component_unit(
        :R,
        :per_length;
        quantity_units = 42
    )
    @test_throws ArgumentError unit_handler.line_component_unit(
        :R,
        :per_length;
        length_unit = :kilometer
    )
    @test_throws ArgumentError unit_handler.line_component_unit(
        :R,
        :per_length;
        quantity_units = (; R = 1)
    )

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
    @test_throws ArgumentError SeriesImpedance{ComplexF64, :invalid}(
        zeros(ComplexF64, 1, 1, 1)
    )
    @test_throws ArgumentError ShuntAdmittance{ComplexF64, 1}(
        zeros(ComplexF64, 1, 1, 1)
    )

    series_frames, shunt_frames = DataFrame(parameters)
    @test DataFrames.metadata(series_frames[1, 1], "units")[:R] == "Ω"
    @test DataFrames.metadata(series_frames[1, 1], "units")[:L] == "mH"
    @test DataFrames.metadata(shunt_frames[1, 1], "units")[:C] == "μF"
    @test_throws ArgumentError DataFrame(parameters; tol = -1.0)

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
    @test_throws DomainError DataFrame(zero_frequency)

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
    import Distributions

    values = [1.0, 2.0, 3.0, 4.0, 5.0]
    summary = @inferred SampleSummary(values)
    @test summary == SampleSummary(3.0, sqrt(2.5), 1.0, 1.2, 3.0, 4.8, 5.0)
    @test SampleSummary([1, 2, 3]).mean === 2.0
    @test_throws ArgumentError SampleSummary(Float64[])
    @test_throws ArgumentError SampleSummary([1.0, Inf])
    @test_throws ArgumentError SampleSummary(2.0, -1.0, 1.0, 1.1, 2.0, 2.9, 3.0)
    @test_throws ArgumentError SampleSummary(2.0, 1.0, 1.0, 2.1, 2.0, 2.9, 3.0)

    model = @inferred HistogramPDF([1.0, 3.0, 5.0], [0.25, 0.25])
    @test HistogramPDF([0, 1], [2]) isa HistogramPDF{Float64}
    @test Distributions.pdf(model, 2.0) == 0.25
    @test Distributions.pdf(model, 6.0) == 0.0
    @test Distributions.cdf(model, 0.0) == 0.0
    @test Distributions.cdf(model, 3.0) == 0.5
    @test Distributions.cdf(model, 6.0) == 1.0
    @test quantile(model, 0.25) == 2.0
    @test isfinite(rand(MersenneTwister(91), model))
    @test_throws DomainError quantile(model, -0.1)
    @test_throws ArgumentError HistogramPDF([0.0, 1.0], Float64[])
    @test_throws ArgumentError HistogramPDF([0.0, 1.0, 2.0], [1.0])
    @test_throws ArgumentError HistogramPDF([0.0, 1.0], [-1.0])
    @test_throws ArgumentError HistogramPDF([0.0, 1.0], [0.0])
    @test_throws ArgumentError HistogramPDF([0.0, 0.0], [1.0])
    @test_throws ArgumentError HistogramPDF([0.0, Inf], [1.0])
    @test_throws ArgumentError HistogramPDF{Float64}([0.0, 0.0], [1.0])
    constants_statistics = CableConstants(summary, summary, summary)
    constants_samples = CableConstants(values, 2values, 3values)
    constants_distributions = CableConstants(model, model, model)
    constants_surrogate = CableConstants(
        measurement(3.0, 1.0),
        measurement(6.0, 2.0),
        measurement(9.0, 3.0)
    )
    constants_result = @inferred CableConstantsMC(
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
    @test @inferred(surrogate(constants_result)) === constants_surrogate
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
    @test_throws ArgumentError CableConstantsMC(
        constants_statistics,
        constants_samples,
        constants_distributions,
        CableConstants(3.0, 6.0, 9.0),
        5,
        0.95
    )

    no_storage = @inferred CableConstantsMC(
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

    constants_samples_only = @inferred CableConstantsMC(
        constants_statistics,
        constants_samples,
        nothing,
        constants_surrogate,
        5,
        0.95
    )
    @test has_samples(constants_samples_only)
    @test !has_distributions(constants_samples_only)
    @test @inferred(samples(constants_samples_only, :R)) === values
    @test @inferred(trial(constants_samples_only, 2)) == CableConstants(2.0, 4.0, 6.0)

    constants_distributions_only = @inferred CableConstantsMC(
        constants_statistics,
        nothing,
        constants_distributions,
        constants_surrogate,
        5,
        0.95
    )
    @test !has_samples(constants_distributions_only)
    @test has_distributions(constants_distributions_only)
    @test @inferred(distribution(constants_distributions_only, :L)) === model
    @test quantile(constants_distributions_only, :L, 0.25) == 2.0
    @test_throws ArgumentError trial(constants_distributions_only, 1)

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
    line_result = @inferred LineParametersMC(
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
    @test @inferred(statistics(line_result, :R, 1, 1, 2)) ==
          line_statistics.R[1, 1, 2]
    @test mean(line_result, :L, 1, 1, :) == getproperty.(line_statistics.L[1, 1, :], :mean)
    @test @inferred(samples(line_result, :C, 1, 1, 2)) ==
          capacitance_samples[1, 1, 2, :]
    @test @inferred(distribution(line_result, :G, 1, 1, 1)) === model
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

    distribution_only = @inferred LineParametersMC(
        line_statistics,
        nothing,
        line_distributions,
        surrogate_parameters,
        sample_count,
        0.90
    )
    @test quantile(distribution_only, :R, 0.25, 1, 1, 1) == 2.0
    @test_throws ArgumentError trial(distribution_only, 1)

    samples_only = @inferred LineParametersMC(
        line_statistics,
        line_samples,
        nothing,
        surrogate_parameters,
        sample_count,
        0.90
    )
    @test has_samples(samples_only)
    @test !has_distributions(samples_only)
    @test @inferred(trial(samples_only, 1)) isa LineParameters
    @test quantile(samples_only, :R, 0.25, 1, 1, 1) ==
          quantile(resistance_samples[1, 1, 1, :], 0.25)

    line_no_storage = @inferred LineParametersMC(
        line_statistics,
        nothing,
        nothing,
        surrogate_parameters,
        sample_count,
        0.90
    )
    @test !has_samples(line_no_storage)
    @test !has_distributions(line_no_storage)
    @test @inferred(surrogate(line_no_storage)) === surrogate_parameters
    @test_throws ArgumentError quantile(line_no_storage, :R, 0.25, 1, 1, 1)
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
    @test_throws ArgumentError LineParametersMC(
        line_statistics,
        line_samples,
        line_distributions,
        LineParameters(
            complex.(mean_R, mean_L .* omega),
            complex.(mean_G, mean_C .* omega),
            frequency;
            basis = :total
        ),
        sample_count,
        0.90
    )
end

@testitem "PlotBuilder: result RenderSpec semantics" setup = [defaults] begin
    plotspec_source = read(
        joinpath(dirname(pathof(LineCableModels)), "engine", "plotspecs.jl"),
        String
    )
    @test !occursin("source.values", plotspec_source)
    @test !occursin("UnitHandler.line_component_values", plotspec_source)

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
    @test all(page -> page.layout.name === :single, render.figures)
    @test all(page -> only(page.views).xaxis.label == "Frequency [Hz]", render.figures)
    @test only(render.figures[1].views).yaxis.label == "Series resistance [Ω/km]"
    @test only(render.figures[2].views).yaxis.label == "Series inductance [mH/km]"
    @test length(only(render.figures[1].views).series) == 4
    @test all(series -> series.kind === :line, only(render.figures[1].views).series)
    first_resistance = first(only(render.figures[1].views).series)
    @test first_resistance.label == "R[1,1]"
    @test first_resistance.xdata == frequency
    @test first_resistance.ydata ≈ vec(resistance_values[1, 1, :]) .* 1.0e3
    @test first_resistance.attributes.linewidth == 2
    @test getproperty.(only(render.figures[1].views).series, :label) ==
          ["R[1,1]", "R[1,2]", "R[2,1]", "R[2,2]"]
    @test render.figures[1].controls.reset
    @test render.figures[1].controls.export_svg
    @test only(render.figures[1].views).xaxis.allowed_scales == (:linear, :log10)
    @test only(render.figures[1].views).yaxis.allowed_scales == (:linear, :log10)
    @test render.figures[1].key.mode === :RLCG
    @test render.figures[1].key.coordinates === :cart
    @test render.figures[1].key.conductors == (1:2, 1:2)
    @test all(page -> page.export_spec.theme === :default, render.figures)
    @test all(page -> page.export_spec.open_file, render.figures)
    conductance_page = render.figures[3]
    first_conductance = first(only(conductance_page.views).series)
    @test only(conductance_page.views).yaxis.exponent == -6
    @test only(conductance_page.views).yaxis.label == "Shunt conductance [S/km]"
    @test first_conductance.ydata == fill(3.0e-6, length(frequency))

    series_render = LineCableModels.PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        parameters.Z;
        frequencies = frequency,
        mode = :ZY,
        coord = :cart,
        freq_unit = :kilo,
        length_unit = :base,
        quantity_units = (; resistance = :milli),
        con = ([1], [2]),
        xscale = :log10,
        yscale = :log10
    )
    @test length(series_render.figures) == 2
    @test all(page -> length(only(page.views).series) == 1, series_render.figures)
    @test first(only(series_render.figures[1].views).series).label == "R[1,2]"
    @test only(series_render.figures[1].views).xaxis.label == "Frequency [kHz]"
    @test only(series_render.figures[1].views).yaxis.label ==
          "Series resistance [mΩ/m]"
    @test only(series_render.figures[1].views).xaxis.scale === :log10
    @test only(series_render.figures[1].views).yaxis.scale === :log10

    publication_render = LineCableModels.PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        parameters;
        export_theme = :publication,
        open_export = false
    )
    @test all(page -> page.export_spec.theme === :publication, publication_render.figures)
    @test all(page -> !page.export_spec.open_file, publication_render.figures)
    @test_throws ArgumentError LineCableModels.PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        parameters;
        export_theme = :unsupported
    )

    shunt_render = LineCableModels.PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        parameters.Y;
        frequencies = frequency,
        mode = :ZY,
        coord = :polar,
        con = (2, :)
    )
    @test length(shunt_render.figures) == 2
    @test all(page -> length(only(page.views).series) == 2, shunt_render.figures)
    @test occursin(
        "admittance magnitude",
        lowercase(only(shunt_render.figures[1].views).yaxis.label)
    )

    polar = LineCableModels.PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        parameters;
        mode = :ZY,
        coord = :polar
    )
    @test length(polar.figures) == 4
    @test occursin("impedance magnitude", lowercase(only(polar.figures[1].views).yaxis.label))
    @test occursin("angle", lowercase(only(polar.figures[2].views).yaxis.label))
    @test_throws ArgumentError LineCableModels.PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        parameters.Z;
        frequencies = [50.0, NaN, 900.0]
    )
    @test_throws DomainError LineCableModels.PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        parameters.Z;
        frequencies = [0.0, 50.0, 900.0],
        xscale = :log10
    )
    @test_throws ArgumentError LineCableModels.PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotSpec,
        parameters.Z;
        frequencies = frequency,
        con = 1
    )

    summary = SampleSummary([1.0, 2.0, 3.0, 4.0])
    model = HistogramPDF([1.0, 3.0, 5.0], [0.25, 0.25])
    mc_result = CableConstantsMC(
        CableConstants(summary, summary, summary),
        CableConstants([1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]),
        CableConstants(model, model, model),
        CableConstants(
            measurement(2.5, 0.0),
            measurement(2.5, 0.0),
            measurement(2.5, 0.0)
        ),
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
        @test only(mc_render.figures).key.mode === mode
        @test only(mc_render.figures).controls.export_svg
        @test only(mc_render.figures).export_spec.theme === :default
        @test only(mc_render.figures).export_spec.open_file
        @test only(only(mc_render.figures).views).xaxis.allowed_scales == (:linear,)
        @test only(only(mc_render.figures).views).yaxis.allowed_scales == (:linear,)
        @test only(only(mc_render.figures).views).xaxis.exponent == 3
        expected_y_exponent = mode in (:hist, :pdf) ? -4 : mode === :qq ? 3 : 0
        @test only(only(mc_render.figures).views).yaxis.exponent == expected_y_exponent
        kinds = getproperty.(only(only(mc_render.figures).views).series, :kind)
        expected_kinds = mode === :hist ? [:histogram, :stairs] :
                         mode === :pdf ? [:stairs] :
                         mode === :ecdf ? [:line, :line] : [:scatter, :line]
        @test kinds == expected_kinds
        if mode === :hist
            histogram_series, pdf_series = only(mc_render.figures).views[1].series
            @test histogram_series.xdata == [1000.0, 2000.0, 3000.0, 4000.0]
            @test histogram_series.attributes.normalization === :pdf
            @test pdf_series.xdata == [1000.0, 3000.0, 5000.0]
            @test pdf_series.ydata == [0.00025, 0.00025, 0.00025]
            @test pdf_series.attributes.color === :red
            @test pdf_series.attributes.linewidth == 2
        elseif mode === :qq
            scatter_series = first(only(mc_render.figures).views[1].series)
            @test scatter_series.attributes.color === :steelblue
            @test scatter_series.attributes.markersize == 6
        end
    end
    @test_throws ArgumentError LineCableModels.PlotBuilder.make_render(
        LineCableModels.UQ.MCDistributionPlotSpec,
        mc_result;
        quantity = :R,
        nbins = 0
    )

    line_samples = reshape(collect(1.0:12.0), 1, 1, 3, 4)
    summarize(values) = map(
        index -> SampleSummary(view(values, index.I..., :)),
        CartesianIndices(size(values)[1:3])
    )
    line_statistics = RLCG(
        summarize(line_samples),
        summarize(line_samples .* 1.0e-3),
        summarize(line_samples .* 1.0e-6),
        summarize(line_samples .* 1.0e-4)
    )
    line_distributions = RLCG(
        fill(model, 1, 1, 3),
        fill(model, 1, 1, 3),
        fill(model, 1, 1, 3),
        fill(model, 1, 1, 3)
    )
    measured_parameters = LineParameters(
        complex.(measurement.(real.(parameters.Z.values), 0.0),
            measurement.(imag.(parameters.Z.values), 0.0)),
        complex.(measurement.(real.(parameters.Y.values), 0.0),
            measurement.(imag.(parameters.Y.values), 0.0)),
        frequency
    )
    line_mc = LineParametersMC(
        line_statistics,
        RLCG(
            line_samples,
            line_samples .* 1.0e-3,
            line_samples .* 1.0e-6,
            line_samples .* 1.0e-4
        ),
        line_distributions,
        LineParameters(
            measured_parameters.Z.values[1:1, 1:1, :],
            measured_parameters.Y.values[1:1, 1:1, :],
            frequency
        ),
        4,
        0.95
    )
    for mode in (:hist, :pdf, :ecdf, :qq)
        line_mc_render = LineCableModels.PlotBuilder.make_render(
            LineCableModels.UQ.MCDistributionPlotSpec,
            line_mc;
            quantity = :R,
            ijk = (1, 1, 2),
            mode,
            data = :both
        )
        @test only(line_mc_render.figures).key.selection == (1, 1, 2)
        @test only(line_mc_render.figures).key.mode === mode
        @test !isempty(only(only(line_mc_render.figures).views).series)
        if mode === :hist
            retained = first(only(line_mc_render.figures).views[1].series)
            @test retained.xdata == collect(line_samples[1, 1, 2, :]) .* 1.0e3
        end
    end
end
