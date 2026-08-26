if lowercase(get(ENV, "LINECABLEMODELS_UPDATE_PLOT_REFERENCES", "false")) == "true"
    using LineCableModels
    using CairoMakie
    using Measurements: measurement
    using LineCableModels.DataModel: CablePosition, LineCableSystem
    using LineCableModels.EarthProps: EarthModel

    include(joinpath(@__DIR__, "..", "support", "golden_fixtures.jl"))
    using .GoldenFixtures: custom_layout_plot

    reference_directory = joinpath(@__DIR__, "..", "fixtures", "golden")
    requested_reference = get(ENV, "LINECABLEMODELS_PLOT_REFERENCE", "")

    function save_reference(name, plot_handle)
        isempty(requested_reference) || name == requested_reference || return nothing
        CairoMakie.save(joinpath(reference_directory, "$name.png"), plot_handle.figure)
        return nothing
    end

    frequency = [50.0, 100.0, 500.0]
    omega = reshape(2π .* frequency, 1, 1, :)
    resistance_values = reshape([1.0, 0.2, 0.2, 2.0], 2, 2, 1) .*
                        ones(1, 1, length(frequency)) .* 1.0e-4
    inductance_values = fill(2.0e-7, 2, 2, length(frequency))
    conductance_values = fill(3.0e-9, 2, 2, length(frequency))
    capacitance_values = fill(4.0e-10, 2, 2, length(frequency))
    parameters = LineParameters(
        complex.(resistance_values, inductance_values .* omega),
        complex.(conductance_values, capacitance_values .* omega),
        frequency
    )

    save_reference(
        "line_rlcg",
        first(plot(parameters, (R, L, G, C); backend = :cairo, display_plot = false))
    )
    save_reference(
        "line_zy_cartesian",
        first(plot(parameters; backend = :cairo, display_plot = false))
    )
    save_reference(
        "line_zy_polar",
        first(plot(parameters, (abs, angle); backend = :cairo, display_plot = false))
    )

    measurement_parameters = LineParameters(
        complex.(measurement.(resistance_values, resistance_values .* 0.05),
            measurement.(inductance_values .* omega, inductance_values .* omega .* 0.05)),
        complex.(measurement.(conductance_values, conductance_values .* 0.05),
            measurement.(capacitance_values .* omega, capacitance_values .* omega .* 0.05)),
        frequency
    )
    save_reference(
        "line_measurements",
        first(plot(
            measurement_parameters,
            (R, L, G, C);
            backend = :cairo,
            display_plot = false
        ))
    )

    summary = SampleSummary([1.0, 2.0, 3.0, 4.0])
    histogram = HistogramDensity([1.0, 3.0, 5.0], [0.25, 0.25])
    representation = CableConstants(2.5, 2.5, 2.5)
    statistics_value = RLC(summary, summary, summary)
    samples_value = RLC(
        [1.0, 2.0, 3.0, 4.0],
        [1.0, 2.0, 3.0, 4.0],
        [1.0, 2.0, 3.0, 4.0]
    )
    histograms_value = RLC(histogram, histogram, histogram)
    mc_formulation = MonteCarlo(
        Formulation(); trials = 4, seed = 1,
        return_samples = true, return_histograms = true
    )
    mc_result = MonteCarloResult(
        mc_formulation,
        [representation],
        [statistics_value],
        [samples_value],
        [histograms_value],
        UInt64(1),
        UInt64[1],
        [4]
    )
    for mode in (:hist, :pdf, :ecdf, :qq)
        save_reference(
            "mc_$mode",
            plot(
                mc_result,
                :R;
                mode,
                data = :both,
                backend = :cairo,
                display_plot = false
            )
        )
    end

    library = CablesLibrary()
    load!(library;
        file_name = joinpath(
            pkgdir(LineCableModels),
            "test",
            "fixtures",
            "data",
            "mv_cable_design.json"
        ))
    design = first(values(library.data))
    save_reference(
        "cable_preview",
        preview(design; backend = :cairo, display_plot = false)
    )
    save_reference(
        "cable_preview_compact",
        preview(
            design;
            size = (900, 350),
            backend = :cairo,
            display_plot = false
        )
    )

    position = CablePosition(
        design,
        0.0,
        -0.20,
        Dict(component.id => (index == 1 ? 1 : 0)
        for
        (index, component) in enumerate(design.components))
    )
    system = LineCableSystem("reference-system", 1000.0, position)
    earth = EarthModel(100.0, 10.0, 1.0)
    save_reference(
        "system_preview",
        preview(
            system;
            earth_model = earth,
            backend = :cairo,
            display_plot = false
        )
    )
    save_reference(
        "material_scale",
        LineCableModels.DataModel.show_material_scale(
            backend = :cairo,
            display_plot = false
        )
    )
    save_reference(
        "custom_layout",
        custom_layout_plot()
    )
end
