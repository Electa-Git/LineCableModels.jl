using LineCableModels
using Measurements: measurement

function build_manual_plot_gallery(
        backend::Symbol;
        display_plot::Bool,
        export_theme::Symbol = :publication
)
    frequency = [10.0, 50.0, 100.0, 500.0, 1_000.0, 10_000.0]
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

    measurement_parameters = LineParameters(
        complex.(
            measurement.(resistance_values, resistance_values .* 0.05),
            measurement.(inductance_values .* omega, inductance_values .* omega .* 0.05)
        ),
        complex.(
            measurement.(conductance_values, conductance_values .* 0.05),
            measurement.(capacitance_values .* omega, capacitance_values .* omega .* 0.05)
        ),
        frequency
    )

    summary = SampleSummary([1.0, 2.0, 3.0, 4.0])
    distribution_model = HistogramPDF([1.0, 3.0, 5.0], [0.25, 0.25])
    mc_result = CableConstantsMC(
        CableConstants(summary, summary, summary),
        CableConstants(
            [1.0, 2.0, 3.0, 4.0],
            [1.0, 2.0, 3.0, 4.0],
            [1.0, 2.0, 3.0, 4.0]
        ),
        CableConstants(distribution_model, distribution_model, distribution_model),
        CableConstants(
            measurement(2.5, 0.5),
            measurement(2.5, 0.5),
            measurement(2.5, 0.5)
        ),
        4,
        0.95
    )

    gallery = Pair{String, UIPlot}[]
    function add_pages!(title, handles)
        for (index, handle) in enumerate(handles)
            push!(gallery, "$title — page $index" => handle)
        end
        return handles
    end

    add_pages!(
        "Line parameters: RLCG",
        Makie.plot(parameters; mode = :RLCG, backend, display_plot, export_theme)
    )
    add_pages!(
        "Line parameters: Z/Y Cartesian",
        Makie.plot(
            parameters;
            mode = :ZY,
            coord = :cart,
            backend,
            display_plot,
            export_theme
        )
    )
    add_pages!(
        "Line parameters: Z/Y polar",
        Makie.plot(
            parameters;
            mode = :ZY,
            coord = :polar,
            backend,
            display_plot,
            export_theme
        )
    )
    add_pages!(
        "Line parameters: measurement error bars",
        Makie.plot(
            measurement_parameters;
            mode = :RLCG,
            backend,
            display_plot,
            export_theme
        )
    )

    for mode in (:hist, :pdf, :ecdf, :qq)
        push!(
            gallery,
            "Monte Carlo: $mode" => Makie.plot(
                mc_result,
                :R;
                mode,
                data = :both,
                backend,
                display_plot,
                export_theme
            )
        )
    end

    library = CablesLibrary()
    load!(library; file_name = joinpath(pkgdir(LineCableModels), "test", "cable_test.json"))
    design = first(values(library.data))
    push!(
        gallery,
        "Cable preview" => preview(design; backend, display_plot, export_theme)
    )

    position = CablePosition(
        design,
        0.0,
        -0.20,
        Dict(component.id => (index == 1 ? 1 : 0)
        for (index, component) in enumerate(design.components))
    )
    system = LineCableSystem("manual-system", 1000.0, position)
    earth = EarthModel(frequency, 100.0, 10.0, 1.0)
    push!(
        gallery,
        "System preview" => preview(
            system;
            earth_model = earth,
            backend,
            display_plot,
            export_theme
        )
    )
    push!(
        gallery,
        "Material scale" => LineCableModels.DataModel.show_material_scale(
            ; backend, display_plot, export_theme
        )
    )

    @assert length(gallery) == 23
    @assert all(pair -> pair.second isa UIPlot, gallery)
    @assert all(pair -> pair.second.context.backend === backend, gallery)
    @assert all(pair -> pair.second.page.export_spec.theme === export_theme, gallery)
    return gallery
end
