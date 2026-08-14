using LineCableModels
using CairoMakie
using Measurements: measurement

const REFERENCE_DIRECTORY = @__DIR__

function save_reference(name, plot_handle)
    CairoMakie.save(joinpath(REFERENCE_DIRECTORY, "$name.png"), plot_handle.figure)
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
    first(plot(parameters; mode = :RLCG, backend = :cairo, display_plot = false))
)
save_reference(
    "line_zy_cartesian",
    first(plot(
        parameters; mode = :ZY, coord = :cart, backend = :cairo, display_plot = false))
)
save_reference(
    "line_zy_polar",
    first(plot(
        parameters; mode = :ZY, coord = :polar, backend = :cairo, display_plot = false))
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
    first(plot(measurement_parameters; mode = :RLCG, backend = :cairo, display_plot = false))
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
        measurement(2.5, 0.0),
        measurement(2.5, 0.0),
        measurement(2.5, 0.0)
    ),
    4,
    0.95
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
load!(library; file_name = joinpath(pkgdir(LineCableModels), "test", "cable_test.json"))
design = first(values(library.data))
save_reference(
    "cable_preview",
    preview(design; backend = :cairo, display_plot = false)
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
earth = EarthModel(frequency, 100.0, 10.0, 1.0)
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
    show_material_scale(backend = :cairo, display_plot = false)
)
