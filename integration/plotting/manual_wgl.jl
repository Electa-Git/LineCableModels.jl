const WGL_MANUAL_ENABLED = lowercase(get(ENV, "LINECABLEMODELS_WGL_MANUAL", "false")) ==
                           "true"

if !WGL_MANUAL_ENABLED
    println("WGL manual gallery is disabled.")
    println("Run with LINECABLEMODELS_WGL_MANUAL=true to opt in.")
    exit()
end

using LineCableModels
using WGLMakie
using CairoMakie
using Measurements: measurement

const Bonito = WGLMakie.Bonito
const DOM = Bonito.DOM
const WGL_SMOKE_ONLY = lowercase(get(ENV, "LINECABLEMODELS_WGL_SMOKE", "false")) == "true"
const WGL_ARTIFACT_DIRECTORY = abspath(get(
    ENV,
    "LINECABLEMODELS_WGL_ARTIFACTS",
    "/tmp/linecablemodels-wgl-artifacts"
))

mkpath(WGL_ARTIFACT_DIRECTORY)
cd(WGL_ARTIFACT_DIRECTORY)
set_backend!(:wgl)

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
    Makie.plot(parameters; mode = :RLCG, backend = :wgl, display_plot = false)
)
add_pages!(
    "Line parameters: Z/Y Cartesian",
    Makie.plot(
        parameters;
        mode = :ZY,
        coord = :cart,
        backend = :wgl,
        display_plot = false
    )
)
add_pages!(
    "Line parameters: Z/Y polar",
    Makie.plot(
        parameters;
        mode = :ZY,
        coord = :polar,
        backend = :wgl,
        display_plot = false
    )
)
add_pages!(
    "Line parameters: measurement error bars",
    Makie.plot(
        measurement_parameters;
        mode = :RLCG,
        backend = :wgl,
        display_plot = false
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
            backend = :wgl,
            display_plot = false
        )
    )
end

library = CablesLibrary()
load!(library; file_name = joinpath(pkgdir(LineCableModels), "test", "cable_test.json"))
design = first(values(library.data))
push!(
    gallery,
    "Cable preview" => preview(design; backend = :wgl, display_plot = false)
)

position = CablePosition(
    design,
    0.0,
    -0.20,
    Dict(component.id => (index == 1 ? 1 : 0)
    for (index, component) in enumerate(design.components))
)
system = LineCableSystem("wgl-manual-system", 1000.0, position)
earth = EarthModel(frequency, 100.0, 10.0, 1.0)
push!(
    gallery,
    "System preview" => preview(
        system;
        earth_model = earth,
        backend = :wgl,
        display_plot = false
    )
)
push!(
    gallery,
    "Material scale" => show_material_scale(backend = :wgl, display_plot = false)
)

@assert length(gallery) == 23
@assert all(pair -> pair.second isa UIPlot, gallery)
@assert all(pair -> pair.second.context.backend === :wgl, gallery)

println("Built $(length(gallery)) WGL inspection panels.")
println("SVG exports are written to $WGL_ARTIFACT_DIRECTORY")

if WGL_SMOKE_ONLY
    println("WGL smoke-only mode complete; browser server was not started.")
    exit()
end

function gallery_card(title, handle)
    return DOM.section(
        DOM.h2(title; style = "margin: 0 0 0.5rem 0; font: 600 1.1rem sans-serif;"),
        handle.figure;
        style = "background: #e5e5e5; padding: 0.75rem; border-radius: 0.4rem;"
    )
end

app = Bonito.App() do
    cards = [gallery_card(title, handle) for (title, handle) in gallery]
    return DOM.main(
        DOM.h1("LineCableModels WGL manual gallery"),
        DOM.p(
            "Inspect reset, SVG export, logarithmic axes, legends, visibility, " *
            "zoom/pan, previews, material toggling, histograms, PDFs, ECDFs, and Q-Q plots."
        ),
        DOM.div(
            cards...;
            style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(820px, 1fr)); gap: 1rem;"
        );
        style = "max-width: 1800px; margin: auto; padding: 1rem; font-family: sans-serif;"
    )
end

host = get(ENV, "LINECABLEMODELS_WGL_HOST", "127.0.0.1")
port = parse(Int, get(ENV, "LINECABLEMODELS_WGL_PORT", "8082"))
server = Bonito.Server(app, host, port)
Bonito.HTTPServer.start(server)
println("Open http://$host:$port in a browser; press Ctrl+C here to stop.")

try
    wait(server)
finally
    close(server)
end
