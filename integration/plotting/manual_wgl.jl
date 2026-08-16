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

include(joinpath(@__DIR__, "_gallery_fixtures.jl"))
gallery = build_manual_plot_gallery(:wgl; display_plot = false)

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
            "Inspect reset, SVG export, line-parameter logarithmic axes, legends, " *
            "visibility, zoom/pan, previews, material toggling, histograms without " *
            "log controls, PDFs, ECDFs, Q-Q plots, and responsive legend restoration."
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
