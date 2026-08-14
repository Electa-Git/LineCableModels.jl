using LineCableModels
using GLMakie
using CairoMakie

const GL_GALLERY_SMOKE_ONLY = lowercase(
    get(ENV, "LINECABLEMODELS_GL_GALLERY_SMOKE", "false")
) == "true"
const GL_GALLERY_ARTIFACT_DIRECTORY = abspath(get(
    ENV,
    "LINECABLEMODELS_GL_ARTIFACTS",
    joinpath(tempdir(), "linecablemodels-gl-artifacts")
))

mkpath(GL_GALLERY_ARTIFACT_DIRECTORY)
cd(GL_GALLERY_ARTIFACT_DIRECTORY)
set_backend!(:gl)

include(joinpath(@__DIR__, "_gallery_fixtures.jl"))
gallery = build_manual_plot_gallery(:gl; display_plot = !GL_GALLERY_SMOKE_ONLY)

println("Built $(length(gallery)) GLMakie inspection panels.")
println("SVG exports are written to $GL_GALLERY_ARTIFACT_DIRECTORY")

if GL_GALLERY_SMOKE_ONLY
    @assert all(pair -> pair.second.context.window === nothing, gallery)
    println("GL gallery smoke-only mode complete; native windows were not opened.")
    exit()
end

@assert all(pair -> pair.second.context.window !== nothing, gallery)
println("Each page is open in its own native GLMakie window.")
println("Close all windows to finish, or press Ctrl+C here.")

try
    while any(pair -> isopen(pair.second.context.window), gallery)
        sleep(0.1)
    end
finally
    GLMakie.closeall()
end
