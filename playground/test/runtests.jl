using Dates
using LineCableModelsPlayground
using LineCableModelsPlaygroundProtocol
using Test
include("published_assets.jl")
using TOML

include("architecture.jl")
include("application_catalogue.jl")
include("artifacts.jl")
include("cli.jl")
include("container_runtime.jl")
include("geographic_map.jl")
include("power_system_canvas.jl")
include("presentation.jl")
include("jobhandle.jl")
include("repeater.jl")
include("ribbon.jl")
include("toolkit.jl")
include("runtime_controls.jl")
include("julia_terminal.jl")
include("scientific_views.jl")
include("runtime_conformance.jl")
include("uploads.jl")
include("visual_contracts.jl")
include("workbench.jl")
include("xray_preview.jl")
include("ui_precompile.jl")

if haskey(ENV, "NATS_TEST_PUBLISHER_URL")
    include("broker_lifecycle.jl")
end
