using Test
using LineCableModels
using LineCableModels.DataModel
using GLMakie

const CABLE_COLLECTION_GALLERY_SMOKE_ONLY = lowercase(
    get(ENV, "LINECABLEMODELS_GL_GALLERY_SMOKE", "false")
) == "true"

# Current construction inputs; no serialized historical model is loaded.
include(joinpath(@__DIR__,"../../../test/support/scenarios.jl"))
designs=[CurrentScenarios.coaxial_design(;scale=1+index/10,name="Cable $index") for index in 1:5]
display_plot = !CABLE_COLLECTION_GALLERY_SMOKE_ONLY

# Omitting `layout` exercises the near-square rule: five designs become a 2×3
# canvas. The one set of material colorbars is calculated from all five designs.
automatic = preview(
    designs;
    size = (1400, 900),
    backend = :gl,
    display_plot,
    open_export = false
)

# Passing `(rows, columns)` exercises the caller-owned layout choice. The same
# preview data is used; only the native grid differs.
explicit = preview(
    designs[1:4];
    layout = (1, 4),
    size = (1600, 500),
    backend = :gl,
    display_plot,
    open_export = false
)

@testset "manual GL cable collection preview" begin
    @test length(automatic.axes) == 5
    @test length(explicit.axes) == 4
    @test automatic.legend === nothing
    @test explicit.legend === nothing
    @test length(automatic.colorbars) == 3
    @test [axis.title[] for axis in automatic.axes] ==
          getproperty.(designs, :cable_id)
end

println("Built automatic 2×3 and explicit 1×4 cable-preview canvases.")
println("Confirm that every subplot title is its cable id and no legend is present.")
println("Confirm that each canvas has one shared, top-aligned set of material colorbars.")

if CABLE_COLLECTION_GALLERY_SMOKE_ONLY
    println("GL gallery smoke-only mode complete; native windows were not opened.")
    exit()
end

println("Inspect the native GLMakie figures, then press Enter to finish.")
readline()
GLMakie.closeall()
