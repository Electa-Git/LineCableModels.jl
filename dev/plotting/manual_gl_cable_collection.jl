using Test
using LineCableModels
using LineCableModels.DataModel
using GLMakie

const CABLE_COLLECTION_GALLERY_SMOKE_ONLY = lowercase(
    get(ENV, "LINECABLEMODELS_GL_GALLERY_SMOKE", "false")
) == "true"

# Load one detailed design and derive renderer-independent variants. `preview`
# accepts an ordinary vector; no gallery or Makie wrapper type is required at
# the call site.
library = CablesLibrary()
load!(library;
    file_name = joinpath(
        pkgdir(LineCableModels),
        "test",
        "fixtures",
        "data",
        "mv_cable_design.json"
    ))
source = only(values(library.data))
designs = [build(
               CableDesign,
               cable_id,
               source.root
           )
           for cable_id in (
    "Detailed cable A",
    "Detailed cable B",
    "Detailed cable C",
    "Detailed cable D",
    "Detailed cable E"
)]
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
