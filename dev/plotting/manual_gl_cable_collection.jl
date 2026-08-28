using Test
using LineCableModels
using LineCableModels.DataModel
using GLMakie

const CABLE_COLLECTION_GALLERY_SMOKE_ONLY = lowercase(
    get(ENV, "LINECABLEMODELS_GL_GALLERY_SMOKE", "false")
) == "true"

set_backend!(:gl)

# Load one detailed design and derive renderer-independent variants. The
# collection recipe accepts an ordinary vector; no gallery or Makie wrapper
# type is required at the call site.
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
) for cable_id in (
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
# feature grammar is used; only the resolved grid differs.
explicit = preview(
    designs[1:4];
    layout = (1, 4),
    size = (1600, 500),
    backend = :gl,
    display_plot,
    open_export = false
)

@testset "manual GL cable collection preview" begin
    @test automatic.page.payload.layout == (2, 3)
    @test explicit.page.payload.layout == (1, 4)
    @test length(automatic.panels) == 5
    @test length(explicit.panels) == 4
    @test automatic.context.legend === nothing
    @test explicit.context.legend === nothing
    @test length(automatic.context.colorbars) == 3
    automatic_valign = automatic.context.shell.side.valign[]
    explicit_valign = explicit.context.shell.side.valign[]
    @test automatic_valign == convert(typeof(automatic_valign), :top)
    @test explicit_valign == convert(typeof(explicit_valign), :top)
    @test [panel.axis.title[] for panel in automatic.panels] ==
          getproperty.(designs, :cable_id)
end

println("Built automatic 2×3 and explicit 1×4 cable-preview canvases.")
println("Confirm that every subplot title is its cable id and no legend is present.")
println("Confirm that each canvas has one shared, top-aligned set of material colorbars.")

if CABLE_COLLECTION_GALLERY_SMOKE_ONLY
    @test automatic.context.window === nothing
    @test explicit.context.window === nothing
    println("GL gallery smoke-only mode complete; native windows were not opened.")
    exit()
end

println("Close both windows to finish, or press Ctrl+C here.")

try
    while isopen(automatic.context.window) || isopen(explicit.context.window)
        sleep(0.1)
    end
finally
    GLMakie.closeall()
end
