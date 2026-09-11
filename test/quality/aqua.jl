@testitem "Quality / Aqua / package hygiene" tags = [:quality] begin
    using Aqua
    # Extensions share the root environment. Pkg is intentionally
    # imported only by LineCableModelsGmshExt, so it is invisible to Aqua's
    # root-module stale-dependency inspection. CairoMakie is likewise loaded
    # only on SVG export; dev/plotting/manual_gl.jl exercises its first-use load.
    Aqua.test_all(LineCableModels; stale_deps=(ignore=[:Pkg, :CairoMakie],))
end
