@testitem "Quality / Aqua / package hygiene" tags = [:quality] begin
    using Aqua
    # Extensions share the root environment. LazyArtifacts is intentionally
    # imported only by LineCableModelsGmshExt, so it is invisible to Aqua's
    # root-module stale-dependency inspection.
    Aqua.test_all(LineCableModels; stale_deps=(ignore=[:LazyArtifacts],))
end
