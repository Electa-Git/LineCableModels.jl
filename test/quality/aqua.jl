@testitem "Quality / Aqua / package hygiene" tags = [:aqua] begin
    using Aqua
    # Pkg is imported by LineCableModelsGmshExt, so the root-module
    # stale-dependency inspection cannot observe its use there.
    Aqua.test_all(LineCableModels; stale_deps=(ignore=[:Pkg],))
end
