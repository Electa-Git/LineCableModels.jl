@testmodule GauntletSupport begin
    using LineCableModels
    using Gmsh
    using Measurements
    include(joinpath(pkgdir(LineCableModels), "gauntlet", "Gauntlet.jl"))
    using .Gauntlet
end
