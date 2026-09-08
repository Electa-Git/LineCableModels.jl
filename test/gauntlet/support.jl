@testmodule GauntletSupport begin
    using Gmsh
    using Measurements
    include(joinpath(@__DIR__, "runtime.jl"))
    include(joinpath(@__DIR__, "campaign.jl"))
end
