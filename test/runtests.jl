using Test
using Aqua
using LineCableModels

@testset "LineCableModels" begin
    @testset "Aqua" begin
        Aqua.test_all(LineCableModels)
    end
    include("commons/commons.jl")
    include("materials/materials.jl")
    include("earth/earth.jl")
    include("datamodel/geometry.jl")
    include("datamodel/design.jl")
    include("datamodel/system.jl")
    include("architecture/architecture.jl")
end
