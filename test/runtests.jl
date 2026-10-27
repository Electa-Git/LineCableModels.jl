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
    include("datamodel/blueprint.jl")
    include("datamodel/builders.jl")
    include("engine/problem.jl")
    include("engine/laws.jl")
    include("engine/impedance.jl")
    include("engine/input.jl")
    include("references/problems.jl")
    include("engine/stages.jl")
    include("engine/reductions.jl")
    include("engine/integral.jl")
    include("engine/unified.jl")
    include("engine/formulation.jl")
    include("engine/compute.jl")
    include("engine/result.jl")
    include("gmsh/options.jl")
    include("gmsh/model.jl")
    include("architecture/architecture.jl")
end
