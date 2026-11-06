using Test
using Aqua
using LineCableModels
using Measurements

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
    include("basegeometry/basegeometry.jl")
    include("engine/stages.jl")
    include("engine/reductions.jl")
    include("engine/integral.jl")
    include("engine/martins2026.jl")
    include("engine/ideal.jl")
    include("engine/gary1976.jl")
    include("engine/wedepohl1973.jl")
    include("engine/saad1996.jl")
    include("engine/lucca1994.jl")
    include("engine/ametani2009.jl")
    include("engine/carson1926.jl")
    include("engine/pollaczek1926.jl")
    include("engine/papadopoulos2010.jl")
    include("engine/wise1934.jl")
    include("engine/wise1948.jl")
    include("engine/xue2018.jl")
    include("engine/formulation.jl")
    include("engine/compute.jl")
    include("engine/result.jl")
    include("measurements/extension.jl")
    include("measurements/propagation.jl")
    include("gmsh/options.jl")
    include("gmsh/model.jl")
    include("gmsh/geometry.jl")
    # The FEM workflow sets LINECABLEMODELS_FEM to run the tests that mesh and solve.
    if haskey(ENV, "LINECABLEMODELS_FEM")
        include("gmsh/mesh.jl")
        include("gmsh/compute.jl")
    end
    include("architecture/architecture.jl")
end
