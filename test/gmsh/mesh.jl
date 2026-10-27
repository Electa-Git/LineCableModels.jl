@testset "Gmsh mesh" begin
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    options = ComputationOptions(GmshBackend; workers = 2)
    @testset "$case" for (case, reference) in REFERENCES[1:2]
        problem = LineParametersProblem(reference.system, reference.earth;
            temperature = reference.temperature, frequencies = reference.frequencies[1:1])
        model = extension.FieldModel(problem, FEM())
        mktempdir() do directory
            mkpath(joinpath(directory, "input"))
            cp(joinpath(@__DIR__, "..", "..", "ext", "LineCableModelsGmshExt", "getdp"),
                joinpath(directory, "input", "getdp"))
            write(joinpath(directory, "input", "model_data.pro"), model, options)
            gmsh.initialize()
            try
                geometry = build(extension.Geometry, model, options)
                write(joinpath(directory, "input", "physical.geo"), geometry)
                paths = extension.mesh(directory, model, options)
                path = joinpath(directory, "mesh", "frequency-0001.msh")
                @test paths == [path]
                @test isfile(path)
                @test extension.validate(path, model) == path
                gmsh.model.add("mesh-counts")
                gmsh.merge(path)
                @test (1, 9001) in gmsh.model.getPhysicalGroups(1)
                @test (0, 9002) in gmsh.model.getPhysicalGroups(0)
                gmsh.model.remove()
            finally
                gmsh.finalize()
            end
        end
    end
end
