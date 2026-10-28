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
            @testset "GetDP" begin
                primitives = extension.scan(directory, model, options)
                n = length(model.label)
                @test size(primitives.Z) == (n, n, 1)
                @test size(primitives.P) == (n, n, 1)
                @test all(isfinite, primitives.Z)
                @test all(isfinite, primitives.P)
                for suffix in (".done", "-Z.tsv", "-P.tsv")
                    @test all(
                        b -> isfile(joinpath(directory, "results", "jobs",
                            "getdp-f0001-b" * lpad(b, 4, '0') * suffix)),
                        1:n)
                end
                header = "frequency_index\tfrequency_hz\tresponse_terminal\tbasis_terminal\treal\timaginary"
                for name in ("Z", "P")
                    lines = readlines(joinpath(directory, "results", name * ".tsv"))
                    @test first(lines) == header
                    @test length(lines) == n^2 + 1
                end
                @test isfile(joinpath(directory, "logs", "frequency-0001.log"))
                @test isempty(readdir(joinpath(directory, "work")))
                executable = realpath(extension.GetdpArtifact())
                @test isfile(executable)
                buffer = IOBuffer()
                run(pipeline(Cmd([executable, "-info"]); stdout = buffer, stderr = buffer))
                output = String(take!(buffer))
                @test occursin("getdp", lowercase(output))
                @test occursin(r"Version\s*:\s*3\.5(?:\.|\s|$)", output)
            end
        end
    end
    @testset "GetDP executable" begin
        reference = last(first(REFERENCES))
        problem = LineParametersProblem(reference.system, reference.earth;
            temperature = reference.temperature, frequencies = reference.frequencies[1:1])
        model = extension.FieldModel(problem, FEM())
        mktempdir() do directory
            absent = ComputationOptions(GmshBackend; getdp_executable = "/none")
            @test_throws Base.IOError extension.scan(directory, model, absent)
            executable = realpath(extension.GetdpArtifact())
            folder = joinpath(directory, repeat("a", 120), repeat("b", 120))
            mkpath(folder)
            path = joinpath(folder, "getdp")
            cp(executable, path)
            oversized = ComputationOptions(GmshBackend; getdp_executable = path)
            @test_throws "getdp_executable" extension.scan(directory, model, oversized)
        end
    end
end
