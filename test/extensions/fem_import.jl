@testitem "FEM files / detached mesh and fields preserve native identities and caller state" tags=[:extension] begin
    using Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
    gmsh = Gmsh.gmsh
    fixtures = joinpath(pkgdir(LineCableModels), "test", "fixtures", "data", "fem")
    mesh_path = joinpath(fixtures, "sparse.msh")
    field_path = joinpath(fixtures, "discontinuous.pos")
    withenv("LINECABLEMODELS_GETDP" => "/unavailable/getdp", "DISPLAY" => "") do
        mesh = import_data(:msh, mesh_path)
        @test !Bool(gmsh.is_initialized())
        @test mesh isa Base.get_extension(LineCableModels, :LineCableModelsGmshExt).FEMMesh
        @test mesh.node_tags == [7, 23, 64, 105]
        @test mesh.coordinates == [0 1 0 1; 0 0 1 1; 0 0 0 0]
        block = only(mesh.blocks)
        @test block.element_tags == [5, 19]
        @test block.connectivity == [1 2; 2 4; 3 3]
        @test block.physical_tags == [8]
        @test mesh.physical_names[(2, 8)] == "sample material"
        @test import_data(:msh, mesh_path; coordinate_scale = 1e-3).coordinates ==
              mesh.coordinates .* 1e-3
        @test_throws ArgumentError import_data(:msh, mesh_path; coordinate_scale = 0)

        views = import_data(:pos, field_path)
        @test !Bool(gmsh.is_initialized())
        @test length(views) == 2
        scalar, vector = views
        @test scalar.representation === :complex
        @test vector.representation === :real # Two steps alone imply no phasor.
        @test length(scalar.times) == 2
        @test vector.times == [0.25, 0.75]
        @test size(only(scalar.blocks).values) == (1, 3, 2, 2)
        @test only(scalar.blocks).coordinates[:, :, 1] == [0 1 0; 0 0 1; 0 0 0]
        @test only(scalar.blocks).values[1, :, 1, 1] == [1, 2, 3]
        @test only(scalar.blocks).values[1, :, 1, 2] == [9, 10, 11]
        @test only(scalar.blocks).values[1, :, 2, 2] == [12, 13, 14]
        @test only(vector.blocks).values[:, 1, 1, 1] == [3, 4, 0]
        @test import_data(:pos, field_path; view = 2, representation = :complex).representation ===
              :complex
        @test_throws ArgumentError import_data(:pos, field_path; view = 3)
        @test_throws ArgumentError import_data(:pos, field_path; view = true)
        @test_throws ArgumentError import_data(:pos, field_path; representation = :bad)
        @test_throws ArgumentError import_data(:msh, field_path)
        @test_throws ArgumentError import_data(:pos, mesh_path)
        @test !Bool(gmsh.is_initialized())

        # A native binary MSH 4.1 file carries multiple physical memberships.
        # The sparse ASCII fixture above provides the independent connectivity.
        gmsh.initialize(String[], false, false)
        try
            gmsh.option.set_number("General.Terminal", 0)
            gmsh.merge(mesh_path)
            gmsh.model.add_physical_group(2, [11], 9)
            gmsh.model.set_physical_name(2, 9, "second membership")
            gmsh.option.set_number("Mesh.MshFileVersion", 4.1)
            gmsh.option.set_number("Mesh.Binary", 1)
            mktempdir() do root
                path=joinpath(root, "binary.msh")
                gmsh.write(path)
                imported=import_data(:msh, path)
                @test imported.node_tags == mesh.node_tags
                @test imported.coordinates == mesh.coordinates
                @test only(imported.blocks).connectivity == only(mesh.blocks).connectivity
                @test Set(only(imported.blocks).physical_tags) == Set([8, 9])
                @test imported.physical_names[(2, 9)] == "second membership"
            end
        finally
            Gmsh.finalize()
        end

        gmsh.initialize(String[], false, false)
        try
            gmsh.model.add("caller")
            caller_view = gmsh.view.add("caller field")
            gmsh.option.set_number("General.Verbosity", 3)
            gmsh.onelab.set_string("LineCableModels/FEM/caller", ["untouched"])
            for path in (mesh_path, field_path)
                import_data(endswith(path, "msh") ? :msh : :pos, path)
                @test gmsh.model.get_current() == "caller"
                @test Set(gmsh.model.list()) == Set(["", "caller"])
                @test gmsh.view.get_tags() == [caller_view]
                @test gmsh.option.get_number("General.Verbosity") == 3
                @test gmsh.onelab.get_string("LineCableModels/FEM/caller") == ["untouched"]
            end
            mktempdir() do root
                bad = joinpath(root, "broken.msh")
                write(bad, "\$MeshFormat\nnot a mesh\n")
                @test_throws Exception import_data(:msh, bad)
                @test gmsh.model.get_current() == "caller"
                @test gmsh.view.get_tags() == [caller_view]
            end
        finally
            Gmsh.finalize()
        end
        @test only(scalar.blocks).values[1, 1, 1, 1] == 1
    end
end
