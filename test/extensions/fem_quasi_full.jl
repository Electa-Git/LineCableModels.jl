@testitem "Manual quasi-full / edge-element voltage paths" tags=[:extension] begin
    using LineCableModels, Gmsh
    include(joinpath(pkgdir(LineCableModels), "dev", "quasi_full_paths.jl"))
    triangles = [((0.,0.), (1.,0.), (1.,1.)), ((0.,0.), (1.,1.), (0.,1.))]
    points = quasi_full_path_points(triangles, [(0.2,0.), (0.8,1.)])
    @test sum(p[3] for p in points) ≈ 0.6
    @test sum(p[4] for p in points) ≈ 1.0
    # A lowest-order edge field is a constant vector plus c*(-y,x).
    @test sum(-p[2]*p[3] + p[1]*p[4] for p in points) ≈ 0.2
    # A shared mesh edge must be integrated once, in either orientation.
    edge = quasi_full_path_points(triangles, [(0.,0.), (1.,1.)])
    @test sum(p[3] for p in edge) ≈ 1.0
    reverse = quasi_full_path_points(triangles, [(0.8,1.), (0.2,0.)])
    @test sum(-p[2]*p[3] + p[1]*p[4] for p in reverse) ≈ -0.2
    loop = quasi_full_path_points(triangles, [(0.,0.), (1.,0.), (1.,1.), (0.,1.), (0.,0.)])
    @test sum(-p[2]*p[3] + p[1]*p[4] for p in loop) ≈ 2.0
    @test_throws ErrorException quasi_full_path_points(triangles, [(0.,0.), (0.,2.)])
end

@testitem "Manual quasi-full / native gauge and current-basis invariance" tags=[:extension, :integration, :fem_numerical] begin
    using LineCableModels, Gmsh, LinearAlgebra, Printf
    include(joinpath(pkgdir(LineCableModels), "dev", "quasi_full_paths.jl"))
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    wire = build(CableDesign, "bare", terminal(:core,
        core(Material(kind=:conductor, rho=1e-12); r=0.0425)))
    positions = [(0.,-1.), (1.,-1.)]
    system = build(LineCableSystem, [wire,wire], positions;
        connections=[Dict(:core=>1),Dict(:core=>2)], line_length=1.)
    problem = LineParametersProblem(system; frequencies=[1e5],
        earth_props=homogeneous(rho=0.1,eps_r=1.,mu_r=1.))
    fem = Formulation(:LineCableModelsFEM;
        options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false),
        fem_options=(gmsh_verbosity=0,getdp_verbosity=3,mesh_policy=:remesh))
    model = FEM._resolved_fem_model(FEM._preflight_fem_problem(problem), fem)
    mktempdir() do root
        info = FEM._create_run(root)
        path = joinpath(info.path, "paths.pro")
        mesh, renumbered = lock(FEM.FEM_SESSION_LOCK) do
            session = FEM._start_gmsh(0)
            try
                geometry = FEM._build_geometry!(model, "quasi-full-test")
                mesh = only(FEM._select_meshes!(info, model, geometry, fem, root))
                FEM._prepare_run_inputs!(info, model)
                write_quasi_full_paths(path, mesh, only(model.mesh_plans), model, positions, 0.0425)
                # Same elements and coordinates; a different node ordering
                # gives GetDP a different admissible spanning tree.
                tags, _, _ = gmsh.model.mesh.get_nodes()
                gmsh.model.mesh.renumber_nodes(tags, reverse(tags))
                renumbered = joinpath(info.path, "mesh", "renumbered.msh")
                gmsh.write(renumbered)
                # GetDP enumerates edges from the elements. Reverse the MSH 4.1
                # element blocks and their numbering as well as the nodes.
                lines = readlines(renumbered)
                first = findfirst(==("\$Elements"), lines)
                last = findfirst(==("\$EndElements"), lines)
                blocks = Tuple{String,Vector{String}}[]
                cursor = first + 2
                while cursor < last
                    count = parse(Int, split(lines[cursor])[4])
                    push!(blocks, (lines[cursor], lines[cursor+1:cursor+count]))
                    cursor += count + 1
                end
                open(renumbered, "w") do io
                    foreach(line -> println(io,line), lines[1:first+1])
                    element = 0
                    for (header, rows) in reverse(blocks)
                        println(io, header)
                        for row in reverse(rows)
                            element += 1
                            println(io, element, " ", split(row;limit=2)[2])
                        end
                    end
                    foreach(line -> println(io,line), lines[last:end])
                end
                (mesh, renumbered)
            finally
                FEM._finish_gmsh(session)
            end
        end
        executable = FEM._getdp_selection(fem).path
        source = joinpath(pkgdir(LineCableModels), "ext", "LineCableModelsGmshExt", "getdp", "quasi-full.pro")
        function solve_case(name, mesh; reuse=true, pec=true)
            directory = joinpath(root, name); mkpath(directory)
            command = FEM._getdp_command(executable, source, mesh, info, fem,
                only(model.mesh_plans), [1,2], directory; reuse_factorization=reuse)
            open(joinpath(directory, "getdp.log"), "w") do io
                run(pipeline(`$command -setstring PathDataPath $path -setnumber PerfectConductors $(Int(pec))`; stdout=io, stderr=io))
            end
            matrices = map(("Z", "P", "Pscalar")) do quantity
                matrix = zeros(ComplexF64, 2, 2)
                for basis in 1:2
                    raw = joinpath(directory, "raw/jobs", @sprintf("getdp-f0001-b%04d-%s.tsv", basis, quantity))
                    @test FEM._valid_job_raw(raw, 2, 1, 1e5, basis)
                    for line in eachline(raw)
                        row = split(line)
                        matrix[parse(Int,row[3]),basis] = complex(parse(Float64,row[5]),parse(Float64,row[6]))
                    end
                end
                matrix
            end
            matrices
        end
        z, p, scalar = solve_case("original", mesh)
        z2, p2, scalar2 = solve_case("renumbered", renumbered)
        @test norm(scalar2-scalar) > 0.001norm(scalar)
        @test z2 ≈ z rtol=1e-7
        @test p2 ≈ p rtol=1e-7
        z3, p3, _ = solve_case("fresh-factors", mesh; reuse=false)
        @test z3 ≈ z rtol=1e-8
        @test p3 ≈ p rtol=1e-8
        @test minimum(eigvals(Symmetric(real.(inv(p))))) > 0
        # The retained finite-metal branch must remain executable as well.
        z4, p4, _ = solve_case("finite-metal", mesh; pec=false)
        @test z4 ≈ z rtol=1e-4
        @test p4 ≈ p rtol=1e-6
    end
end
